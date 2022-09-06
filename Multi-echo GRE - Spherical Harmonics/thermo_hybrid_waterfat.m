function [m,wts,Ac,c,A,phi] = thermo_hybrid_waterfat(imgs, algp, scanp, lib, mswitch,phiswitch,minit,phiinit)
%% Multi-echo MR Thermometry in the upper leg at 7T using near-harmonic 2D reconstruction for initialization
%% HIMM: Harmonic Initialized Model-based Multi-echo
%Mrmonic Initialized Model-based Multi-echo
% Function to solve for heating, field shifts, and select a baseline in the
% presence of fat. 
%
% Algorithm code for the iterative solver was developed by Megan Poorman 
% and William Grissom. The code presented in this version does not use the 
% regularization term anymore as it is not focused on applications using
% localized heating. Alternarively, the drift field (Ac) was initialized 
% using near-harmonic 2D reconstruction (proposed by Salomir et al.)
%
% Creators: Megan Poorman, William Grissom
% Updated by: Mathijs Kikken
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% Inputs:
%   imgs: Dynamic multi-echo images (Nx x Ny x Necho)
%   ----- algorithm params
%     algp.order = order;
%     algp.nmiter = number of m fit iters;
%     algp.nciter = number of c fit iters;
%     algp.stopthresh = fit tolerance;
%     algp.suppress_info = whether updates will be visualized on the command window
%
%   ----- scan parameters
%     scanp.dim = [rows cols];
%     scanp.b0 = b0;
%     scanp.tes = tes;
%     scanp.prc = prc;
% 
%   ----- baseline library
%     lib.Wlib = Wlib;
%     lib.Flib = Flib;
%     lib.dw0lib = dw0lib;
%     lib.R2starlib = R2starlib;
%     lib.fatmodel = fatmodel;
%
% Output: 
%   m:      frequency shit due to heating
%   wts:    baseline weight
%   c:      polynomial coefficients fit
%   A:      polynomial basis functions
%   Ac:     estimated drift field
%   phi:    Tx/Rx gain estimated


%% --- Initialization --- %%
%--- construct 2nd order finite difference roughness penalty object
% keyboard
if isfield('algp','beta')
  R = Robject(ones(scanp.dim(1),scanp.dim(2)),'order',2,'beta',algp.beta,'type_denom','matlab');
else
  R = [];
end

%--- Initialize vars for backward compatibility
if ~isfield(algp,'usequadprog')
   algp.usequadprog = 1;
end
if ~isfield(algp,'nlib')
    algp.nlib = size(lib.Wlib,3);
end

%--- Create polynomial matrix of given order (and degree)
if strcmp(algp.method,'poly')
    A = get_polynomial_function(algp.order,scanp.dim,2);
elseif strcmp(algp.method,'spher')
    A = get_spherical_function(algp.order, scanp.dim);
else
    message = {'The function with coefficients that are to be '
        'optimized was specified incorrection, '
        'please set algo.method to the correct string. '};
    errorStruct.message = sprintf('%s\n',message{:});
    errorStruct.identifier = 'MyFunction:parameterIncorrecctlyDefined';
    error(errorStruct)
end

%--- Initialize Ac, Ac_salomir, m and phi estimates
if isfield(lib,'cinit')
    c = lib.cinit;
else
    c = zeros(size(A,2),1);
end

if scanp.dim(1)==scanp.dim(2)
    Ac = reshape(A*c,[scanp.dim(1) scanp.dim(2)]);
else
    Ac = reshape(A*c,[scanp.dim(2) scanp.dim(1)]).';
end

if ~exist('minit','var')
    m = zeros(scanp.dim(1),scanp.dim(2));
else
    m = minit;
end

if ~exist('phiinit','var')
    phi = zeros(scanp.dim(1),scanp.dim(2));
else
    phi = phiinit;
end 

% Define Ac initialization: to be used in salomir approach
if isfield(lib,'Acinit')
    Ac_salomir = lib.Acinit;
else
    Ac_salomir = zeros(scanp.dim);
end

if length(scanp.dim) == 1 %assume square
    scanp.dim = [scanp.dim scanp.dim];
end


%% --- Start HIMM iterations --- %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Summary of the current algorithm code

% The first two iterations will be used to define the drift and temperature
% fields in regions with considerable fat content. The mask that defines
% this region is calculated using the function 'get_fatty_regions'.
% The third iteration will start with near-harmonic 2D reconstruction to
% determine the drift field in the entire field of view. This estimated
% field will be expressed in coefficients (c) using QR decomposition. These
% coefficients will be further updated in subsequent iterations.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normalization factor (median of heating dynamics * num echoes)
algp.norm = abs(median(imgs,'all') * length(scanp.tes));

% calculate a terrible cost that is guaranteed to get iterations going
wts_temp = zeros(size(lib.Wlib,3),1); wts_temp(1) = 1;
cost = cost_eval(imgs,scanp,lib,wts_temp,Ac,m,phi,algp); 
costOld = 2*cost;
ii = 0; % iteration counter

% The optimization algorithm is run for at least 4 iterations:
% 2 iterations without Salomir and 2 iterations with Salomir fitting
while costOld - cost > algp.stopthresh*costOld || ii < 6     
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % run f_update with current m, Ac, phi
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Build Z matrix
    for jj = 1:size(lib.Wlib,3)
        % Compute model values 
        Z(:,:,:,jj) = calcmeimgs(lib.Wlib(:,:,jj),lib.Flib(:,:,jj), ...           
                scanp.tes,scanp.b0,lib.dw0lib(:,:,jj)+Ac, ... 
                lib.R2starlib(:,:,jj),m,scanp.prc,lib.fatmodel);
    end
    Z = Z.*repmat(exp(1i*phi),[1 1 size(Z,3) size(Z,4)]); 
    
    % First iteration
    if ii == 0   % in first iteration its best to use magnitude images for the fit,
                 % in case there is a large phase shift between baselines
                 % and dynamic image
        if algp.nlib < size(lib.Wlib,3)
            % compress the library down to the nlib entries that best 
            % correlate with the dynamic image
            [wts,algp.libinds] = f_update(abs(imgs),abs(Z),algp.usequadprog,algp.nlib);
            lib.Wlib = lib.Wlib(:,:,algp.libinds);
            lib.Flib = lib.Flib(:,:,algp.libinds);
            lib.dw0lib = lib.dw0lib(:,:,algp.libinds);
            lib.R2starlib = lib.R2starlib(:,:,algp.libinds);
        else
            wts = f_update(imgs,Z,algp.usequadprog);
        end
        
        % Update normalization factor: median of heating dynamics * num echoes
        [~, body_mask, ~] = get_fatty_regions(lib.Wlib,lib.Flib,wts,scanp.dim,algp.th_fat,algp.no_signal,algp.SE,algp.body_factor,algp.body_erode);
        body_masks = repmat(body_mask,1,1,length(scanp.tes),size(imgs,4));
        imgsdyn_masked = imgs; imgsdyn_masked(body_masks ~= 1) = nan; 
        imgsdyn_masked = imgsdyn_masked(:);
        imgsdyn_masked = imgsdyn_masked(~isnan(imgsdyn_masked));
        algp.norm = abs(median(imgsdyn_masked) * length(scanp.tes));
    end
    
    % For the first few iterations, Ac and m were not used for updating of baseline 
    % coefficients; but after Salomir fitting, f_update is called again
    if ii > 2
        wts = f_update(imgs,Z,algp.usequadprog);
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % run c_update with current f, m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % Evaluate in voxels where the body is present, such that 
    % fitting is insensitive to background noise
    [fat_idx, body_mask, ~] = get_fatty_regions(lib.Wlib,lib.Flib,wts,scanp.dim,algp.th_fat,algp.no_signal,algp.SE,algp.body_factor,algp.body_erode);
    fat_idx = repmat(fat_idx,1,1,length(scanp.tes));
    body_masks = repmat(body_mask,1,1,length(scanp.tes));
    
    for kk = 1:algp.nciter
        
        % Calculate and apply weights to baseline images
        Ztot = zeros([scanp.dim(1) scanp.dim(2) length(scanp.tes)]);
        
        for jj = 1:size(lib.Wlib,3)
            % Compute model values
            Z = calcmeimgs(lib.Wlib(:,:,jj),lib.Flib(:,:,jj), ...
                scanp.tes,scanp.b0,lib.dw0lib(:,:,jj)+Ac, ...
                lib.R2starlib(:,:,jj),m,scanp.prc,lib.fatmodel);
            Ztot = Ztot + Z*wts(jj);
        end
        Ztot = Ztot.*repmat(exp(1i*phi),[1 1 size(Z,3)]);
        
        % Select voxels on which coefficient fitting will be determined:
        %   We evaluate only on fat when Salomir fitting has not yet occured
        %   After Salomir fitting, we are more confident on proper
        %   initialization and we dare to fit on the entire body region
        if ii > 2
            % Evaluate on pixels in the body
            % Voxels in the background will not be evaluated         
            Ztot(body_masks ~= 1) = 0;
            imgs_temp = imgs;
            imgs_temp(body_masks ~= 1) = 0;
        else
            % Evaluate on pixels in fatty regions
            Ztot(fat_idx ~= 1) = 0;
            imgs_temp = imgs;
            imgs_temp(fat_idx ~= 1) = 0; 
        end
        
        % Update c
        c = c_update(imgs_temp,Ztot,A,c,scanp.tes,scanp.prc);
        
        % Recalculate frequency shift map
        if scanp.dim(1)==scanp.dim(2)
            Ac = reshape(A*c,[scanp.dim(1) scanp.dim(2)]);
        else
            Ac = reshape(A*c,[scanp.dim(2) scanp.dim(1)]).';
        end
 
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % run Salomir with recently acquired Ac
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Get combined water and fat images to determine relative water/fat contents in each pixel    
    [Sal_mask, body_mask, signal_mask] = get_fatty_regions(lib.Wlib,lib.Flib,wts,scanp.dim,algp.th_fat,algp.no_signal,algp.SE,algp.body_factor,algp.body_erode);
    
    % The Salomir algorithm is only performed after 2 iterations
    if ii == 2
        % Describe the background field as a harmonic function 
        Ac = Salomir_fit(Ac,Ac_salomir,Sal_mask,algp.Nstd_ol,algp.neighbours,algp.neighbours);
        
        % Make sure the next Salomir fitting algorithm is initialized with the previous
        Ac_salomir = Ac;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % run c_fit to update c to approximate the Salomir map
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Create input parameters
        % indepvar: the independent variable is the physical mask with ones in
        %           regions where off-resonance is assumed to be estimated correctly
        % depvar:   the estimated off-resonace in the masked body region (dependent variable)
        indepvar = body_mask;
        depvar = Ac(indepvar == 1); 

        % Perform the fit with the desired order (and degree) as provided in
        % the polynomial / spherical harmomics array (A)
        % This fit makes use of QR decomposition to decompose into the
        % basic functions of polynomials / spherical harmonics
        c = c_fit(indepvar,depvar,scanp.dim,A);

        if scanp.dim(1)==scanp.dim(2)
            Ac_fit = reshape(A*c,[scanp.dim(1) scanp.dim(2)]);
        else
            Ac_fit = reshape(A*c,[scanp.dim(2) scanp.dim(1)]).';
        end
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % run phi_update with current f, m, Ac
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % This accounts for any universal DC offsets due to reciever gain etc.
    if phiswitch
        for kk = 1:algp.nciter
            
            % Calculate and apply weights to baseline images
            Ztot = zeros([scanp.dim(1) scanp.dim(2) length(scanp.tes)]);
 
            for jj = 1:size(lib.Wlib,3)
                % Compute model values
                Z = calcmeimgs(lib.Wlib(:,:,jj),lib.Flib(:,:,jj), ...
                    scanp.tes,scanp.b0,lib.dw0lib(:,:,jj)+Ac, ...
                    lib.R2starlib(:,:,jj),m,scanp.prc,lib.fatmodel);
                Ztot = Ztot + Z*wts(jj);
            end
            Ztot = Ztot.*repmat(exp(1i*phi),[1 1 size(Z,3) size(Z,4)]); 

            % Update phi
            phi = phi_update(imgs,Ztot,phi(:),scanp.tes,scanp.prc);

            % Recalculate frequency shift map
            phi = reshape(phi,[scanp.dim(1) scanp.dim(2)]);

        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % run m_update with current f, Ac
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if mswitch
        for kk = 1:algp.nmiter

            % Build Z matrix using true values of Ac and m
            Zw = zeros([scanp.dim(1) scanp.dim(2) length(scanp.tes)]);
            Zf = zeros([scanp.dim(1) scanp.dim(2) length(scanp.tes)]);
            
            for jj = 1:size(lib.Wlib,3)
                % Compute model values of water and fat separately
                [~,foow,foof] = calcmeimgs(lib.Wlib(:,:,jj).*exp(1i*phi),lib.Flib(:,:,jj).*exp(1i*phi), ...
                    scanp.tes,scanp.b0,lib.dw0lib(:,:,jj)+Ac, ... 
                    lib.R2starlib(:,:,jj),m,scanp.prc,lib.fatmodel);
                Zw = Zw + foow*wts(jj);
                Zf = Zf + foof*wts(jj);
            end

            % Update m
            m = m_update(imgs-Zf,Zw,m,scanp.tes,scanp.prc,algp,[],body_mask,signal_mask,Sal_mask,R);
            
        end
    end
        
    % report cost + iteration
    costOld = cost;ii = ii + 1;
    cost = cost_eval(imgs,scanp,lib,wts,Ac,m,phi,algp);
    if algp.suppress_info ~= 1
        fprintf('Optimization iteration %d: Cost = %0.2d | difference = %0.2d\n',ii,cost,cost-costOld);
    end
end        
    
%% Cost function
% Function to evaluate the overall cost function
function cost = cost_eval(imgs,scanp,lib,wts,Ac,m,phi,algp)
    Ztot = zeros([scanp.dim(1) scanp.dim(2) length(scanp.tes)]);
    for jj = 1:size(lib.Wlib,3)
        Z = calcmeimgs(lib.Wlib(:,:,jj),lib.Flib(:,:,jj), ...
            scanp.tes,scanp.b0,lib.dw0lib(:,:,jj)+Ac, ...
            lib.R2starlib(:,:,jj),m,scanp.prc,lib.fatmodel);
        Ztot = Ztot + Z*wts(jj);
    end

    Ztot = Ztot.*repmat(exp(1i*phi),[1 1 size(Z,3) size(Z,4)]);
    cost = norm((imgs(:)-Ztot(:))/algp.norm)^2;