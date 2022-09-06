%% Multi-echo MR Thermometry in the upper leg at 7T using near-harmonic 2D reconstruction for initialization
% 
%  HIMM: Harmonic Initialized Model-based Multi-echo  
%  Main script to extract temperature and drift fields from measured data
%
% Creator: Mathijs Kikken (University Medical Center Utrecht)
% The basis of this algorithm was created by Megan Poorman & William Grissom
% (https://github.com/poormanme/waterFatSeparated_MRThermometry)
% The code provided here does not use the l1 norm for regularization, but 
% rather initializes the drift field (spatial B0 distortion) using the
% algorithm of near-harmonic 2D reconstruction (proposed by Salomir et al.)
%
% Updated from: Megan Poorman, William Grissom (Vanderbilt University Institute of Imaging Science)
% Do not reproduce, distribute, or modify without proper citation according
% to license file

clear all; close all;

% Add the software (HIMM algorithm code and ISMRM toolbox) to the directory
addpath(genpath('.\Code\MRT optimization\HIMM MRT\Multi-echo GRE - Spherical Harmonics'))
addpath(genpath('.\Code\ISMRM Water Fat Toolbox'))

%% Define scanner and sequence parameters

b0 = 7;                     % Tesla
frequency_system = 298*1e6; % Frequency measured by used MR system [Hz]
prc = 1;                    % Direction of precession (scanner-dependent) 
                            % +1 = heating induces negative apparent freq shift
                            % -1 = heating induces positive apparent freq shift
phi0shift = 0;              % Tx/Rx gain 
alpha = -0.01;              % ppm/deg C;
gamma = 42.5778;            % MHz/T
rad2degC = 1/(2*pi*b0*alpha*gamma);  % Factor to convert radians to temperature


%% Load baseline (imgslib) data and dynamic (imgsdyn) data from which temperature must be extracted

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Add code to load the MRT images
% The following variables must be defined
% 1. scan          - 4D Matrix with image data, dimensions comprise: [nx,ny,TE,dyns]
% 2. imDataParams  - Image parameters used for the water/fat reconstruction
% 3. algoParams    - Algorithm parameters used for the water/fat reconstruction
%                    Check the ISMRM Fat/Water toolbox on how to construct
%                    these parameters
% imgslib is the dynamics which will be used in the library
% Usually I only have 1 image in the library: the first dynamic
% -> imgslib = scan(:,:,:,[1]);
% imgslib must subsequently be inserted in the struct 'imDataParams'

% 2. Use the ISMRM Fat-Water toolbox to separate into water, fat, R2*, and dw0 
% % Selection of fat model
% fatmodel = get_fatmodel(num_fat_peaks);
% 
% disp('Separating baseline WF - using ISMRM toolbox');
% 
% % WF separation is done for each dynamic individually with image size of [ nx | ny | 1 | ncoils | nTE ]
% for ii = 1:size(imgslib,4)
%     % Perform Mixed-Magnitude method initialized with the Graph Cut method
%     disp(['.. Executing baseline library dynamic ' num2str(ii) ' out of ' num2str(size(imgslib,4))]);
%     % Call the function 'fw_i2cm1i_3pluspoint_hernando_graphcut' from the fat-water toolbox
%     outParams = fw_i2cm1i_3pluspoint_hernando_graphcut(imDataParams,algoParams);
%     
%     % Combine data (water, fat, off-resonance and R2*) of individual dynamics
%     Wlib(:,:,ii) = outParams.species(1).amps;
%     Flib(:,:,ii) = outParams.species(2).amps;
%     dw0lib(:,:,ii) = outParams.fieldmap*2*pi;
%     R2starlib(:,:,ii) = outParams.r2starmap;
%     
% end

% 3. Load the dynamic data (imgsdyn), from which temperature must be extracted
% The following variables must be defined
% 1. imgsdyn     - 4D Matrix with image data, dimensions comprise: [nx,ny,TE,dyns]
% 2. tes_heat    - echo times (in seconds) corresponding to the third dimension in 'imgsdyn'
% imgslib is the dynamics which will be used in the library.
% Usually I only have 1 image in the library: the first dynamic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For now I created arteficial MRT data, to show the principles of the algorithm
% The parameters Wlib, Flib, R2starlib, and dw0lib should normally be
% generated from imgslib using the ISMRM Fat-Water toolbox. Here, the data is
% arteficially generated...
dim = 60;               % dimensionality of the arteficial data
fatPercentage = 0.5;    % percentage of fat present in regions that have both water and fat
hssig2 = 0.005;         % sigma for the hotspot
order = 2;              % polynomial order of the arteficially generated drift field
dw0shift = 2*pi;        % background polynomial frequency shift (will be replicated to all polynomial coeffs) (rad/sec)   
num_fat_peaks = 6;                          % Choose an x-peak model, options are:
                                            % [ 3 | 4 | 5 | 6 | 7 | 9 ]                                 
tes_heat = [0.0148 0.0155 0.0161 0.0168 0.0174 0.0181 0.0188 0.0194 0.0201 0.0207 0.0214]; % echo time in seconds
noiselevel = 0.05;      % add in noise (0.2 for SNR of 40)

% Generate the hotspot used for heating
[x,y] = meshgrid(-dim/2:dim/2-1);
hotspot = exp(-hssig2*(x.^2 + y.^2)); % normalized hot spot

maxtemps = linspace(0,1,25); % temperature increase (in degC) as the dynamic count increases
for dyn = 1:length(maxtemps)
    hotspot_data(:,:,dyn) = hotspot*alpha*maxtemps(dyn)*b0*2*pi*gamma;
end

% Generate the drift field
[yc,xc] = meshgrid(linspace(-1/2,1/2,dim));
yc = yc(:);
xc = xc(:);
A = [];
for yp = 0:order
    for xp = 0:(order-yp)
        A = [A (xc.^xp).*(yc.^yp)];
    end
end
coefficients = zeros(size(A,2),1) + dw0shift;  
drift = reshape(A*coefficients,[dim dim]);    % polynomial background phase shift

% Generate the true water and fat images
[x,y] = meshgrid(-dim/2:dim/2-1);
oval_large = x.^2 + 1.5*y.^2 <= (0.8*dim/2)^2; % larger oval
oval_small = x.^2 + 1.5*y.^2 <= (0.7*dim/2)^2; % smaller oval

% Fat part of the image is in the outer circle
% Water part the region within the outer circle, which constitutes fat
F = oval_large - oval_small;
W = oval_large;

% Apply WF ratio and initialize R2* and dw0
Wlib = W*(1-fatPercentage);
Flib = F*fatPercentage;
R2starlib = zeros(dim);
dw0lib = zeros(dim)+dw0shift;
fatmodel = get_fatmodel(num_fat_peaks);

% Generate the dynamic images
for dyn = 1:length(maxtemps)
    noise = noiselevel*randn(dim,dim,length(tes_heat));
    imgsdyn(:,:,:,dyn) = calcmeimgs(Wlib,Flib,tes_heat,b0,dw0lib+drift,...
        R2starlib,hotspot_data(:,:,dyn),prc,fatmodel)+noise;
end

% Provide a visualization
subplot(2,2,1)
imagesc(Wlib); colorbar; title('water'); axis('off')
subplot(2,2,2)
imagesc(Flib); colorbar; title('fat'); axis('off')
subplot(2,2,3)
imagesc(drift); colorbar; title('drift'); axis('off')
subplot(2,2,4)
imagesc(hotspot_data(:,:,end)*rad2degC); colorbar; title('temperature'); axis('off')


%% Masking
% The drift field will be initialized using near harmonic 2D reconstruction, 
% this algorithm requires a mask of fatty regions from which the harmonic
% reconstruction can be initialized.
% In addition to the fat mask, we also create masks for the body and for
% regions with negligable signal
% Three parameters are to be optimized for the fat, body and no-signal masks:
%   1. th_fat    - threshold in fat/water image to assign voxels to fat mask 
%   2. no_signal - threshold in respectively water and fat image to define low signal region
%   3. SE        - structure element to erode thick fatty regions
th_fat = 0.4;               % Higher value results in smaller 'fat' mask
no_signal = [0.05, 0.05]; 	% Higher value results in larger 'no signal' mask (first corresponds to water and second value to fat)
SE = strel('diamond', 0);   % Larger SE means thinner fat mask, but if too large we keep the original non-eroded mask
body_factor = 0.9;          % Larger values result in larger body mask
body_erode = 0;             % Integer value that determined how much the body mask is eroded

% Compute the fat, body and no signal mask with current parameters
[fat_mask,body_mask,signal_mask] = get_fatty_regions(Wlib(:,:,1), Flib(:,:,1), 1, [size(imgsdyn,1) size(imgsdyn,2)], th_fat, no_signal, SE, body_factor, body_erode);

% Provide a visualization
subplot(1,3,1); imagesc(fat_mask); axis image; set(gca,'XTick',[]); set(gca,'YTick',[]); title('fat mask');
subplot(1,3,2); imagesc(body_mask); axis image; set(gca,'XTick',[]); set(gca,'YTick',[]); title('body mask');
subplot(1,3,3); imagesc(signal_mask); axis image; set(gca,'XTick',[]); set(gca,'YTick',[]); title('no signal mask');


%% Algorithm parameters

% Only algorithm parameters may be altered, which does not hold for scanner 
% and sequence since those parameters are constant

% Parameter structure: algorithm parameters
algp.order = 2;                 % Polynomial / spherical background phase shift order 
algp.method = 'spher';           % 'poly' for polynomial and 'spher' for spherical harmonics
algp.nmiter = 5;                % # m iterations
algp.nciter = 5;                % # c iterations
algp.stopthresh = 1e-1;         % algorithm stopping threshold
algp.suppress_info = 0;         % Whether command window information is shown per iteration or not
initWithPreviousEstimate = 1;   % Parameter to speed up the algorithm

% Parameter structure: algorithm masking parameters
algp.th_fat = th_fat;           % threshold for fat definition for when l1 norm is only calculated in fats
algp.no_signal = no_signal; 	% threshold for percentage of voxel intensity in water and fat image to specify as voxel without signal
algp.SE = SE;                   % Stucturing Element used as final step to erode the thick fatty regions
algp.body_factor = body_factor; % Value used to vary the threshold for the in-phase image
algp.body_erode = body_erode;   % Value used to make the body mask slightly smaller
algp.Nstd_ol = 2;               % Filter out local extreme values in the fat border
algp.neighbours = 0;            % Define local neighbourhood
algp.min_nb = 0;                % Eliminate voxels at the fat-tissue interface from the fat mask (value between 0 and 8)

% Parameter structure: scan parameters
scanp.dim = [size(imgsdyn,1) size(imgsdyn,2)];
scanp.b0 = b0;
scanp.tes = tes_heat;
scanp.prc = prc;

% Parameter structure: baseline library
lib.Wlib = Wlib;
lib.Flib = Flib;
lib.dw0lib = dw0lib;
lib.R2starlib = R2starlib;
lib.fatmodel = fatmodel;        % x-peak fat model (defined prior to water/fat recon)

% Visualize the mask that is used for the salomir fit
mask = remove_outer_border_pixels(fat_mask,abs(imgsdyn(:,:,1,1)),algp.min_nb);
imagesc(mask); axis off; axis image;


%% Perform iterative optimization to separate the temperature and drift fields
clear lib cfunc mfunc wts phi;
disp('Optimization of multi-echo fat-suppressed MR thermometry has started..')

% Provide which dynamics are to be implemented in the optimization
show_ix = size(imgsdyn,4);   % number of dynamics that will be evaluated | size(imgsdyn,4)
init_ix = 1;                 % dynamic which will be evaluated first
step = 1;                    % step size in dynamics

% Updates onto the parameter structures
lib.Wlib = Wlib;
lib.Flib = Flib;
lib.dw0lib = dw0lib;
lib.R2starlib = R2starlib;
lib.fatmodel = fatmodel;
scanp.tes = tes_heat;               	% Update echo times to the ones used for the optimization part

tic
for ii = init_ix:step:init_ix+show_ix-1
    % Evaluate a single dynamic
    dyn_img = imgsdyn(:,:,:,ii);

    disp(['dyn ' num2str(ii)]);
    if initWithPreviousEstimate && ii == init_ix
        % Optimization of the first dynamic
        [mfunc(:,:,ii),wts(:,ii),Acfunc(:,:,ii),cfunc(:,ii),A,phi(:,:,ii)] = ...
            OptimizeHIMM(dyn_img, algp, scanp, lib,1,0,zeros(scanp.dim));
    elseif initWithPreviousEstimate && ii > init_ix
        % Optimization of remainder dynamics (with nonzero initialization of mfunc and drift)
        [mfunc(:,:,ii),wts(:,ii),Acfunc(:,:,ii),cfunc(:,ii),A,phi(:,:,ii)] = ...
            OptimizeHIMM(dyn_img, algp, scanp, lib,1,0,mfunc(:,:,ii-step));
    else
        % Phi switch off, m switch on (solve for temperature without Tx/Rx gain)
        [mfunc(:,:,ii),wts(:,ii),Acfunc(:,:,ii),cfunc(:,ii),A,phi(:,:,ii)] = ...
            OptimizeHIMM(dyn_img, algp, scanp, lib,1,0); 
    end
    lib.cinit = cfunc(:,ii);
    lib.Acinit = Acfunc(:,:,ii);
    
end
toc

%% Plot using step size defined in the optimization

figure; plot_ix = 1;
step_plot = floor((init_ix + show_ix - 1)/25);
if step_plot == 0
    step_plot = 1;
end
for ii = init_ix:step_plot:init_ix+show_ix-1
    subplot(5,5,plot_ix);
    imagesc(mfunc(:,:,ii).*body_mask*rad2degC, [-0.5 1]); axis image; axis off; colorbar; title(['dynamic ' num2str(ii)]);
    %imagesc(Acfunc(:,:,ii).*body_mask); axis image; axis off; colorbar; title(['dynamic ' num2str(ii)]);
    plot_ix = plot_ix + 1;
end
  

%% Lineplot of temperature over time in a specfic pixel
y_pixels = round([3*size(mfunc,1)/7, ...
                  3*size(mfunc,1)/4, ...
                  2*size(mfunc,1)/7, ...
                  3*size(mfunc,1)/4]);
x_pixels = round([1*size(mfunc,1)/4, ...
                  1*size(mfunc,1)/3, ...
                  4*size(mfunc,1)/7, ...
                  3*size(mfunc,1)/4]);
figure;
for ii = 1:length(x_pixels)
    subplot(2,length(x_pixels),ii); im_add = zeros(scanp.dim); im_add(y_pixels(ii),x_pixels(ii)) = 1;
    imagesc(mat2gray(abs(lib.Wlib(:,:,1)))+im_add); axis image;
    % Calculate average in 3x3 region that is to be plotted
    avg = 0;
    for xi = [y_pixels(ii)-1,y_pixels(ii),y_pixels(ii)+1]
        for yi = [x_pixels(ii)-1,x_pixels(ii),x_pixels(ii)+1]
            avg = avg + squeeze(mfunc(xi,yi,:))*rad2degC;
        end
    end
    subplot(2,length(x_pixels),ii+length(x_pixels));
    %plot(squeeze(mfunc(x_pixels(ii),y_pixels(ii),:))*rad2degC);
    to_plot = avg/9;
    plot(to_plot(:,:));
end