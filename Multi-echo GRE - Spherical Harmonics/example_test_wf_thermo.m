%% Multi-echo fat-suppressed MR thermometry using iterative separation 
%% of baseline water and fat images
%
%  Sample function to run on simulated data
%
% Creators: Megan Poorman, William Grissom
% Location: Vanderbilt University Institute of Imaging Science
% Created: 13/2025
% Updated: 05/2030
% Do not reproduce, distribute, or modify without proper citation according
% to license file

clear all; close all;

%% %%%%%%%%%%%%%%%%%%%%%%%%
% algorithm parameters
%%%%%%%%%%%%%%%%%%%%%%%%%
libsource = 'true'; % Switch to use true Water/Fat baseline images or WF separation output ('toolbox')
nciter = 5; % # c iterations
nmiter = 5; % # m iterations
niter = 10; % # outer iterations for script version
lam = 10^-5; % l1 regularization parameter for script version
mthresh = 0.001; % heat-induced frequency threshold for 'significant' heat
stopthresh = 0.0001; % algorithm stopping threshold
initWithPreviousEstimate = 1;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scanner+sequence parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b0 = 7; % Tesla
dim = 128; % image size
fatppm = 3.4; % mean fat offset in ppm, used to calculate echo times
fatfrq = b0*42.5778*fatppm; % fat frequency in Hz for echo time calc

% File location of the measurement that is used for simulations
baseline_test_path = 'Data/Processed_DIXON/M0009/baseline_test.mat';

% Provide which echo times are used
te1 = 14.1445*1e-3;        % imgslib: 1.56247*1e-3 | imgsdyn: 14.1445*1e-3 | test: 0.97*1e-3
dt = 2/fatfrq/3;        
num_echoes = 12;        % imgslib: 6 | imgsdyn: 10

tes = te1 + [0:dt:(num_echoes-1)*dt]; % echo time calculation
prc = 1;                % direction of precession (scanner-dependent) +1 = heating induces negative apparent freq shift; -1 = heating induces positive apparent freq shift
order = 3;              % polynomial background phase shift order 
hssig2 = 0.1;           % hot spot sigma^2
dw0shift = 2*pi*0.008;
phi0shift = 0;          % Tx/Rx gain 
alpha = -0.01;          % ppm/deg C;
dynmotind = 2;          % index of motion vector to use for dynamic image
gamma = 42.5778;        % MHz/T
snr = 40;               % add noise
fatPercentage = 0.05;    % 0.5 = 50% fat

disp(tes*1e3);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate baseline images 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

water_fat_image = 4;        % Options are:
                            % 1. Original water oval with 3 fat circles inside
                            % 2. Fat shell that surrounds a water oval (shell and oval overap near the edges)
                            % 3. Fat shell that surrounds a water oval + some fatty structures inside the water oval
                            % 4. Fat shell that surrounds a water oval + some fatty structures inside the water oval (regions of sat and bone are small)
                            % 5. Load water / fat / dw0 / R2star image from a .mat file
                            
%--- generate the true W,F images
if water_fat_image == 1
    
    disp('Generating true Water and Fat images');
    [x,y] = meshgrid(-dim/2:dim/2-1);
    WTrue = x.^2 + 1.5*y.^2 <= (0.8*dim/2)^2; % water oval

    f = (x.^2 + y.^2) <= (0.2*dim/2)^2; % small fat circle
    FTrue = circshift(f,[0 -dim/4]) + f + circshift(f,[0 dim/4]); % replicate fat circle
    R2starTrue = zeros(dim);
    dw0True_base = zeros(dim)+dw0shift;

    %apply WF ratio
    WTrue = WTrue*(1-fatPercentage);
    FTrue = FTrue*fatPercentage;
    
elseif water_fat_image == 2
    
    disp('Generating true Water and Fat images');
    [x,y] = meshgrid(-dim/2:dim/2-1);
    WTrue = x.^2 + 1.5*y.^2 <= (0.8*dim/2)^2; % water oval

    FTrue = 1.3*x.^2 + 2.0*y.^2 <= (0.8*dim/2)^2; % Slightly smaller sized fat oval
    FTrue = abs(FTrue - WTrue);
    R2starTrue = zeros(dim);
    dw0True_base = zeros(dim)+dw0shift;

    %apply WF ratio
    WTrue = WTrue*(1-fatPercentage);
    FTrue = FTrue*fatPercentage;
    
elseif water_fat_image == 3
    
    disp('Generating true Water and Fat images');
    [x,y] = meshgrid(-dim/2:dim/2-1);
    WTrue = x.^2 + 1.5*y.^2 <= (0.8*dim/2)^2; % water oval

    FTrue_part1 = 1.3*x.^2 + 2.0*y.^2 <= (0.8*dim/2)^2; % Slightly smaller sized fat oval: SAT
    FTrue_part1 = abs(FTrue_part1 - WTrue);
    FTrue_part2 = (x.^2 + y.^2) <= (0.2*dim/2)^2; % small bone region
    FTrue = FTrue_part1 + FTrue_part2;
    R2starTrue = zeros(dim);
    dw0True_base = zeros(dim)+dw0shift;

    %apply WF ratio
    WTrue = WTrue*(1-fatPercentage);
    FTrue = FTrue*fatPercentage;
    
elseif water_fat_image == 4
    
    disp('Generating true Water and Fat images');
    [x,y] = meshgrid(-dim/2:dim/2-1);
    WTrue = x.^2 + 1.5*y.^2 <= (0.8*dim/2)^2; % water oval

    FTrue_part1 = 1.25*x.^2 + 1.75*y.^2 <= (0.8*dim/2)^2; % Slightly smaller sized fat oval: SAT
    FTrue_part1 = abs(FTrue_part1 - WTrue);
    FTrue_part2 = (x.^2 + y.^2) <= (0.08*dim/2)^2; % small bone region
    FTrue = FTrue_part1 + FTrue_part2;
    R2starTrue = zeros(dim);
    dw0True_base = zeros(dim)+dw0shift;
    
    % Additionally add some fat markers near the coil
    %f = (x.^2 + y.^2) <= (0.05*dim/2)^2; % small fat circle
    %FTrue = FTrue + circshift(f,[-3*dim/8 -3*dim/8]) + circshift(f,[3*dim/8 -3*dim/8]) + ...
    %    circshift(f,[-3*dim/8 3*dim/8]) + circshift(f,[3*dim/8 3*dim/8]);

    %apply WF ratio
    %WTrue = WTrue*(1-fatPercentage);
    FTrue = FTrue*fatPercentage;
    WTrue = WTrue-FTrue;
    
    %FTrue = FTrue * 1e5;
    %WTrue = WTrue * 1e5;
    
elseif water_fat_image == 5
    
    disp('Generating true Water and Fat images');
    load(baseline_test_path);
    WTrue = Wlib(:,:,1);
    FTrue = Flib(:,:,1);
    dim = size(WTrue,1);
    
    R2starTrue = R2starlib(:,:,1);
    dw0True_base = dw0lib(:,:,1);
        
else
    disp('ERROR: the water fat image number that was provided is not yet defined');
end

subplot(2,1,1);
imagesc(abs(WTrue)); colorbar; axis image; axis off;
subplot(2,1,2);
imagesc(abs(FTrue)); colorbar; axis image; axis off;
caxis([0 1]);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose a specific hotspot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxtemps = 0:0.0018:0.25; % KEEP 0 as first dyn to get lib with 0 (0:4:40;)

hotspot_ix = 9;       % Options are:
                      % 1. Very small HIFU hotspot in the center of the phantom
                      % 2. Uniform hotspot that is slightly smaller than the water oval
                      % 3. Gradient hotspot that is slightly smaller than the water oval
                      % 4. Uniform hotspot that has the same radius as the water oval
                      % 5. Gradient hotspot that has the same radius as the water oval
                      % 6. Uniform hotspot that has the same radius as the water oval, but the central bone region is excluded
                      % 7. Gradient hotspot that has the same radius as the water oval, but the central bone region is excluded
                      % 8. (small fat + bone) uniform hotspot that has the same radius as the water oval, but the central bone region is excluded
                      % 9. (small fat + bone) uniform hotspot over entire leg
                      % 10. Uniform hotspot positioned in the water part of the loaded leg
                      
%--- generate the hot spot
if hotspot_ix == 1
    
    disp('Generating hotspot');
    [x,y] = meshgrid(-dim/2:dim/2-1);
    hotspot = exp(-hssig2*(x.^2 + y.^2)); % normalized hot spot
    
elseif hotspot_ix == 2
    
    disp('Generating hotspot');
    [x,y] = meshgrid(-dim/2:dim/2-1);
    smaller_radius = 2.5*x.^2 + 3.5*y.^2 <= (0.8*dim/2)^2; % slightly smaller size than water oval
    hotspot = mat2gray(smaller_radius); % Normalize between 0 and 1
    
elseif hotspot_ix == 3
    
    [x,y] = meshgrid(-dim/2:dim/2-1);
    smaller_radius = 2.5*x.^2 + 3.5*y.^2 <= (0.8*dim/2)^2; % same size as water oval
    hotspot = bwdist(abs(smaller_radius-1));
    hotspot = mat2gray(hotspot); % Normalize between 0 and 1
    
elseif hotspot_ix == 4
    
    disp('Generating hotspot');
    [x,y] = meshgrid(-dim/2:dim/2-1);
    same_radius = 1.3*x.^2 + 2.0*y.^2 <= (0.8*dim/2)^2; % same size as water oval
    hotspot = mat2gray(same_radius); % Normalize between 0 and 1
    
elseif hotspot_ix == 5
    
    disp('Generating hotspot');
    [x,y] = meshgrid(-dim/2:dim/2-1);
    same_radius = 1.3*x.^2 + 2.0*y.^2 <= (0.8*dim/2)^2; % same size as water oval
    hotspot = bwdist(abs(same_radius-1));
    hotspot = mat2gray(hotspot); % Normalize between 0 and 1
    
elseif hotspot_ix == 6
    
    disp('Generating hotspot');
    [x,y] = meshgrid(-dim/2:dim/2-1);
    same_radius = 1.3*x.^2 + 2.0*y.^2 <= (0.8*dim/2)^2; % same size as water oval
    bone = (x.^2 + y.^2) <= (0.2*dim/2)^2; % small bone region
    same_radius(bone ~= 0) = 0;
    hotspot = mat2gray(same_radius); % Normalize between 0 and 1
    
elseif hotspot_ix == 7
    
    disp('Generating hotspot');
    [x,y] = meshgrid(-dim/2:dim/2-1);
    same_radius = 1.3*x.^2 + 2.0*y.^2 <= (0.8*dim/2)^2; % same size as water oval
    bone = (x.^2 + y.^2) <= (0.2*dim/2)^2; % small bone region
    same_radius(bone ~= 0) = 0;
    hotspot = bwdist(abs(same_radius-1));
    hotspot = mat2gray(hotspot); % Normalize between 0 and 1
    
elseif hotspot_ix == 8
    
    disp('Generating hotspot');
    [x,y] = meshgrid(-dim/2:dim/2-1);
    %same_radius = 1.15*x.^2 + 1.65*y.^2 <= (0.8*dim/2)^2; % same size as water oval
    same_radius = 1.25*x.^2 + 1.75*y.^2 <= (0.8*dim/2)^2; % Slightly smaller sized fat oval: SAT
    %bone = (x.^2 + y.^2) <= (0.05*dim/2)^2; % small bone region
    bone = (x.^2 + y.^2) <= (0.08*dim/2)^2; % small bone region
    same_radius(bone ~= 0) = 0;
    hotspot = mat2gray(same_radius); % Normalize between 0 and 1    
    
elseif hotspot_ix == 9
    
    disp('Generating hotspot');
    [x,y] = meshgrid(-dim/2:dim/2-1);
    %same_radius = 1.15*x.^2 + 1.65*y.^2 <= (0.8*dim/2)^2; % same size as water oval
    same_radius = x.^2 + 1.5*y.^2 <= (0.8*dim/2)^2; % same size as water oval
    %bone = (x.^2 + y.^2) <= (0.05*dim/2)^2; % small bone region
    %bone = (x.^2 + y.^2) <= (0.08*dim/2)^2; % small bone region
    %same_radius(bone ~= 0) = 0;
    %hotspot = bwdist(abs(same_radius-1));
    hotspot = mat2gray(same_radius); % Normalize between 0 and 1
    
elseif hotspot_ix == 10
    
    disp('Generating hotspot');
    
    summed_WF = abs(WTrue) + abs(FTrue);
    
    WTrue_frac = (abs(WTrue) - min(summed_WF,[],'all')) / (max(summed_WF,[],'all') - min(summed_WF,[],'all'));
    FTrue_frac = (abs(FTrue) - min(summed_WF,[],'all')) / (max(summed_WF,[],'all') - min(summed_WF,[],'all'));
    
    heating_region = WTrue_frac;
    heating_region(heating_region<0.25 & FTrue_frac>0.09) = 0;
    
    heating_region(heating_region>0.05) = 1;
    heating_region(heating_region<=0.05) = 0;
    hotspot = heating_region;
    
else
    disp('ERROR: the hotspot_ix number that was provided is not yet defined');
end

for ii = 1:length(maxtemps)
    hotspotTrue(:,:,ii) = hotspot*alpha*maxtemps(ii)*b0*2*pi*gamma;
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the remaining factors present in the baseline image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---build polynomial matrix A
disp('Forming A polynomial Matrix');

[yc,xc] = meshgrid(linspace(-1/2,1/2,dim));
yc = yc(:);
xc = xc(:);
A = [];
for yp = 0:order
    for xp = 0:(order-yp)
        A = [A (xc.^xp).*(yc.^yp)];
    end
end

%--- get true polynomial coefficients and spatial phase shift
cTrue = zeros(size(A,2),1) + dw0shift;
AcTrue = reshape(A*cTrue,[dim dim]);% polynomial background phase shift

%--- apply gradient to the initial off resonance map
if water_fat_image == 5 
    % Use the off resonance as obtained from the .mat file
    dw0True = dw0True_base;
else
    %cGrad = ones(size(A,2),1) .* (-1+2*rand(size(A,2),1));
    cGrad = [0.0963 ; -0.7295 ; -0.5856 ; 0.7446 ; 0.6307 ; 0.5737 ; 0.4227 ; 0.6904 ; -0.9726 ; 0.8639];
    dw0True = dw0True_base + 1e-10*(reshape(A*cGrad,[dim dim]));
end

%--- create phi static shift
phiTrue = phi0shift; 

%--- generate library
disp('Generating image library');
motionvec = [0 0 0]; % vector of object positions (integer) [-dim/8 0 dim/8]
for ii = 1:length(motionvec)
    imgslib(:,:,:,ii) = awgn(calcmeimgs(WTrue,FTrue,tes,b0,dw0True,R2starTrue,0,prc),snr,'measured');
end

for ii = 1:length(motionvec) 
    
    % shift the images in first spatial dimension to simulate motion
    imgslib(:,:,:,ii) = circshift(imgslib(:,:,:,ii),motionvec(ii));
    
    % if we will not use Graphcut (usetruewf == 1), we will need the following:
    WTruebase(:,:,ii) = circshift(WTrue,motionvec(ii));
    FTruebase(:,:,ii) = circshift(FTrue,motionvec(ii));
    R2starTruebase(:,:,ii) = circshift(R2starTrue,motionvec(ii));
    dw0Truebase(:,:,ii) = circshift(dw0True,motionvec(ii));
    
    % to test the functions using the true values for all other variables,
    % we need the following:
    hotspotTrueshift(:,:,:,:,ii) = circshift(hotspotTrue,motionvec(ii));
        
    AcTrueshift(:,:,ii) = circshift(AcTrue,motionvec(ii));
    phiTrueshift(:,:,ii) = circshift(phiTrue,motionvec(ii));%MEP 20150727
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%
% generate dynamic images
%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Generating Dynamic images');

% Generate the dynamic images using the hotspot as defined in the section 'Choose a specific hotspot'
for ii = 1:length(maxtemps)
    imgsdyn(:,:,:,ii) = awgn(calcmeimgs(WTruebase(:,:,dynmotind).*exp(1i*angle(phiTrueshift(:,:,dynmotind))),FTruebase(:,:,dynmotind).*exp(1i*angle(phiTrueshift(:,:,dynmotind))),...
        tes,b0,dw0Truebase(:,:,dynmotind) + ii*AcTrueshift(:,:,dynmotind),...%MEP 20150727
        R2starTruebase(:,:,dynmotind),hotspotTrueshift(:,:,ii,:,dynmotind),prc),snr,'measured');
end

AcTrue_base = AcTrue;
AcTrue = [];
for ii = 1:length(maxtemps)
    AcTrue(:,:,ii) = ii*AcTrue_base;
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%
% apply averaging
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_averages = 3;

% Apply averaging to the original dynamic images
averaged_idx = 1;
imgsdyn_avg = [];
%for ii = 1:floor(length(maxtemps)/num_averages)
for ii = 1:floor(size(imgsdyn,4)/num_averages)
    imgsdyn_avg(:,:,:,ii) = mean(imgsdyn(:,:,:,averaged_idx:averaged_idx+num_averages-1),4);
    averaged_idx = averaged_idx + num_averages;
end

% Also correct the true hotspot and off resonance for averaging (for evaluation)
averaged_idx = 1;
hotspotTrueshift_avg = [];
for ii = 1:floor(length(maxtemps)/num_averages)
    % Hotspot
    %hotspotTrueshift_avg(:,:,ii,1,:) = mean(hotspotTrueshift(:,:,averaged_idx:averaged_idx+num_averages-1,:,:),3);
    hotspotTrueshift_avg(:,:,ii,1,:) = hotspotTrueshift(:,:,averaged_idx+num_averages-1,:,:);
    
    % Off resonance
    AcTrue_avg(:,:,ii) = mean(AcTrue(:,:,averaged_idx:averaged_idx+num_averages-1),3);
    %AcTrue_avg(:,:,ii) = AcTrue(:,:,averaged_idx+num_averages-1);
    
    averaged_idx = averaged_idx + num_averages;
end

subplot(1,5,1);
imagesc(angle(imgsdyn(:,:,1,1))); colorbar;
subplot(1,5,2);
imagesc(angle(imgsdyn(:,:,1,2))); colorbar;
subplot(1,5,3);
imagesc(angle(imgsdyn(:,:,1,3))); colorbar;
subplot(1,5,4);
imagesc(angle(imgsdyn_avg(:,:,1,1))); colorbar; title('averaged');
subplot(1,5,5);
imagesc(angle(imgsdyn_avg(:,:,1,1)) - angle(imgsdyn(:,:,1,3))); colorbar; title('difference');


%% Determine SNR: defined as the mean / std over dynamics

scan = imgsdyn; % imgsdyn | imgsdyn_avg
SNR = squeeze(mean(permute(scan,[4,1,2,3]),1)./std(permute(scan,[4,1,2,3])));
for i = 1:size(SNR,3)
    subplot(2,ceil(size(SNR,3)/2),i);
    imagesc(abs(SNR(:,:,i))); title(['echo ' num2str(i)]); colorbar;
end


%% Provide a visualization on the phantom that is to be evaluated

% Original image
subplot(2,3,1);
imagesc(abs(WTrue)); axis image; axis off; colorbar; title('Water');
subplot(2,3,2);
imagesc(abs(FTrue)); axis image; axis off; colorbar; title('Fat');
%subplot(2,4,3);
%imagesc(abs(imgslib(:,:,3,2))); axis image; colorbar; title('Combined acc. to equation');
subplot(2,3,4);
imagesc(R2starTruebase(:,:,2)); axis image; axis off; colorbar; title('R2*');
subplot(2,3,5);
imagesc(dw0Truebase(:,:,2)); axis image;axis off;  colorbar; title('Initial off resonance');

% Factors applied
subplot(2,3,3);
imagesc(AcTrueshift(:,:,2)); axis image; axis off; colorbar; title('Drift');
subplot(2,3,6);
imagesc(hotspotTrue(:,:,2)); axis image; axis off; colorbar; title('Heated area');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply water/fat separation to baseline images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Separating baseline WF');
%libsource = 'toolbox';
switch libsource
    
    case 'toolbox'
        disp('using toolbox');
        % Set recon parameters
        R2star_range = [10 100 91];                  % [minimum | maximum | number_for_quantization]
        off_resonance_range = [-900 900 750];     % [minimum | maximum | number_for_quantization]
        itr = 100;                                  % Number of graph cut iterations
        lambda = 0.05;                              % Regularization parameter
        theshold = 0.05;                            % Threshold
        %fatmodel.frequency = [0.59 0.49 -0.5 -1.95 -2.46 -2.68 -3.10 -3.40 -3.80];
        %fatmodel.relAmps = [0.0538 0.01 0.0399 0.014 0.0598 0.0798 0.0598 0.5932 0.0897];
        if ~exist('fatmodel','var')
            algoParams = getReconParams(R2star_range,off_resonance_range,itr,lambda,theshold);
        else
            algoParams = getReconParams(R2star_range,off_resonance_range,itr,lambda,theshold,fatmodel);
        end
        info = ['Data\Processed_DIXON\M',num2str(7,'%04u'),'\DIXON_info.mat']; info = load(info);

        % WF separation is done for each dynamic individually with image size of [ nx | ny | 1 | ncoils | nTE ]
        for ii = 1:size(imgslib,4)
        %for ii = 1:1
            % Image data parameters
            GRE_info.imgdef.slice_orientation_TRA_SAG_COR.uniq = 1;
            GRE_info.imgdef.echo_time.uniq = tes*1e3;
            GRE_info.pardef.Water_Fat_shift_pixels = str2double(info.pardef.Water_Fat_shift_pixels);
            GRE_info.imgdef.pixel_spacing_x_y.uniq = info.imgdef.pixel_spacing_x_y.uniq;
            GRE_info.imgdef.slice_thickness_in_mm.uniq = info.imgdef.slice_thickness.uniq;
            frequency_system = 298.034563*1e6;
            imDataParams = getImDataParams(GRE_info,imgslib(:,:,:,ii),'test',frequency_system,b0,prc);

            % Perform Mixed-Magnitude method initialized with the Graph Cut method
            disp(['.. Executing baseline library dynamic ' num2str(ii) ' out of ' num2str(size(imgslib,4))]);
            outParams = fw_i2cm1i_3pluspoint_hernando_graphcut(imDataParams,algoParams);
            %outParams = fw_i2xm1c_3pluspoint_hernando_mixedfit(imDataParams,algoParams);
            %outParams = fw_i2cm0c_3pluspoint_tsaojiang(imDataParams,algoParams);
            
            % Combine data (water, fat, off-resonance and R2*) of individual dynamics
            Wlib(:,:,ii) = outParams.species(1).amps;
            Flib(:,:,ii) = outParams.species(2).amps;
            dw0lib(:,:,ii) = outParams.fieldmap*2*pi;
            R2starlib(:,:,ii) = outParams.r2starmap;
            
        end
        
        % W_voxel = ones(1,1); F_voxel = ones(1,1); dw0_voxel = zeros(1,1); R2star_voxel = 80*ones(1,1);
        % t = linspace(0,20e-3,500);
        % signal_over_time = calcmeimgs(W_voxel,F_voxel,t,b0,dw0_voxel,R2star_voxel,0,prc);
        % plot(t*1e3,abs(squeeze(signal_over_time)), 'LineWidth',3);
        % xlabel('time [ms]');
        % ylabel('signal');
        % a = get(gca,'XTickLabel');
        % set(gca,'XTickLabel',a,'fontsize',18)
        
    case 'true'
        disp('using true water/fat');
        Wlib = WTruebase;
        Flib = FTruebase;
        dw0lib = dw0Truebase;
        R2starlib = R2starTruebase;
        
end

% Provide a visualization (magnitude and phase) of acquired water/fat separation
for i = 1:size(Wlib,3)
    subplot(4,size(Wlib,3),i);
    imagesc(abs(Wlib(:,:,i))); axis off; title(['Water | dyn ' num2str(i)]); axis image;
    subplot(4,size(Wlib,3),size(Wlib,3)+i); 
    imagesc(abs(Flib(:,:,i))); axis off; title(['Fat | dyn ' num2str(i)]); axis image;
    subplot(4,size(Wlib,3),2*size(Wlib,3)+i);
    imagesc(angle(Wlib(:,:,i))); axis off; title(['Water | dyn ' num2str(i)]); axis image;
    subplot(4,size(Wlib,3),3*size(Wlib,3)+i); 
    imagesc(angle(Flib(:,:,i))); axis off; title(['Fat | dyn ' num2str(i)]); axis image;
end


%% Difference between true and reconstructed W/F

% figure(1);
% subplot(2,3,1);
% imagesc(abs(Wlib(:,:,1))); colorbar; axis image; axis off; title('recon');
% subplot(2,3,2);
% imagesc(abs(WTruebase(:,:,1))); colorbar; axis image; axis off; title('true');
% subplot(2,3,3);
% imagesc(abs(Wlib(:,:,1)).*body_mask - abs(WTruebase(:,:,1)).*body_mask); colorbar; axis image; axis off; title('difference');
% subplot(2,3,4);
% imagesc(abs(Flib(:,:,1))); colorbar; axis image; axis off; title('recon');
% subplot(2,3,5);
% imagesc(abs(FTruebase(:,:,1))); colorbar; axis image; axis off; title('true');
% subplot(2,3,6);
% imagesc(abs(Flib(:,:,1)).*body_mask - abs(FTruebase(:,:,1)).*body_mask); colorbar; axis image; axis off; title('difference');
% 
% figure(2);
% subplot(2,3,1);
% imagesc(angle(Wlib(:,:,1))); colorbar; axis image; axis off; title('recon');
% subplot(2,3,2);
% imagesc(angle(WTruebase(:,:,1))); colorbar; axis image; axis off; title('true');
% subplot(2,3,3);
% imagesc(angle(Wlib(:,:,1)).*body_mask - angle(WTruebase(:,:,1)).*body_mask); colorbar; axis image; axis off; title('difference');
% subplot(2,3,4);
% imagesc(angle(Flib(:,:,1))); colorbar; axis image; axis off; title('recon');
% subplot(2,3,5);
% imagesc(angle(FTruebase(:,:,1))); colorbar; axis image; axis off; title('true');
% subplot(2,3,6);
% imagesc(angle(Flib(:,:,1)).*body_mask - angle(FTruebase(:,:,1)).*body_mask); colorbar; axis image; axis off; title('difference');
% 
% %dw0lib = lib.dw0lib;
% figure(3);
% subplot(2,3,1);
% imagesc(dw0lib(:,:,1)); colorbar; axis image; axis off; title('recon');
% subplot(2,3,2);
% imagesc(dw0Truebase(:,:,1)); colorbar; axis image; axis off; title('true');
% subplot(2,3,3);
% imagesc(dw0Truebase(:,:,1).*body_mask - dw0lib(:,:,1).*body_mask); colorbar; axis image; axis off; title('difference');
% subplot(2,3,4);
% imagesc(R2starlib(:,:,1)); colorbar; axis image; axis off; title('recon');
% subplot(2,3,5);
% imagesc(R2starTruebase(:,:,1)); colorbar; axis image; axis off; title('true');
% subplot(2,3,6);
% imagesc(R2starTruebase(:,:,1).*body_mask - R2starlib(:,:,1).*body_mask); colorbar; axis image; axis off; title('difference');

        
%% Optimalisation of masking parameters - keep on updating parameters until satisfied
% Three parameters are to be optimized for the fat, body and no-signal masks:
%   1. th_fat    - threshold in fat/water image to assign voxels to fat mask 
%   2. no_signal - threshold in respectively water and fat image to define low signal region
%   3. SE        - structure element to erode thick fatty regions

th_fat = 0.01;              % Higher value results in smaller 'fat' mask
no_signal = [0.03 0.03]; 	% Higher value results in larger 'no signal' mask
SE = strel('diamond', 0);   % Larger SE means thinner fat mask, but if too large we keep the original non-eroded mask
body_factor = 1;            % Larger values result in larger body mask

% Compute the fat, body and no signal mask with current parameters
[fat_mask,body_mask,signal_mask] = get_fatty_regions(Wlib(:,:,1), Flib(:,:,1), 1, [dim dim], th_fat, no_signal, SE);

% Provide a visualization
set(0,'DefaultFigureWindowStyle','docked'); % make visualization docked
subplot(1,3,1); imagesc(fat_mask); axis image; set(gca,'XTick',[]); set(gca,'YTick',[]); title('fat mask');
subplot(1,3,2); imagesc(body_mask); axis image; set(gca,'XTick',[]); set(gca,'YTick',[]); title('body mask');
subplot(1,3,3); imagesc(signal_mask); axis image; set(gca,'XTick',[]); set(gca,'YTick',[]); title('no signal mask');
set(0,'DefaultFigureWindowStyle','normal'); % back to undocked visualization


%% Algorithm parameters

% Only algorithm parameters may be altered, which does not hold for scanner 
% and sequence since those parameters are constant

% Selection of fat model and construction of algorithm parameters
num_fat_peaks = 9;                          % Choose an x-peak model, options are:
                                            % [ 3 | 4 | 5 | 6 | 7 | 9 ] 
fatmodel = get_fatmodel(num_fat_peaks);

% Parameter structure: algorithm parameters
algp.order = 5;                 % Polynomial / spherical background phase shift order 
algp.method = 'poly';           % 'poly' for polynomial and 'spher' for spherical harmonics
algp.nmiter = 5;                % # m iterations
algp.nciter = 5;                % # c iterations
algp.mthresh = 0.001;           % heat-induced frequency threshold for 'significant' heat
algp.stopthresh = 1e-5;         % algorithm stopping threshold
algp.modeltest = 0;             % feature where heating is not optimized for and a residual is computed (no heating is expected)
algp.lam = 0;                   % l1 regularization parameter for script version
algp.fat_l1 = 1;                % 0: l1 norm calculated in whole slice | 1: l1 norm only calculated in fats
algp.max = 0.15;              	% Percentage value (between 0 and 0.5) to cap outliers in the heating process (m_postprocess)
algp.sigma = 3;               	% Value for gaussian filter as a final post-processin step of heating
algp.suppress_info = 0;         % Whether command window information is shown per iteration or not
algp.Ac_const = 0;              % Whether we evaluate with (0) or without (1) accounting for any drift
initWithPreviousEstimate = 1;   % Parameter to speed up the algorithm

% Parameter structure: algorithm masking parameters
algp.th_fat = th_fat;           % threshold for fat definition for when l1 norm is only calculated in fats
algp.no_signal = no_signal; 	% threshold for percentage of voxel intensity in water and fat image to specify as voxel without signal
algp.SE = SE;                   % Stucturing Element used as final step to erode the thick fatty regions
algp.body_factor = body_factor; % Value used to vary the threshold for the in-phase image
algp.Nstd_ol = 2;               % Filter out local extreme values in the fat border
algp.neighbours = 0;            % Define local neighbourhood
algp.min_nb = 0;                % Eliminate voxels at the fat-tissue interface from the fat mask (value between 0 and 8)

% Parameter structure: scan parameters
scanp.dim = [size(imgsdyn,1) size(imgsdyn,2)];
scanp.b0 = b0;
scanp.tes = tes;
scanp.prc = prc;

% Parameter structure: baseline library
lib.Wlib = Wlib;
lib.Flib = Flib;
lib.dw0lib = dw0lib;
lib.R2starlib = R2starlib;
lib.fatmodel = fatmodel;        % x-peak fat model (defined prior to water/fat recon)

% Visualize the mask that is used for the salomir fit
mask = remove_outer_border_pixels(fat_mask,abs(imgslib(:,:,1,1)),algp.min_nb);
imagesc(mask); axis off; axis image;


%% Perform iterative optimization to separate the temperature and drift fields
clear lib cfunc mfunc;
disp('Optimization of multi-echo fat-suppressed MR thermometry has started..')

% Updates onto the parameter structures
lib.Wlib = Wlib;
lib.Flib = Flib;
lib.dw0lib = dw0lib;                 	% Update static off resonance | Options: 'dw0lib' or 'dw0corr' 
lib.R2starlib = R2starlib;
lib.fatmodel = fatmodel;
scanp.tes = tes;                       	% Update echo times to the ones used for the optimization part
rad2degC = 1/(2*pi*b0*alpha*gamma);     % Factor to convert temperature to radians

hsmask = (hotspotTrueshift(:,:,end,:,dynmotind)*rad2degC) > 0.01;
hsmask = imclose(hsmask,strel('disk',20));

mfunc = [];

% Provide which dynamics are to be implemented in the optimization
show_ix = 1;    % size(imgsdyn,4)
init_ix = 138;
averaging = 0;  % Specify of averaging of dynamics is used (1) or not (0)

tic
% lam_params = [1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1];
% for ilam = 1:length(lam_params)
%     disp(['lam ' num2str(ilam)]);
%     algp.lam = lam_params(ilam);
for ii = init_ix:init_ix+show_ix-1
    % Select if averaging is adopted or not
    if averaging
        dyn_img = imgsdyn_avg(:,:,:,ii);
        TrueAc = AcTrue_avg;
        TrueT = hotspotTrueshift_avg;
    else
        dyn_img = imgsdyn(:,:,:,ii);
        TrueAc = AcTrue;
        TrueT = hotspotTrueshift;
    end

    disp(['dyn ' num2str(ii)]);
    if initWithPreviousEstimate && ii == init_ix
        [mfunc(:,:,ii),wts(:,ii),Acfunc(:,:,ii),cfunc(:,ii),A,phi(:,:,ii)] = ...
            thermo_hybrid_waterfat(dyn_img, algp, scanp, lib,1,0,zeros(dim,dim));
%         [mfunc(:,:,ii,ilam),wts(:,ii),Acfunc(:,:,ii),cfunc(:,ii),A,phi(:,:,ii)] = ...
%             thermo_hybrid_waterfat(dyn_img, algp, scanp, lib,1,0,zeros(dim,dim));
    elseif initWithPreviousEstimate && ii > init_ix
        [mfunc(:,:,ii),wts(:,ii),Acfunc(:,:,ii),cfunc(:,ii),A,phi(:,:,ii)] = ...
            thermo_hybrid_waterfat(dyn_img, algp, scanp, lib,1,0,mfunc(:,:,ii-1));
%         [mfunc(:,:,ii,ilam),wts(:,ii),Acfunc(:,:,ii),cfunc(:,ii),A,phi(:,:,ii)] = ...
%             thermo_hybrid_waterfat(dyn_img, algp, scanp, lib,1,0,mfunc(:,:,ii-1));
    else
        %phi switch off, m switch on (solve for temperature without Tx/Rx gain)
        [mfunc(:,:,ii),wts(:,ii),Acfunc(:,:,ii),cfunc(:,ii),A,phi(:,:,ii)] = ...
            thermo_hybrid_waterfat(dyn_img, algp, scanp, lib,1,0); 
%         [mfunc(:,:,ii,ilam),wts(:,ii),Acfunc(:,:,ii),cfunc(:,ii),A,phi(:,:,ii)] = ...
%             thermo_hybrid_waterfat(dyn_img, algp, scanp, lib,1,0); 
    end
    lib.cinit = cfunc(:,ii);
    lib.Acinit = Acfunc(:,:,ii);

%         found_temp = mfunc(:,:,ii); found_temp = found_temp(hsmask == 1);
%         true_temp = TrueT(:,:,ii,:,dynmotind); true_temp = true_temp(hsmask == 1);
%         RMSE_func(ilam,ii) = sqrt(abs(mean(found_temp*rad2degC - true_temp*rad2degC)));
%     end
end
toc


%%
dyn = 138;

% Acquire masks
fat_mask_rmse = imbinarize(Flib(:,:,1));
water_mask_rmse = hsmask - imbinarize(Flib(:,:,1));

% 1D array for body, fats and water
found_temp = mfunc(:,:,dyn); 
true_temp = TrueT(:,:,dyn,1,dynmotind); 
found_temp_body = found_temp(body_mask == 1); 
true_temp_body = true_temp(body_mask == 1); 
found_temp_fats = found_temp(fat_mask_rmse == 1); 
true_temp_fats = true_temp(fat_mask_rmse == 1); 
found_temp_water = found_temp(water_mask_rmse == 1);
true_temp_water = true_temp(water_mask_rmse == 1);

% Compute root mean squared error
RMSE_body = sqrt(mean((found_temp_body(:)*rad2degC - true_temp_body(:)*rad2degC).^2));
RMSE_fats = sqrt(mean((found_temp_fats(:)*rad2degC - true_temp_fats(:)*rad2degC).^2));
RMSE_water = sqrt(mean((found_temp_water(:)*rad2degC - true_temp_water(:)*rad2degC).^2));

disp(['body: ' num2str(RMSE_body)]);
disp(['fats: ' num2str(RMSE_fats)]);
disp(['water: ' num2str(RMSE_water)]);


%%
dyn = 138;
fig = figure;
subplot(1,3,1);
imagesc(mfunc(:,:,dyn)*rad2degC); colorbar; title('predicted'); axis image; axis off;
subplot(1,3,2);
imagesc(TrueT(:,:,dyn,:,dynmotind)*rad2degC); colorbar; title('true'); axis image; axis off;
subplot(1,3,3);
imagesc(mfunc(:,:,dyn)*rad2degC - TrueT(:,:,dyn,:,dynmotind)*rad2degC); colorbar; title('difference'); axis image; axis off;

saveas(gcf,['D:\Thesis\Data\Patient water-fat reconstructed data\Water-Fat content evaluation\WaterFatContent_' num2str(fatPercentage) '.png']);

%% Root Mean Squared Error (RMSE) evaluation of water/fat content

composition = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.90,0.95,0.96,0.97,0.98,0.99,1];
body_rmse = [0.079983,0.086386,0.09419,0.096392,0.095243,0.096836,0.098895,0.097531,0.098608,0.095441,0.096757,0.094697,0.092773,0.090841,0.091705,0.089743,0.087967,0.086896,0.1013,0.12074,0.1748,0.32653,1.575,0.13686];
fat_rmse = [0.065227,0.071627,0.079032,0.081123,0.079399,0.080554,0.081634,0.080176,0.080137,0.075586,0.075639,0.071832,0.067748,0.06204,0.057439,0.04539,0.03023,0.0071216,0.12294,0.19575,0.35194,0.72877,3.6379,0.2466];
water_rmse = [0.083004,0.089434,0.09734,0.09957,0.098524,0.1002,0.10245,0.1011,0.10238,0.099445,0.10099,0.099209,0.097621,0.096248,0.097898,0.097113,0.096469,0.096303,0.095644,0.095518,0.095476,0.095309,0.097788,0.095191];

set(gcf,'color','w');
plot(composition(1:end-2),body_rmse(1:end-2),'LineWidth',3);
hold on;
plot(composition(1:end-2),fat_rmse(1:end-2),'LineWidth',3);
hold on;
plot(composition(1:end-2),water_rmse(1:end-2),'LineWidth',3);
hold off;
ax = gca;
ax.FontSize = 16; 
ylabel('RMSE [^0C]','FontSize',20);
xlabel('fat / water composition','FontSize',20);
legend({'body','fat','water'},'Location','northwest','FontSize',14);


%%
% subplot(1,3,1);
% imagesc(Acfunc(:,:,dyn)); colorbar; title('predicted'); axis image; axis off;
% subplot(1,3,2);
% imagesc(TrueAc(:,:,dyn)*rad2degC); colorbar; title('true'); axis image; axis off;
% subplot(1,3,3);
% imagesc(Acfunc(:,:,dyn) - TrueAc(:,:,dyn)); colorbar; title('difference'); axis image; axis off;
% 
% %% Plot the estimated off resonance and temperature frequency for each dynamic
% dynamics = 5;
% plot_ix2 = 1;
% 
% body_mask = body_mask(:,:,1);
% 
% for ii = init_ix:init_ix+dynamics-1
%     subplot(6,dynamics,plot_ix2)
%     imagesc(Acfunc(:,:,ii).*body_mask); axis image; colorbar; title(['Predicted Ac | dynamic ' num2str(ii)]);
%     subplot(6,dynamics,plot_ix2+dynamics)
%     imagesc(TrueAc(:,:,ii).*body_mask); axis image; colorbar; title(['True Ac | dynamic ' num2str(ii)]);
%     subplot(6,dynamics,plot_ix2+2*dynamics)
%     imagesc((Acfunc(:,:,ii)-TrueAc(:,:,ii)).*body_mask); axis image; colorbar; title(['Difference Ac | dynamic ' num2str(ii)]);
%     subplot(6,dynamics,plot_ix2+3*dynamics)
%     min_value = min(TrueT(:,:,ii,:,dynmotind),[],'all'); % [min_value+min_value/10 0.01]
%     imagesc(mfunc(:,:,ii).*body_mask); axis image; colorbar; title(['Predicted hotspot | dynamic ' num2str(ii)]);
%     subplot(6,dynamics,plot_ix2+4*dynamics)
%     imagesc(TrueT(:,:,ii,:,dynmotind).*body_mask); axis image; colorbar; title(['True hotspot | dynamic ' num2str(ii)]);
%     subplot(6,dynamics,plot_ix2+5*dynamics)
%     imagesc((mfunc(:,:,ii) - TrueT(:,:,ii,:,dynmotind)).*body_mask); axis image; colorbar; title(['Difference hotspot | dynamic ' num2str(ii)]);
%     
%     plot_ix2 = plot_ix2 + 1;
% end
% 
% %% Provide a visualization with a step size
% num_plots = [2 5];
% step_size = round(size(mfunc,3)/(num_plots(1)*num_plots(2)));
% 
% figure; plot_ix = 1;
% for ii = init_ix:step_size:size(mfunc,3)
%     subplot(num_plots(1),num_plots(2),plot_ix);
%     imagesc(mfunc(:,:,ii).*body_mask*rad2degC, [-1,2]); axis image; axis off; colorbar; title(['dynamic ' num2str(ii)]);
%     plot_ix = plot_ix + 1;
% end
% 
% figure; plot_ix = 1;
% for ii = init_ix:step_size:size(mfunc,3)
%     subplot(num_plots(1),num_plots(2),plot_ix);
%     imagesc(Acfunc(:,:,ii).*body_mask); axis image; axis off; colorbar; title(['dynamic ' num2str(ii)]);
%     plot_ix = plot_ix + 1;
% end
% 
% %% Lineplot of temperature over time in a specfic pixel
% x_pixels = [floor(size(lib.Wlib,1)/4),floor(3*size(lib.Wlib,1)/4),floor(  size(lib.Wlib,1)/4),floor(3*size(lib.Wlib,1)/4)];
% y_pixels = [floor(size(lib.Wlib,2)/4),floor(  size(lib.Wlib,2)/4),floor(3*size(lib.Wlib,2)/5),floor(3*size(lib.Wlib,2)/4)];
% figure;
% for ii = 1:length(x_pixels)
%     subplot(2,length(x_pixels),ii); im_add = zeros(scanp.dim); im_add(x_pixels(ii),y_pixels(ii)) = 1;
%     imagesc(mat2gray(abs(lib.Wlib(:,:,1)))+im_add); axis image;
%     subplot(2,length(x_pixels),ii+length(x_pixels));
%     plot(squeeze(mfunc(x_pixels(ii),y_pixels(ii),:))*rad2degC);
% end
% 
% 
% %% Phase unwrapping
% 
% % Apply phase unwrapping with ROMEO
% addpath(genpath('D:\ROMEO\matlab'))   % Add the ROMEO software to the directory
% addpath(genpath('D:\Thesis'))         % Add the data folder to the directory
% 
% % Phase unwrapping of the water signal at t=0
% % Set parameter settings and define the output path
% parameters.output_dir = fullfile(['D:\Thesis\Data\Processed_DIXON\M0007\','wrapping_phase_correction']);
% parameters.mask = 'nomask';
% mkdir(parameters.output_dir);
% % Load the magnitude and phase image with dimensions [ x y z echoes ]
% echo = 1;
% dyn = 1;
% disp(['Dynamic ' num2str(dyn) ' of ' num2str(1)]);
% parameters.mag = abs(imgslib(:,:,echo,dyn));
% phase = angle(imgslib(:,:,echo,dyn));
% [unwrapped_water(:,:,dyn), ~] = ROMEO(phase, parameters);
% 
% imagesc(unwrapped_water)
