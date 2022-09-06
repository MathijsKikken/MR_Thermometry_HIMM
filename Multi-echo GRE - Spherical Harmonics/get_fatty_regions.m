function [mask_F,body_mask,mask_no_signal] = get_fatty_regions(Wlib, Flib, wts, dims, th_f_w, th_no_signal, SE, body_factor, body_erode)
%% Multi-echo MR Thermometry in the upper leg at 7T using near-harmonic 2D reconstruction for initialization
%% HIMM: Harmonic Initialized Model-based Multi-echo
% get_fatty_regions: function to create the masks used in the HIMM algorithm
%
% Creator: Mathijs Kikken (University Medical Center Utrecht)
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% Inputs:
%   Wlib:           water image
%   Flib:           fat image
%   wts:            coefficients to combine the library imaging into the image 
%                   best suited to the current dynamic
%   dims:           dimensions of the current data
%   th_f_w:         threshold applied to the relative water/fat image used to 
%                   construct the fat mask
%   th_no_signal:   list with 2 values used to construct the low signal mask 
%                   (through thresholding)
%   SE: structure   element used to erode make the fat mask
%   body_factor:    factor used on top of the Otsu's threshold to define the body mask
%   body_erode:     factor used to construct the structure element to erode the body mask
%
% Output: 
%   mask_F:         fat mask
%   body_mask:      mask of the regions that have tissues
%   mask_no_signal: regions with low signal have value 1


    % Body segmentation of pixels with high Fat / Water content
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Identify region that can be seen as being human body (separate background)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create 'in phase' (IP) image by summing water and fat images over the image library
    IP = zeros([dims(1) dims(2)]);
    for jj = 1:size(Wlib,3)
        IP = IP + Wlib(:,:,jj)*wts(jj) + Flib(:,:,jj)*wts(jj);
    end

    % Otsu's method
    IP = mat2gray(abs(IP));             % Rescale to [0 1] range
    [th_body] = multithresh(IP,2);      % Use Otsu's method to select a number of thresholds

    % Post processing 
    BW = IP; BW(BW < th_body(1)/body_factor) = 0;   % Pixels with value lower than lowest theshold are background
    BW_filled = imfill(BW,'holes');                 % Filling of voxels inside body that were assigned background
    body_mask = bwareaopen(BW_filled, ceil(dims(1)*dims(2)*0.0005)); % Remove objects containing fewer than dim(1)*dim(2)*0.0005 pixels
    body_mask = logical(body_mask);                 % Set nonzero pixels to 1

    % mask_II has more constraints on background pixels as the slightly 
    % larger foreground voxels in background are also assigned background 
    % (note that this effect might be negligable for certain body part such as the leg)
    body_mask = bwareaopen(body_mask, ceil(dims(1)*dims(2)*0.02));
    
    % Apply an erosion operation onto the body_mask
    body_mask = imerode(body_mask,strel('diamond', body_erode));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get fat / water content and apply processing to get fatty region mask
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    water_content = zeros([dims(1) dims(2)]);
    fat_content = zeros([dims(1) dims(2)]);
    for jj = 1:size(Wlib,3)
        water_content = water_content + Wlib(:,:,jj)*wts(jj);
        fat_content = fat_content + Flib(:,:,jj)*wts(jj);
    end

    % Get fat/water content
    rel_fw = abs(fat_content) ./ (abs(water_content) + abs(fat_content)); 

    % Post processing
    rel_fw(isnan(rel_fw)) = 0;          % Set NaN values to background
    rel_fw(body_mask ~= 1) = 0;         % Set all non-human-body values to background
    mask_F = zeros([dims(1) dims(2)]);  % Create an empty mask that is to be filled with fatty regions
    mask_F(rel_fw > th_f_w) = 1;        % All pixels with fat/water content above threshold are set to 1 

    % Remove masked pixels in thick regions and keep in thin regions
    mask_F = bwareaopen(mask_F, ceil(dims(1)*dims(2)*0.005)); % Remove objects containing fewer than dim(1)*dim(2)*0.005 pixels
    mask_F_erode = imerode(mask_F, SE);         % Removal of superficial layers
    mask_F_dilate = imdilate(mask_F_erode, SE);	% Grow regions back
    mask_F_thin = mask_F - mask_F_dilate;       % Regions that did not grow back are thin
    mask_F = mask_F_erode + mask_F_thin;        % Keep the eroded and thin regions 
    
    % Get a mask of regions that do not provide (a lot of) signal
    mask_no_signal = zeros([dims(1) dims(2)]);
    water_scale = mat2gray(abs(water_content)); % Rescale to [0 1] range
    fat_scale = mat2gray(abs(fat_content));     % Rescale to [0 1] range
    mask_no_signal((water_scale < th_no_signal(1)) & (fat_scale < th_no_signal(2))) = 1;

end