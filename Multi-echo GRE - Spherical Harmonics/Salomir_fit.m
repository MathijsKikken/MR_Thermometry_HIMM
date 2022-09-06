function [image] = Salomir_fit(image,prev,mask,Nstd_ol,neighbours,min_nb)
%% Multi-echo MR Thermometry in the upper leg at 7T using near-harmonic 2D reconstruction for initialization
%% HIMM: Harmonic Initialized Model-based Multi-echo
% Salomir_fit: function that performs the near-harmonic 2D reconstruction
%
% Creator: Mathijs Kikken (University Medical Center Utrecht)
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% Method for performing a background fit
% Adapted from Reference-Free PRFS.. Salomir et al Feb 2012
% 1) Define the starting point of iteration: 
%    - At the fat border, experimental values are used (at each iteration)
%    - At other voxels, the result of the fitting at the previous time step is used
% 2) Iterate to find background using eq. 13 from paper Salomir 
%
% Inputs:
%   time_step:              index of the image
%   image:                  2D unwrapped (subtracted) phase image
%   prev:                   result of the fit at the previous time step
%   mask:                   fat mask
%   Nstd_ol and neighbours: parameters used to filter out local extreme values in the fat border
%   min_nb:                 parameter used to eliminate voxels at the fat-tissue 
%                           interface from the fat mask) is one at voxels where the body is
%
% Outputs:
%   image:                  result of the background fit 


%% Adapt mask such that voxels on the fat-tissue interface are not considered ground truth 
mask = remove_outer_border_pixels(mask,image,min_nb);

%% Correct for values that differ significantly from their local neighbourhood 
image = image.*mask; % within fat layer, use the experimentally determined values 
image = correct_local_outliers(image,mask,neighbours,Nstd_ol);

%% Provide initialization
% Initialize reference region (provided by mask) and initialize remaining region (by means of the previous fit)
image = image.*mask;
image = image + prev.*(1-mask); % Use previous frame as initial guess

%% Iteration using the discretized formula 
iterations = 0;
change = 1;

% Code needed to create a gif
% h=figure;
% filename = 'Dirichelet.gif';

% Start the background fitting
while iterations < 20000 && change > 1e-9
    
    % Plot the off resonance map that follows a near-harmonic function as gif
%     if rem(iterations,10) == 0
%         imagesc(image); colorbar;
%         axis image off
%         set(h,'Color',[1,1,1]);
%         drawnow
% 
%         % Capture the plot as an image 
%         frame = getframe(h); 
%         im = frame2im(frame); 
%         [imind,cm] = rgb2ind(im,256); 
%         % Write to the GIF File 
%         if iterations == 0 
%             imwrite(imind,cm,filename,'gif', 'Loopcount',0); 
%         else 
%             imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%         end 
%     end
    
    iterations = iterations + 1;
    oldimage = image;
    % Use compiled C code to speed up the algorithm
    % Documentation in the .c files in the MeanAdjacent folder
    [image, change] = MeanAdjacentWithMaskIncChange(oldimage, mask); 
end
end