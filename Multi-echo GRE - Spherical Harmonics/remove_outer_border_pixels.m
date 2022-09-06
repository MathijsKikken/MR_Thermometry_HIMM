function mask = remove_outer_border_pixels(mask, image, min_nb)
%% Multi-echo MR Thermometry in the upper leg at 7T using near-harmonic 2D reconstruction for initialization
%% HIMM: Harmonic Initialized Model-based Multi-echo
% remove_outer_border_pixels: function that evaluates the fat mask and
% removes voxels that do not have the required number of neighbors. 
%
% Creator: Mathijs Kikken (University Medical Center Utrecht)
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% Inputs:
%   image:       drift field
%   mask:        fat mask
%   min_nb:      eliminates voxels at the fat-tissue interface from the 
%                fat mask (value between 0 and 8)
%
% Outputs:
%   image:       fat mask from which outliers are removed


filter = [1,1,1;1,0,1;1,1,1];
convolution = conv2(filter,mask);
new_mask = zeros(size(mask));
for row_nr = 1:size(image,1)
    for col_nr = 1:size(image,2)
        if mask(row_nr,col_nr)==1
            if min_nb-convolution(row_nr+1,col_nr+1) <= 0.01
                new_mask(row_nr,col_nr)=1; 
            end
        end
    end
end
mask = logical(new_mask); 

end