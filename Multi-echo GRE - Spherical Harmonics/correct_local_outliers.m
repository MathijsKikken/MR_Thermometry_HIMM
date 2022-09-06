function image = correct_local_outliers(image,mask,neighbours,Nstd_ol)
%% Multi-echo MR Thermometry in the upper leg at 7T using near-harmonic 2D reconstruction for initialization
%% HIMM: Harmonic Initialized Model-based Multi-echo
% correct_local_outliers: function that evaluates the drift field in the
% fat mask. The neighboring voxels are evaluated based on mean and std. If
% the condition is not met, the value will be replaced by the average.
%
% Creator: Mathijs Kikken (University Medical Center Utrecht)
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% Inputs:
%   image:       drift field
%   mask:        fat mask
%   neighbours:  defines local neighbourhood
%   Nstd_ol:     filter out local extreme values in the fat border
%
% Outputs:
%   image:       drift field from which outliers are corrected for


original = image;
for row_nr = 1:size(image,1)
    for col_nr = 1:size(image,2)
        if mask(row_nr,col_nr)==1 && row_nr > neighbours && col_nr > neighbours && row_nr <= size(image,1)-neighbours && col_nr <= size(image,2)-neighbours
            local_region = original(row_nr-neighbours:row_nr+neighbours,col_nr-neighbours:col_nr+neighbours);
            local_region(neighbours+1, neighbours+1)= 0; 
            local_region = local_region(local_region~=0);
            mean_l = mean(local_region,'all');
            std_l = std(local_region,[],'all');
            if original(row_nr,col_nr)<mean_l-Nstd_ol*std_l || original(row_nr,col_nr)>mean_l+Nstd_ol*std_l
                image(row_nr,col_nr)=mean_l;
            end
        end
    end
end

end