
function TP = comp_TP(im_out,ref,pos,background_indices)
% 
% Compute texture preservation (TP) quality metrics as described in the
% following paper:
% 
%     Pizurica, Aleksandra, et al. "Multiresolution denoising for optical 
%     coherence tomography: a review and evaluation." 
%     Current Medical Imaging Reviews 4.4 (2008): 270-284.
% 
%     "TP ranges between 0 and 1, and remains close to 0 for filters that 
%     severely flatten the image structures."
% 
% 
% INPUTS:
%   im_out: denoised image
%   ref: reference image can be either noisy image / noise-free image
%   pos: cell array of ROI positions. The positions are saved in [x y w h]
%   format. See "get_n_roi".
%   background_indices: vector indicating which ROIs are background ROI.
%   e.g., background_indices = [1] means that the first region is the only
%   background region.
% 
% OUTPUT:
%   TP: Texture Preservation
% 
% Ashkan
% 

N=numel(pos);

% foreground_indices=setdiff(1:N,background_indices); % desired regon indexes

% Nf=numel(foreground_indices);% number of desired regions

roi=get_roi_pos(im_out,pos);
roi_ref=get_roi_pos(ref,pos);

var_ratio=comp_var_ratio(roi,roi_ref,background_indices);

TP=(1/N)*sqrt(mean(im_out(:))./ mean(ref(:)))*var_ratio;



% ------------------------------------------------------------------------
function var_ratio=comp_var_ratio(roi,roin,background_indices)
n=numel(roi);
var_ratio=0;

% foreground_indices=setdiff(1:n,background_indices); % desired regon indexes
% for j = foreground_indices % desired regions
%     var_ratio=var_ratio+var(roi{j})/var(roin{j});
% end

for j = 1:n
    var_ratio=var_ratio+var(roi{j})/var(roin{j});
end