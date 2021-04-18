function EP=comp_EP(im_out,ref,pos,background_indices)
% 
% Compute edge preservation (EP) quality metrics as described in the
% following paper:
% 
%     Pizurica, Aleksandra, et al. "Multiresolution denoising for optical 
%     coherence tomography: a review and evaluation." 
%     Current Medical Imaging Reviews 4.4 (2008): 270-284.
% 
%     "This EP measure ranges between 0 and 1, having smaller
%     values when the edges inside the ROI are more blurred."
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
%   EP: Edge Preservation
% 
% Ashkan
% 

h= fspecial('laplacian');

N=numel(pos);

foreground_indices=setdiff(1:N,background_indices); % desired regon indexes

Nf=numel(foreground_indices);

[~,mask]=get_roi_pos(im_out,pos);

EPm=zeros(Nf,1);

for i=foreground_indices
    ROI_N = ref.* mask{i};  %noisy image
    ROI_N_Lap = imfilter(ROI_N, h);

    ROI = im_out.* mask{i};   % denoised image
    ROI_Lap = imfilter(ROI, h);
    
    num = corr2((ROI_N_Lap-mean(ROI_N_Lap(:))),(ROI_Lap-mean(ROI_Lap(:))));
    %             denum = sqrt(num.*corr2((ROI1_Lap-mean(ROI1_Lap(:))),(ROI1_Lap-mean(ROI1_Lap(:)))));
    denum = sqrt(corr2((ROI_N_Lap-mean(ROI_N_Lap(:))),(ROI_N_Lap-mean(ROI_N_Lap(:)))).*...
        corr2((ROI_Lap-mean(ROI_Lap(:))),(ROI_Lap-mean(ROI_Lap(:)))));
    EPm(i)=num./denum; 
end

EPm=EPm(foreground_indices);
EP=mean(EPm,1);