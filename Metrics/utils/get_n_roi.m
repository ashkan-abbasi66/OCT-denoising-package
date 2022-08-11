function [roi,pos,roi_mask]=get_n_roi(img,n)
%
% Display the input image and select "n" regions of interest (ROIs).
%
% INPUT
%   img: the input image
%   n: number of ROIs
%
% OUTPUT
%   roi: cell array of ROI vectors.
%   pos: cell array of ROI positions. The positions are saved in [x y w h]
%   format.
%   roi_mask: cell array of ROI masks. ROI mask is an image where the
%   selected pixels are indicated by 1 and the others are set to 0.
%
% USAGE
%   figure
%   imshow(img);
%   [roi,pos,roi_mask]=get_n_roi(img_noisy,n);
% 
% 

if nargin<2
    n=1;
end

% fh=figure;
% imshow(img);


msg='select ROI number %d & double click';
roi=cell(n,1);
roi_mask=cell(n,1);
pos=cell(n,1);
for i=1:n
    title(sprintf(msg,i));
    [roi{i},pos{i},roi_mask{i}]=get_roi(img);
end

