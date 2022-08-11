function [roi,mask]=get_roi_pos(img,pos)
% 
% returns ROI vectors from their positions
% 
% INPUT:
%   pos: cell array of ROI positions. The positions are saved in [x y w h]
%   format.
%   img: ROI vectors are extracted based on the given positions from this
%   image.
% 
% OUTPUT:
%   roi: cell array of ROI vectors.
%   roi_mask: cell array of ROI masks. ROI mask is an image where the
%   selected pixels are indicated by 1 and the others are set to 0.
% 
% 
if ~iscell(pos)
    pos={pos};
end
[R,C]=size(img);
Nroi=numel(pos);% Number of ROIs
roi=cell(Nroi,1);
mask=cell(Nroi,1);
for i=1:Nroi
    mask{i}=false([R,C]);
    pos{i}=uint16(pos{i});
    x=pos{i}(1);
    y=pos{i}(2);
    h=pos{i}(4);
    w=pos{i}(3);
    rows=y+1:y+h;
    columns=x+1:x+w;
    rows(rows<1 & rows>R)=[];
    columns(columns<1 & columns>C)=[];
    mask{i}(rows,columns)=1;
    %
    roi{i}=mask{i}.*img;
    roi{i}(~mask{i})=[];
end