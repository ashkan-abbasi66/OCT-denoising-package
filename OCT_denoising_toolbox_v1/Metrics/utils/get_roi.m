% INPUTS:
% imf: image (a filtered image)
% OUTPUTS:
% roi: a vector (1 by d)
% pos: a vector of the form [x y w h] that determines ROI position
% roi_mask: an image of size(imf)
%
% see 'MSR_CNR_metrics_demo.m'
% Ashkan.
function [roi,pos,roi_mask]=get_roi(imf,sz)
if nargin < 2
    h = imrect;
else
    h = imrect(gca,sz);
end
% Constrain drag of imrect within axes limits.
api = iptgetapi(h);
fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),...
   get(gca,'YLim'));
api.setPositionConstraintFcn(fcn);
%
wait(h); % comment this to create mask immediately
roi_mask = createMask(h);
roi=roi_mask.*imf;
roi(~roi_mask)=[];
pos=h.getPosition;
