
function [PSNR,SSIM]=comp_psnr_ssim(gt,im_out,maxI)
% 
% Computes PSNR and SSIM metrics
% 
% Note: 
%   MATLAB's "ssim" function was introduced in MATLAB R2014a.
% 
% INPUTS:
% gt: ground-truth image
% im_out: reconstructed image
% maxI: maximum intensity value (OPTIONAL)
%
% Usage:
%   PSNR=comp_psnr_ssim(gt,im_out,255)
%   [PSNR,SSIM]=comp_psnr_ssim(gt,im_out,255)
%
% 
if ~exist('MaxI','var') || isempty(maxI)
    maxI = -1;
end
    
if maxI == -1
    if max(gt(:))<2
        maxI=1;
    else
        maxI=255;
    end
end

% MSE=mean(mean((gt(:)-im_out(:)).^2));
im_out = im_out(1:size(gt,1),1:size(gt,2));
% gt = gt(1:size(im_out,1),1:size(im_out,2));
MSE=mean(mean((gt(:)-im_out(:)).^2));

PSNR=10*log10((maxI^2)/MSE);

if nargout>1
    SSIM=ssim(im_out,gt,'DynamicRange',maxI); % [ssimval,ssimmap] = ssim(A,ref)
end