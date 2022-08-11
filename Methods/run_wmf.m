function [denoised_imgs,run_time] = run_wmf_test(noisy_imgs,params)
% 
% Runs wavelet multi-frame denoising (WMF)
% 
% INPUT
%   noisy_imgs: input noisy OCT volume
%   
% OUTPUT
%   denoised_imgs: denoised OCT volume
%   run_time: running time of the method for denoising the volume.
% 
% USAGE:
%   [denoised_imgs,run_time] = run_wmf(noisy_imgs,config);
% 


% Wavelet is used here. Therefore, the input image size must be 
% divisible by a certain number.
if isfield(params,'sizeDivisibleBy')
    d = params.sizeDivisibleBy;
else
    d = -1;
end

tic

% pad image if needed
if d ~= -1
    [noisy_imgs,valid_rows,valid_cols] = symextend_3d(noisy_imgs,d);
end

denoised_imgs = waveletMultiFrame(...
    noisy_imgs, ...
    'k', params.k, ...
    'p', params.p, ...
    'maxLevel', params.maxLevel, ...
    'weightMode', params.weightMode, ... 
    'basis', params.basis);

if d ~= -1
    denoised_imgs = denoised_imgs(valid_rows,valid_cols,:);  
end

run_time = toc;

