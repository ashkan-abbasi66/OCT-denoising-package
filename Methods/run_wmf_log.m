function [denoised_imgs,run_time] = run_wmf_log_test(noisy_imgs,params)
% 
% Runs wavelet multi-frame denoising (WMF) in the logarithm domain
% 
% INPUT
%   noisy_imgs: input noisy OCT volume
%   config: function handle used for obtaining the required parameters
%   
% OUTPUT
%   denoised_imgs: denoised OCT volume
%   run_time: running time of the method for denoising the volume.
% 
% USAGE:
%   config = @config_wmf_log_dt1;
%   [denoised_imgs,run_time] = run_wmf_log(noisy_imgs,config);
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

normalized_noisy_imgs = normalize_bands(noisy_imgs);
noisy_imgs_log = take_log(normalized_noisy_imgs);


denoised_imgs = waveletMultiFrame(...
    noisy_imgs_log, ...
    'k', params.k, ...
    'p', params.p, ...
    'maxLevel', params.maxLevel, ...
    'weightMode', params.weightMode, ... 
    'basis', params.basis);

if d ~= -1
    denoised_imgs = denoised_imgs(valid_rows,valid_cols,:);  
end

denoised_imgs = take_ilog(denoised_imgs);
denoised_imgs = denoised_imgs.*255;

run_time = toc;
