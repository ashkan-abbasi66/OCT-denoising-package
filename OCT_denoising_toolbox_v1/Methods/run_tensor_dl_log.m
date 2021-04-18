function [denoised_imgs,run_time] = run_tensor_dl_log_test(noisy_imgs,params)
% 
% Runs tensor dictionary learning (Tensor DL / TDL) in the logarithm
% transform domain
% 
% see "run_tensor_dl"
% 
% 

% if you want to determine sigma, you need to also determine peak_value
params.nsigma = params.sigma_value;
params.nsigma = take_log(params.nsigma/255);
fprintf('Applied noise level in log-domain: %0.4f\n',params.nsigma)
params.peak_value = 1; % the upper bound of dynamic range. (required)

tic

d = 4;
noisy_imgs = crop_image(noisy_imgs, d);

% **
% Normalizing each band and taking LOG transformation
% **

normalized_imgs = normalize_bands(noisy_imgs);
noisy_imgs_log = take_log(normalized_imgs);


% **
% Denoising method
% **

[ denoised_imgs, ~, ~, ~, ~ ] = ...
    TensorDL(noisy_imgs_log, params);
denoised_imgs = take_ilog(denoised_imgs);
denoised_imgs = denoised_imgs.*255;

run_time = toc;


