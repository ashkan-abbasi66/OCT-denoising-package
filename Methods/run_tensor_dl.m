function [denoised_imgs,run_time] = run_tensor_dl_test(noisy_imgs,params)
% 
% Runs tensor dictionary learning (Tensor DL / TDL)
% 
% 
% Reference:
% Peng, Yi, et al. "Decomposable nonlocal tensor dictionary learning for 
% multispectral image denoising." Proceedings of the IEEE Conference on 
% Computer Vision and Pattern Recognition. 2014.
% 

% if you want to determine sigma, you need to also determine peak_value
params.nsigma = params.sigma_value;
params.peak_value = 255; % the upper bound of dynamic range. (required)

tic

d = 4;
noisy_imgs = crop_image(noisy_imgs, d);


% **
% Denoising method
% **

[ denoised_imgs, ~, ~, ~, ~ ] = TensorDL(noisy_imgs, params);
% [ denoised_imgs, ~, ~, ~, ~ ] = TensorDL(noisy_imgs);
run_time = toc;


