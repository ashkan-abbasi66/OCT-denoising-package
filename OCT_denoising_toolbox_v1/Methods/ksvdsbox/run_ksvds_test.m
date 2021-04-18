function [denoised_imgs,run_time] = run_ksvds_test(noisy_imgs,params)
% 
% Runs Sparse K-SVD method
% 
  
% OUTPUT
%   denoised_imgs: denoised OCT volume
%   run_time: running time of the method for denoising the volume.
% 

% [params_ksvds,rt] = get_params_ksvds(noisy_imgs,params.noise_estimator);

% **
% Denoising method
% **

tic

denoised_imgs = ksvdsdenoise(params);

run_time = toc;

% figure,imshow(denoised_imgs(:,:,1)/255)


