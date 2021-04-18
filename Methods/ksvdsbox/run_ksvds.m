function [denoised_imgs,run_time] = run_ksvds(noisy_imgs,config)
% 
% Runs Sparse K-SVD method
% 
  
% OUTPUT
%   denoised_imgs: denoised OCT volume
%   run_time: running time of the method for denoising the volume.
% 


[params,rt1] = config(noisy_imgs);
params3d = params.params3d;


% **
% Denoising method
% **

tic

denoised_imgs = ksvdsdenoise(params3d);

run_time = toc + rt1;

% figure,imshow(denoised_imgs(:,:,1)/255)


