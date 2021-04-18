function [denoised_imgs,run_time] = run_ksvds_log_test(noisy_imgs,params)
% 
% Runs Sparse K-SVD method
% 
% 

% [params_ksvds,rt] = get_params_ksvds(noisy_imgs,params.noise_estimator);


% **
% Normalizing each band and taking LOG transformation
% **

normalized_imgs = normalize_bands(noisy_imgs);
noisy_imgs_log = take_log(normalized_imgs);

params.x = noisy_imgs_log;
params.sigma = take_log(params.sigma/255);

fprintf('Applied noise level in log-domain: %0.4f\n',params.sigma)

% **
% Denoising method
% **

tic

denoised_imgs = ksvdsdenoise(params);

denoised_imgs = take_ilog(denoised_imgs);
denoised_imgs = denoised_imgs.*255;

run_time = toc;

% figure,imshow(denoised_imgs(:,:,1)/255)


