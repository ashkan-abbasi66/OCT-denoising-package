function [denoised_imgs,run_time] = run_bm4d_iidnoise_log(noisy_imgs,params)
% 
% Runs BM4D with one noise level which is selected for the whole volume in
% the logarithm domain.
% 
% INPUT
%   noisy_imgs: input noisy OCT volume
%   config: preparation method function handle
%       e.g., config = @config_bm4d_iidnoise_log_dt1
% 
% OUTPUT
%   denoised_imgs: denoised OCT volume
%   run_time: running time of the method for denoising the volume.
% 

n_frames = size(noisy_imgs,3);
if n_frames == 3 
    msg = ['The input volume has 3 frames. \n'...
        'Unfortunately, BM4D cannot process it.' ...
        'We pad one frame to make the processing feasible for BM4D.'];
    warning(sprintf(msg))
    noisy_imgs(:,:,4) = noisy_imgs(:,:,2);
end

sigma_value = params.sigma_value;
sigma_value = take_log(sigma_value/255);
fprintf('Applied noise level in log-domain: %0.4f\n',sigma_value)

tic

% **
% Normalizing each band and taking LOG transformation
% **

normalized_imgs = normalize_bands(noisy_imgs);
noisy_imgs_log = take_log(normalized_imgs);


% **
% Denoising method
% **

denoised_imgs = bm4d(noisy_imgs_log,'Gauss',sigma_value); % newer (v3p2)

denoised_imgs = take_ilog(denoised_imgs);
denoised_imgs = denoised_imgs.*255;

run_time = toc;

denoised_imgs = denoised_imgs(:,:,1:n_frames);



