function [denoised_imgs,run_time] = run_bm4d_iidnoise_test(noisy_imgs,params)
% 
% Runs BM4D method with intra-frame noise estimation. One noise level is
% selected for the whole volume.
% 
% INPUT
%   noisy_imgs: input noisy OCT volume
%   config: preparation method function handle
%       e.g., config = @config_bm4d_iidnoise_dt1
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

% **
% Denoising method
% **
tic

denoised_imgs = bm4d(noisy_imgs,'Gauss',sigma_value); % newer (v3p2)

run_time = toc;

denoised_imgs = denoised_imgs(:,:,1:n_frames);


