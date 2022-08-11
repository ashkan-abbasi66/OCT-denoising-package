function [sigma_value,run_time] = estimate_noise_dt1_max(noisy_imgs)
% 
% Estimate noise for images from "dt1" 
% 
% Noise is estimated in a frame-by-frame approach and MAXIMUM value is
% returned.
% 
% 

Nframes = size(noisy_imgs,3);
sigma_vals = zeros(Nframes,1);

tic

% **
% Noise estimation for each frame
% **

for j = 1:Nframes
    sigma_vals(j) = ceil(function_stdEst(noisy_imgs(:,:,j)));
end

sigma_value = max(sigma_vals);

run_time = toc;

fprintf('Estimated noise level: %4.0f \n',sigma_value);

