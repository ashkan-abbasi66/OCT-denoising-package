function [sigma_value,run_time] = estimate_noise_dt1_min(noisy_imgs)
% 
% Estimate noise for images from "dt1" 
% 
% Noise is estimated in a frame-by-frame approach and MINIMUM value is
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
    sigma_vals(j) = floor(function_stdEst(noisy_imgs(:,:,j)));
end

sigma_value = min(sigma_vals);

run_time = toc;

fprintf('Estimated noise level: %4.0f \n',sigma_value);

