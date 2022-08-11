function sigma_value = estimate_noise_dt2_max(noisy_imgs)
% 
% Estimate noise for images from "dt2" 
% 
% Noise is estimated in a frame-by-frame approach and MAXIMUM value is
% returned.
% 
% 

Nframes = size(noisy_imgs,3);
sigma_vals = zeros(Nframes,1);

% **
% Noise estimation for each frame
% **

for j = 1:Nframes

    sigma_vals(j) = ceil(evar(noisy_imgs(:,:,j)));% evar + 5
    
end

sigma_value = max(sigma_vals);

fprintf('Estimated noise level: %4.0f \n',sigma_value);

