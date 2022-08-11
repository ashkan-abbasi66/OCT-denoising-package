function sigma_value = estimate_noise_dt2_min(noisy_imgs)
% 
% Estimate noise for images from "dt2" 
% 
% Noise is estimated in a frame-by-frame approach and MINIMUM value is
% returned.
% 
% 

Nframes = size(noisy_imgs,3);
sigma_vals = zeros(Nframes,1);

% **
% Noise estimation for each frame
% **

for j = 1:Nframes 
    
    sigma_vals(j) = floor(evar(noisy_imgs(:,:,j)));
    
end

sigma_value = min(sigma_vals);

fprintf('Estimated noise level: %4.0f \n',sigma_value);

