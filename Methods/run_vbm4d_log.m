function [denoised_imgs,run_time] = run_vbm4d_log(noisy_imgs,params)

% V-BM4D parameter profile
%  'lc' --> low complexity
%  'np' --> normal profile
if isfield(params,'profile')
    profile = params.profile;
else
    profile = 'lc';  
end


if isfield(params,'sigma_value')
    sigma = take_log(params.sigma_value/255);
    fprintf('Applied noise level in log-domain: %0.4f\n',sigma)
else
    sigma = -1;
end

do_wiener = 1;       % Wiener filtering
                     %   1 --> enable Wiener filtering
                     %   0 --> disable Wiener filtering
sharpen = 1;         % Sharpening
                     %   1 --> disable sharpening
                     %  >1 --> enable sharpening
deflicker = 1;       % Deflickering
                     %   1 --> disable deflickering
                     %  <1 --> enable deflickering
verbose = 1;         % Verbose mode



% **
% Normalizing each band and taking LOG transformation
% **

normalized_imgs = normalize_bands(noisy_imgs);
noisy_imgs_log = take_log(normalized_imgs);


% **
% Denoising method
% **


tic

denoised_imgs = ...
    vbm4d( noisy_imgs_log, sigma, ...
    profile, do_wiener, sharpen, deflicker, verbose );

denoised_imgs = take_ilog(denoised_imgs);
denoised_imgs = denoised_imgs.*255;

run_time = toc;

