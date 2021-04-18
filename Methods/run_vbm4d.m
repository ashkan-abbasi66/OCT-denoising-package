function [denoised_imgs,run_time] = run_vbm4d_test(noisy_imgs,params)

% V-BM4D parameter profile
%  'lc' --> low complexity
%  'np' --> normal profile
if isfield(params,'profile')
    profile = params.profile;
else
    profile = 'lc';  
end


if isfield(params,'sigma_value')
    sigma = params.sigma_value;
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

tic

denoised_imgs = ...
    vbm4d( noisy_imgs, sigma, ...
    profile, do_wiener, sharpen, deflicker, verbose );

run_time = toc;


