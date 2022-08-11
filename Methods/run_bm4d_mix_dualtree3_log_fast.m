function [denoised_imgs,run_time] = ...
    run_bm4d_mix_dualtree3_log_fast(noisy_imgs,params)



sigma_value = params.sigma_value;
params.sigma_value = take_log(sigma_value/255);
fprintf('Applied noise level in log-domain: %0.4f\n',params.sigma_value)

n_levels = params.n_levels;
filter_bank = params.filter_bank;


if ~exist('filter_bank','var') || isempty('filter_bank') 
    filter_bank = 'antonini';
    % Allowed filter banks are:
    % ("nearsym5_7"), 'nearsym13_19', 'antonini', or 'legall'
end


tic

% **
% Normalizing each band and taking LOG transformation
% **

normalized_imgs = normalize_bands(noisy_imgs);
noisy_imgs_log = take_log(normalized_imgs);

noisy_imgs_log = make_size_even(noisy_imgs_log);

% **
% Decompose noisy image and its LP subbands
% **

lp_coefs = cell(n_levels + 1,1);

lp_coefs{1} = noisy_imgs_log;
bad_dims = 0;
k = 1;
for jj = 1:n_levels
    
    lp_noisy = lp_coefs{jj};
    if size(lp_noisy ,3) ==2
        bad_dims(k) = jj;
        k = k + 1;
    end
    lp_noisy = make_depth_power_2_fast(lp_noisy); %%%% FAST %%% LESS PADDING
    [lo,~] = dualtree3(lp_noisy,2,'LevelOneFilter',filter_bank);
    lp_coefs{jj+1} = lo;
    
end

% denoised LP subbands will be saved here:
lp_coefs_hat = lp_coefs;


% **
% Denoise noisy_imgs and its LP subbands
% **

% Denoise the last LP subband
ind = n_levels + 1;
msg = '\n ~~ Denoising lowpass subband at scale %d ... \n';
fprintf(msg,ind - 1)

% lp_coefs_hat{ind} = ...
%     bm4d(lp_coefs{ind},'Gauss',sigma_value);
lp_coefs_hat{ind} = run_bm4d_iidnoise(lp_coefs{ind},params);



% Denoise the remaining LP subbands

for jj = n_levels:-1:1
    
    % get the latest denoised LP subband
    ind = jj+1;
    tmpA = lp_coefs_hat{ind};
    
    % deonise the LP subband at the the finer level.
    ind = jj;
    msg = '\n ~~ Denoising lowpass subband at scale %d ... \n';
    fprintf(msg,ind - 1)
    
%     lp_coefs_hat{ind} =...
%         bm4d(lp_coefs{ind},'Gauss',sigma_value);
    lp_coefs_hat{ind} = run_bm4d_iidnoise(lp_coefs{ind},params);

    
    % **
    % Replace the LP subband at finer scale with the latest denoised one
    % **
    
    % Get the highpass subband (ho_hat) from the finer scale 
    lp_denoised = lp_coefs_hat{ind};
    %lp_denoised = make_depth_power_2(lp_denoised);
    [~,ho_hat] = dualtree3(lp_denoised,2,'LevelOneFilter',filter_bank);
    
    
    % Reconstruct the finer scale with (tmpA) and (ho_hat)
    if ind == 1 
        if any(ind + 1 == bad_dims)
            lp_coefs_hat{ind} = idualtree3(tmpA(:,:,1:2),ho_hat);
        else
            lp_coefs_hat{ind} = idualtree3(tmpA,ho_hat);
        end   
%         if ind + 1 == bad_dims
%             lp_coefs_hat{ind} = idualtree3(tmpA(:,:,1:2),ho_hat);
%         else
%             lp_coefs_hat{ind} = idualtree3(tmpA,ho_hat);
%         end
    else       
        if mod(size(lp_denoised,3)/2,2) == 0
            tmpA = tmpA(:,:,1:size(lp_denoised,3)/2);  % ## DEFUALT ##        
        end
        lp_coefs_hat{ind} = idualtree3(tmpA,ho_hat);
    end
    
end


denoised_imgs = lp_coefs_hat{1};

denoised_imgs = take_ilog(denoised_imgs);
denoised_imgs = denoised_imgs.*255;

run_time = toc;



