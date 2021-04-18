function [denoised_imgs,run_time] = ...
    run_bm4d_mix_log_test(noisy_imgs,params)
% 
% Runs Multiscale BM4D in logarithm domain
% 
% 

params.sigma_value = take_log(params.sigma_value/255);
fprintf('Applied noise level in log-domain: %0.4f\n',params.sigma_value)

n_levels = params.n_levels;
wname = params.wname;


tic

% **
% Normalizing each band and taking LOG transformation
% **

normalized_imgs = normalize_bands(noisy_imgs);
noisy_imgs_log = take_log(normalized_imgs);

% **
% Decompose noisy image
% **

% store noisy_imgs and its lowpass (LP) subbands at all n_levels
lp_coefs = cell(n_levels + 1,1);

lp_coefs{1} = noisy_imgs_log;

for jj = 1:n_levels
    w = dwt3(lp_coefs{jj},wname,'mode','sym');
    %w = dwt3(lp_coefs{jj},wname,'mode','per');
    lp_coefs{jj+1} = w.dec{1,1,1};
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
% lp_coefs_hat{ind} = bm4d(lp_coefs{ind},'Gauss',sigma_value);
lp_coefs_hat{ind} = run_bm4d_iidnoise_test(lp_coefs{ind},params);



% Denoise remaining LP subbands
% Example:
%   Assume that n_levels = 1, then:
%   the first LP subband is denoised before this loop.
%   then, this loop just denoises the input noisy image (lp_coefs{1}) and
%   stores the denoising result in lp_coefs_hat{1}. Next, the denoised LP 
%   subband is replaced with the LP subband of the denoised input image.
% 
for jj = n_levels:-1:1

    % get the latest denoised LP subband
    ind = jj+1;
    tmpA = lp_coefs_hat{ind};

    % deonise the LP subband at the the finer level.
    ind = jj;
    msg = '\n ~~ Denoising lowpass subband at scale %d ... \n';
    fprintf(msg,ind - 1)
% %     lp_coefs_hat{ind} = ...
% %         bm4d(lp_coefs{ind},'Gauss',sigma_value);
    lp_coefs_hat{ind} = run_bm4d_iidnoise_test(lp_coefs{ind},params);



    % replace LP subband at finer scale with the latest denoised one
    wt = dwt3(lp_coefs_hat{ind},wname,'mode','sym');
    %wt = dwt3(lp_coefs_hat{ind},wname,'mode','per');
    
    wt.dec{1,1,1} = tmpA;
    
    lp_coefs_hat{ind} = idwt3(wt);
end
    

denoised_imgs = lp_coefs_hat{1};

denoised_imgs = take_ilog(denoised_imgs);
denoised_imgs = denoised_imgs.*255;

run_time = toc;


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~Another Implementation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This is not good when n_levels > 1. But it is faster due to padding.
% % % % for jj = 1:n_levels
% % % % 
% % % %         lp_noisy = lp_coefs{jj};
% % % %         lp_noisy = make_depth_power_2(lp_noisy);
% % % %         w = dwt3(lp_noisy,wname);
% % % %         lp_coefs{jj+1} = w.dec{1,1,1};
% % % % 
% % % %     end
% % % % 
% % % %     % denoised LP subbands will be saved here:
% % % %     lp_coefs_hat = lp_coefs;
% % % % 
% % % % 
% % % % 
% % % %     % **
% % % %     % Denoise noisy_imgs and its LP subbands
% % % %     % **
% % % % 
% % % %     % Denoise the last LP subband
% % % %     ind = n_levels + 1;
% % % %     msg = '\n ~~ Denoising lowpass subband at scale %d ... \n';
% % % %     fprintf(msg,ind - 1)
% % % %     lp_coefs_hat{ind} = bm4d(lp_coefs{ind},'Gauss',take_log(sigma_value/255));
% % % % 
% % % % 
% % % %     % Denoise remaining LP subbands
% % % %     for jj = n_levels:-1:1
% % % % 
% % % %         % get the latest denoised LP subband
% % % %         ind = jj+1;
% % % %         tmpA = lp_coefs_hat{ind};
% % % % 
% % % %         % deonise the LP subband at the the finer level.
% % % %         ind = jj;
% % % %         msg = '\n ~~ Denoising lowpass subband at scale %d ... \n';
% % % %         fprintf(msg,ind - 1)
% % % % 
% % % %         lp_coefs_hat{ind} = ...
% % % %             bm4d(lp_coefs{ind},'Gauss',take_log(sigma_value/255));
% % % % 
% % % % 
% % % %         % replace LP subband at finer scale with the latest denoised one
% % % %         lp_denoised = lp_coefs_hat{ind};
% % % %         lp_denoised = make_depth_power_2(lp_denoised);
% % % %         wt = dwt3(lp_denoised,wname);
% % % % 
% % % % 
% % % %         % Reconstruct the finer scale with (tmpA) and (details subbands of wt)
% % % %         if ind == 1 
% % % %             wt.dec{1,1,1} = tmpA;
% % % %             lp_coefs_hat{ind} = idwt3(wt);
% % % %         else       
% % % %             tmpA = tmpA(:,:,1:size(lp_denoised,3)/2);
% % % %             wt.dec{1,1,1} = tmpA;
% % % %             lp_coefs_hat{ind} = idwt3(wt);
% % % %         end
% % % %     end
% % % %     
% % % %     
% % % % % ========================================================================
% % % % % ========================================================================
% % % % 
% % % % function imgs = make_depth_power_2(imgs)
% % % % % 
% % % % % Make the 3rd dimension (depth) to be a power of 2
% % % % % 
% % % % 
% % % % nf = size(imgs,3);
% % % % if mod(nf,2) ~= 0 || nf ~= 2^nextpow2(nf) 
% % % %     imgs(:,:,nf + 1:2*nf) = imgs(:,:,end:-1:1);
% % % %     imgs = imgs(:,:,1:pow2(ceil(log2(nf))));
% % % % elseif nf<4
% % % %     imgs(:,:,nf + 1:2*nf) = imgs(:,:,end:-1:1);
% % % % end