function [denoised_imgs,run_time] = run_bm4d_mix(noisy_imgs,params)
% 
% Runs Multiscale BM4D via denoising and mixing lowpass(LP) subbands across
% scales.
% 
% 

% Note: "params.sigma_value" is used later by "run_bm4d_iidnoise"
n_levels = params.n_levels;
wname = params.wname;

tic 

% **
% Decompose noisy image
% **

% store noisy_imgs and its lowpass (LP) subbands at all n_levels
lp_coefs = cell(n_levels + 1,1);

lp_coefs{1} = noisy_imgs;
for jj = 1:n_levels
    w = dwt3(lp_coefs{jj},wname);
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
lp_coefs_hat{ind} = run_bm4d_iidnoise(lp_coefs{ind},params);


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
%     lp_coefs_hat{ind} = bm4d(lp_coefs{ind},'Gauss',sigma_value);
    lp_coefs_hat{ind} = run_bm4d_iidnoise(lp_coefs{ind},params);


    
    % replace LP subband at finer scale with the latest denoised one
    wt = dwt3(lp_coefs_hat{ind},wname);
    wt.dec{1,1,1} = tmpA;
    lp_coefs_hat{ind} = idwt3(wt);
end


denoised_imgs = lp_coefs_hat{1};

run_time = toc;

end



