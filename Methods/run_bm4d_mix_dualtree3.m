function [denoised_imgs,run_time] = ...
    run_bm4d_mix_dualtree3_test(noisy_imgs,params)
% 
% Runs Multiscale BM4D
% Implemented using dualtree3
% 
% INPUT
%   noisy_imgs: input noisy OCT volume - all three dimensions must be even
%   and greater than or equal to 4.
% 
%   n_levels (optional): number of decomposition level. It is different
%   than the LEVEL argument in DUALTREE3. We apply the dualtree3 multiple
%   (n_levels) times. In each time, LEVEL argument in dualtree3 is set to
%   2.
% 
%   filter_bank (optional): any compatible filterbank name in DUALTREE3.
% 
% OUTPUT
%   denoised_imgs: denoised OCT volume
%   run_time: running time of the method for denoising the volume.
% 



% "params.sigma_value" will be used later in this function
n_levels = params.n_levels;
filter_bank = params.filter_bank;


if ~exist('filter_bank','var') || isempty('filter_bank') 
    filter_bank = 'antonini';
    % Allowed filter banks are:
    % ("nearsym5_7"), 'nearsym13_19', 'antonini', or 'legall'
end


tic

% **
% Decompose noisy image and its LP subbands
% **

lp_coefs = cell(n_levels + 1,1);

lp_coefs{1} = noisy_imgs;
bad_dims = 0; 
k = 1;
for jj = 1:n_levels
    
    lp_noisy = lp_coefs{jj};
    if size(lp_noisy ,3) ==2
        bad_dims(k) = jj;
        k = k + 1;
    end
    lp_noisy = make_depth_power_2(lp_noisy);
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
    
% %     lp_coefs_hat{ind} =...
% %         bm4d(lp_coefs{ind},'Gauss',sigma_value);
    lp_coefs_hat{ind} = run_bm4d_iidnoise(lp_coefs{ind},params);

    
    % **
    % Replace the LP subband at finer scale with the latest denoised one
    % **
    
    % Get the highpass subband (ho_hat) from the finer scale 
    lp_denoised = lp_coefs_hat{ind};
    lp_denoised = make_depth_power_2(lp_denoised);
    [~,ho_hat] = dualtree3(lp_denoised,2,'LevelOneFilter',filter_bank);
    
    
    % Reconstruct the finer scale with (tmpA) and (ho_hat)
    if ind == 1 
        if any(ind + 1 == bad_dims)
            lp_coefs_hat{ind} = idualtree3(tmpA(:,:,1:2),ho_hat);
        else
            lp_coefs_hat{ind} = idualtree3(tmpA,ho_hat);
        end            
    else       
        tmpA = tmpA(:,:,1:size(lp_denoised,3)/2);
        lp_coefs_hat{ind} = idualtree3(tmpA,ho_hat);
    end
    
end


denoised_imgs = lp_coefs_hat{1};


run_time = toc;




% ========================================================================
% ========================================================================

% % function imgs = make_depth_power_2(imgs)
% % % 
% % % Make the 3rd dimension (depth) to be a power of 2
% % % 
% % 
% % nf = size(imgs,3);
% % if mod(nf,2) ~= 0 || nf ~= 2^nextpow2(nf) 
% %     imgs(:,:,nf + 1:2*nf) = imgs(:,:,end:-1:1);
% %     imgs = imgs(:,:,1:pow2(ceil(log2(nf))));
% % elseif nf<4
% %     imgs(:,:,nf + 1:2*nf) = imgs(:,:,end:-1:1);
% % end
