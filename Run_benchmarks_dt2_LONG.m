%{
for i = 1: 8: 128
    fprintf('%d-%d\n', i, i+ 8)
end
...
17-25
...
65-73
...
105-113

for i = 1: 16: 128
    fprintf('%d-%d\n', i, i+ 16)
end
...
17-33
...
65-81
...
97-113
%}

% TODO
% For each volume, we may reconstruct some frames. Save frames seperately.
% 
% 
% IMPORTANT NOTES:
% 
% HINT: DEALING WITH THE WHOLE VOLUME
% 
% Each volume in this dataset has 128 framse (slices). 3D denoising of this
% whole volume is very time-consuming for some methods. 
% We can instad use batchwise approach to denoise the volume. 
% 
% Example: 
%   With batch size = 16, the method is required to denoise 8 sub-volumes 
%   (8*16=128).
%   25, 75, 105 ==> to compute metrics
% 
%   With batch size = 8, the method is required to denoise 16 sub-volumes
%   21, 69, 109
%   13-29 ==> 1 - 16 (21) 29-21 
% 
% Example: (Batchwise denoising of a whole volume)
%       params.start_frame = 1;
%       params.n_frames = 128;
%       params.batch_size = 16;
%   
%   Assume that you want to save slice #21 as a separate TIF file.
%       params.frame_numbers = ones(length(params.test_indices),1)*[21]% 
% 
% 
% Example:
%   Assume that you want to denoise a subset of slices, say 17 to 32
%   and you want to save slice number 21 as separate TIF file. Set the
%   followings:
%     params.start_frame = 17;
%     params.n_frames = 16; % 17 + 16 - 1 = 32 ==> frame #17 to #32 is used.
%     params.batch_size = -1;
%     params.frame_numbers = ones(length(params.test_indices),1)*[21]
% 
% HINT: NOISE Estimation for "dt2"
% 
% We have provided two functions for noise estimation on images from this
% dataset:
%   estimate_noise_dt2_max
%   estimate_noise_dt2_min
% 
% Both of these methods cannot correctly estimate the noise level for this
% dataset. Therefore, researchers usually increase the estimated noise
% level by adding or multipling constants and experimentally and viauslly
% select an approporiate noise level.
% 

clear
clc

addpath(genpath('./Methods'));
addpath('./Metrics/utils');

% ----------------------------------------------------------
% --- Indices of noisy volumes and frames to be denoised ---

% To evaluate on the whole dataset, use the followings:
load('./Metrics/dt2_random_images')
% % % params.test_indices  = 1:13; % indicates volumes to be denoised
% % % params.frame_numbers = random_images; % each row indicates frame numbers for each volume

% To conduct a fast experiment, use these:
params.test_indices = 1:4; % [3,5]
params.frame_numbers = ones(length(params.test_indices),1)*[21,69, 109]; % [10, 60] for all images
params.frame_numbers = ones(length(params.test_indices),1)*[21];
% ----------------------------------------------------------

% params.valid_rows = 150:512 - 1; % 150:380 - 1 -- TEST 100:
params.valid_rows = 51:562;

dataset_name = 'dt2';

% ===========================================
% ===========================================
% Examples:
%   Denoise the whole volume at once:
%       start_frame = 1;
%       n_frames = 128;
%       batch_size = -1;
%   Denoise part of a volume:
%       start_frame = 50;
%       n_frames = 55;
%       batch_size = -1;
%   Denoise the whole volume in batches:
%       start_frame = 1;
%       n_frames = 100;
%       batch_size = 5;

params.start_frame = 17;%1;
params.n_frames = 16;%128; % >= 2 
params.batch_size = -1;

% ===========================================
% ===========================================

params.save_mat = true;

% Directly used by "benchmark_X_on_..." and "evaluate_metrics_..."
common_params = params;



%% Select a denoising method

% To run a method, set its corresponding variable to 1
KSVDS = 0;
KSVDS_log = 0;
WMF = 0;
WMF_log = 0;
VBM4D = 0;
VBM4D_log = 0;
TENSOR_DL = 0; 
TENSOR_DL_log = 0;
BM4D_org = 0; 
BM4D_org_log = 0; 
BM4D_iid = 0; %%%%
BM4D_iid_log = 0; %%%%
MS_BM4D_dwt3 = 0; 
MS_BM4D_dwt3_log = 1; %%%% 
MS_BM4D_dualtree3 = 0;
MS_BM4D_dualtree3_log = 0; %%%%



%% Sparse K-SVD

if KSVDS == 1
    
    params = common_params;
   
    X = @run_ksvds;
    params.get_params = @get_params_ksvds;
    
%     params.noise_estimator = @(y) estimate_noise_dt2_max(y)*2 - 1;
    params.noise_estimator = @(y) estimate_noise_dt2_max(y) + 2;
    
    output_folder_name = 'benchmark_ksvds_dt2';
    benchmark_X_on_dt2_new(output_folder_name,X,params)
    evaluate_metrics_dt2(output_folder_name,params);

end



%% Sparse K-SVD in the logarithm domain

if KSVDS_log == 1
    
    params = common_params;   
    
    X = @run_ksvds_log;
    params.get_params = @get_params_ksvds;
    
    params.noise_estimator = @(y) estimate_noise_dt2_max(y)*2 - 1;
    
    output_folder_name = 'benchmark_ksvds_log_dt2';
    benchmark_X_on_dt2_new(output_folder_name,X,params)
    evaluate_metrics_dt2(output_folder_name,params);

end



%% WMF 
% **NOTE: Parameters are not good for Topcon**

if WMF == 1

    params = common_params;

    X = @run_wmf;
    params.get_params = @get_params_wmf_dt1; %%%%% Can we use better params?
    
    output_folder_name = 'benchmark_wmf_dt2';
    benchmark_X_on_dt2_new(output_folder_name,X,params)
    params.preprocess = @make_size_even;
    evaluate_metrics_dt2(output_folder_name,params);

end



%% WMF in the logarithm domain
% **NOTE: Parameters are not good for Topcon**

if WMF_log == 1
    
    params = common_params;

    output_folder_name = 'benchmark_wmf_log_dt2';

    X = @run_wmf_log;
    params.get_params = @get_params_wmf_dt1;
    
    benchmark_X_on_dt2_new(output_folder_name,X,params)
    params.preprocess = @make_size_even;
    evaluate_metrics_dt2(output_folder_name,params);

end




%% V-BM4D

if VBM4D == 1
    
    params = common_params;
    
    X = @run_vbm4d;
    % Profile parameter: 'np' or 'lc'
    % 'lc' is fast and we have set it as default.
    % params.method.profile = 'np';
    
    params.noise_estimator = @(y) estimate_noise_dt2_max(y)*2 - 1;
    
    output_folder_name = 'benchmark_vbm4d_dt2';  
    benchmark_X_on_dt2_new(output_folder_name,X,params)
    evaluate_metrics_dt2(output_folder_name,params);

end



%% V-BM4D in the logarithm domain

if VBM4D_log == 1
    
    params = common_params;
    
    X = @run_vbm4d_log;
    % Profile parameter: 'np' or 'lc'
    % 'lc' is fast and we have set it as default.
    % params.method.profile = 'np';
    
    params.noise_estimator = @(y) estimate_noise_dt2_max(y)*2 - 1;
    
    output_folder_name = 'benchmark_vbm4d_log_dt2';
    benchmark_X_on_dt2_new(output_folder_name,X,params)
    evaluate_metrics_dt2(output_folder_name,params);
end



%% Tensor DL
% NOTE: This method does not work on "dt2". I don't know why!. Even with
% bigger sigma_value it does not work.
% You may want to DELETE it for this dataset.

if TENSOR_DL == 1
    
    params = common_params;
    
    X = @run_tensor_dl;
    
    params.noise_estimator = @(y) estimate_noise_dt2_max(y)*2 - 1;
    
    output_folder_name = 'benchmark_tensor_dl_dt2';
    benchmark_X_on_dt2_new(output_folder_name,X,params)
    d = 4;
    params.preprocess = @(x) crop_image(x,d);
    evaluate_metrics_dt2(output_folder_name,params);
    
end



%% Tensor DL in the logarithm domain
% NOTE: This method does not work on "dt2". I don't know why!. Even with
% bigger sigma_value it does not work.
% You may want to DELETE it for this dataset.

if TENSOR_DL_log == 1
    
    params = common_params;  
    
    X = @run_tensor_dl_log;
    
    params.noise_estimator = @(y) estimate_noise_dt2_max(y)*2 - 1;
    
    output_folder_name = 'benchmark_tensor_dl_log_dt2';
    benchmark_X_on_dt2_new(output_folder_name,X,params)
    d = 4;
    params.preprocess = @(x) crop_image(x,d);
    evaluate_metrics_dt2(output_folder_name,params);

end



%% Original BM4D: BM4D with voxel-wise noise estimation

if BM4D_org == 1
    
    params = common_params;  
    
    X = @run_bm4d;
    
    output_folder_name = 'benchmark_bm4d_dt2';
    benchmark_X_on_dt2_new(output_folder_name,X,params)
    evaluate_metrics_dt2(output_folder_name,params);
    
end



%% Original BM4D in the logarithm domain

if BM4D_org_log == 1
    
    params = common_params;
    
    X = @run_bm4d_log;
    
    output_folder_name = 'benchmark_bm4d_log_dt2';  
    benchmark_X_on_dt2_new(output_folder_name,X,params)
    evaluate_metrics_dt2(output_folder_name,params);
    
end


%% BM4D with iid noise assumption -- %%%%

noise_factors = {1};
if BM4D_iid == 1
    for coeff = noise_factors
        
        params = common_params;

        X = @run_bm4d_iidnoise;

%         params.noise_estimator = @(y) estimate_noise_dt2_max(y)*2 - 1;
        params.noise_estimator = @(y) estimate_noise_dt2_max(y)*coeff{1};

        output_folder_name = sprintf('benchmark_bm4d_iidnoise_dt2_times%g', coeff{1});
        benchmark_X_on_dt2_new(output_folder_name,X,params)
        evaluate_metrics_dt2(output_folder_name,params); 
        
    end
end


%% BM4D with iid noise assumption in the logarithm domain -- %%%%

noise_factors = {1};
if BM4D_iid_log == 1
    for coeff = noise_factors
        
        params = common_params;

        X = @run_bm4d_iidnoise_log;

        params.noise_estimator = @(y) estimate_noise_dt2_max(y)*2 - 1;

        output_folder_name = sprintf('benchmark_bm4d_iidnoise_log_dt2_times%g', coeff{1});
        benchmark_X_on_dt2_new(output_folder_name,X,params)
        evaluate_metrics_dt2(output_folder_name,params); 
        
    end
end



%% MS BM4D (dwt3)
% NOTE: wavelet name (wname), number of levels, and noise estimator should
% be tuned in such a way that the method can reasonably outperform KSVDS in
% terms of some metrics such as ENL, TP, and EP.

if MS_BM4D_dwt3 == 1
    for n_levels = 3     % number of decomposition level

        params = common_params;
        
        X = @run_bm4d_mix;
        params.method.n_levels = n_levels;
        params.method.wname = 'sym4';
        
        params.noise_estimator = @(y) estimate_noise_dt2_max(y)*2.3 - 1;
        
        output_folder_name = sprintf('benchmark_bm4d_mix_dt2_%s_%d',...
            params.method.wname, n_levels);
        benchmark_X_on_dt2_new(output_folder_name,X,params)
        evaluate_metrics_dt2(output_folder_name,params);

    end
end



%% MS BM4D (dwt3) in the logarithm domain -- %%%%
% for this dataset, it seems that the log version is not better than the
% spatial domain.
% 

wnames = {'sym4'};
noise_factors = {1, 1.25, 1.5};
if MS_BM4D_dwt3_log == 1
    for coeff = noise_factors
        for n_levels = 1     % number of decomposition level

            params = common_params;

            X = @run_bm4d_mix_log;
            params.method.n_levels = n_levels;
            params.method.wname = wnames{1};

%             params.noise_estimator = @(y) estimate_noise_dt2_max(y)*2.3;
            params.noise_estimator = @(y) estimate_noise_dt2_max(y)*coeff{1};

            output_folder_name = sprintf('benchmark_bm4d_mix_log_dt2_%s_%d_times%g',...
                params.method.wname, n_levels, coeff{1});
            benchmark_X_on_dt2_new(output_folder_name,X,params)
            evaluate_metrics_dt2(output_folder_name,params);

        end
    end
end



%% MS BM4D (dualtree3) 
% Allowed filter banks are:
% ("nearsym5_7"), 'nearsym13_19', 'antonini', or 'legall'

% TODO: 
% We can improve the result with the following command in less than one
% second:
% par.win = 4; par.nsig = 19; par.Beta = 1;
% im_out2 = averaging_nearby_slices_org(denoised_imgs,par);
% 

if MS_BM4D_dualtree3 == 1
    for n_levels = 1     % number of decomposition level

        params = common_params;
        
        X = @run_bm4d_mix_dualtree3;
        params.method.n_levels = n_levels;
        params.method.filter_bank = 'nearsym5_7'; % antonini legall
        
        params.noise_estimator = @(y) estimate_noise_dt2_max(y)*2 - 1;
        
        output_folder_name = sprintf('benchmark_bm4d_mix_dualtree3_dt2_%s_%d',...
            params.method.filter_bank, n_levels);
        benchmark_X_on_dt2_new(output_folder_name,X,params)
        evaluate_metrics_dt2(output_folder_name,params);

    end
end



%% MS BM4D (dualtree3) in the logarithm domain -- %%%%
% Allowed filter banks are:
% ("nearsym5_7"), 'nearsym13_19', 'antonini', or 'legall'

% noise_factors = {1, 1.25, 1.5};
noise_factors = {1};
if MS_BM4D_dualtree3_log == 1
    for coeff = noise_factors
        for n_levels = 1     % number of decomposition level

            params = common_params;

            X = @run_bm4d_mix_dualtree3_log; % run_bm4d_mix_dualtree3_log_fast
            params.method.n_levels = n_levels;
            params.method.filter_bank = 'nearsym5_7'; 

    %         params.noise_estimator = @(y) estimate_noise_dt2_max(y)*2.8;
%             params.noise_estimator = @(y) estimate_noise_dt2_max(y)*2 - 1;
            params.noise_estimator = @(y) estimate_noise_dt2_max(y)*coeff{1};

            output_folder_name = sprintf('benchmark_bm4d_mix_dualtree3_log_dt2_%s_%d_times%g',...
                params.method.filter_bank, n_levels, coeff{1});
            benchmark_X_on_dt2_new(output_folder_name,X,params)
            evaluate_metrics_dt2(output_folder_name,params);

        end
    end
end



