% TODO
% For each volume, we may reconstruct some frames. Save frames seperately.
% 
% 
% IMPORTANT NOTES:
% 
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
params.test_indices  = 1:13; % indicates volumes to be denoised
params.frame_numbers = random_images; % each row indicates frame numbers for each volume

% To conduct a fast experiment, use these:
params.test_indices = [3]; % [3,5]
params.frame_numbers = ones(length(params.test_indices),1)*[10]; % [10, 60]
% ----------------------------------------------------------

% Number of frames to be denoised.
params.n_frames = 5; % >= 2 - Defualt: 5

params.valid_rows = 150:512 - 1; % 150:380 - 1

params.save_mat = false;

% Directly used by "benchmark_X_on_..." and "evaluate_metrics_..."
common_params = params;



%% Select a denoising method

% To run a method, set its corresponding variable to 1
KSVDS = 1;
KSVDS_log = 1;
WMF = 0;
WMF_log = 0;
VBM4D = 0;
VBM4D_log = 0;
TENSOR_DL = 0; 
TENSOR_DL_log = 0;
BM4D_org = 0; 
BM4D_org_log = 0; 
BM4D_iid = 1;
BM4D_iid_log = 0;
MS_BM4D_dwt3 = 1; 
MS_BM4D_dwt3_log = 0;
MS_BM4D_dualtree3 = 1;
MS_BM4D_dualtree3_log = 0;



%% Sparse K-SVD

if KSVDS == 1
    
    params = common_params;
   
    X = @run_ksvds;
    params.get_params = @get_params_ksvds;
    
    params.noise_estimator = @(y) estimate_noise_dt2_max(y)*2 - 1;
    
    output_folder_name = 'benchmark_ksvds_dt2';
    benchmark_X_on_dt2(output_folder_name,X,params)
    evaluate_metrics_dt2(output_folder_name,params);

end



%% Sparse K-SVD in the logarithm domain

if KSVDS_log == 1
    
    params = common_params;   
    
    X = @run_ksvds_log;
    params.get_params = @get_params_ksvds;
    
    params.noise_estimator = @(y) estimate_noise_dt2_max(y)*2 - 1;
    
    output_folder_name = 'benchmark_ksvds_log_dt2';
    benchmark_X_on_dt2(output_folder_name,X,params)
    evaluate_metrics_dt2(output_folder_name,params);

end



%% WMF 
% **NOTE: Parameters are not good for Topcon**

if WMF == 1

    params = common_params;

    X = @run_wmf;
    params.get_params = @get_params_wmf_dt1; %%%%% Can we use better params?
    
    output_folder_name = 'benchmark_wmf_dt2';
    benchmark_X_on_dt2(output_folder_name,X,params)
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
    
    benchmark_X_on_dt2(output_folder_name,X,params)
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
    benchmark_X_on_dt2(output_folder_name,X,params)
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
    benchmark_X_on_dt2(output_folder_name,X,params)
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
    benchmark_X_on_dt2(output_folder_name,X,params)
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
    benchmark_X_on_dt2(output_folder_name,X,params)
    d = 4;
    params.preprocess = @(x) crop_image(x,d);
    evaluate_metrics_dt2(output_folder_name,params);

end



%% Original BM4D: BM4D with voxel-wise noise estimation

if BM4D_org == 1
    
    params = common_params;  
    
    X = @run_bm4d;
    
    output_folder_name = 'benchmark_bm4d_dt2';
    benchmark_X_on_dt2(output_folder_name,X,params)
    evaluate_metrics_dt2(output_folder_name,params);
    
end



%% Original BM4D in the logarithm domain

if BM4D_org_log == 1
    
    params = common_params;
    
    X = @run_bm4d_log;
    
    output_folder_name = 'benchmark_bm4d_log_dt2';  
    benchmark_X_on_dt2(output_folder_name,X,params)
    evaluate_metrics_dt2(output_folder_name,params);
    
end


%% BM4D with iid noise assumption

if BM4D_iid == 1
    
    params = common_params;
    
    X = @run_bm4d_iidnoise;
    
    params.noise_estimator = @(y) estimate_noise_dt2_max(y)*2 - 1;
    
    output_folder_name = 'benchmark_bm4d_iidnoise_dt2';
    benchmark_X_on_dt2(output_folder_name,X,params)
    evaluate_metrics_dt2(output_folder_name,params); 
    
end


%% BM4D with iid noise assumption in the logarithm domain

if BM4D_iid_log == 1
    
    params = common_params;
    
    X = @run_bm4d_iidnoise_log;
    
    params.noise_estimator = @(y) estimate_noise_dt2_max(y)*2 - 1;
    
    output_folder_name = 'benchmark_bm4d_iidnoise_log_dt2';
    benchmark_X_on_dt2(output_folder_name,X,params)
    evaluate_metrics_dt2(output_folder_name,params); 
    
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
        benchmark_X_on_dt2(output_folder_name,X,params)
        evaluate_metrics_dt2(output_folder_name,params);

    end
end



%% MS BM4D (dwt3) in the logarithm domain
% for this dataset, it seems that the log version is not better than the
% spatial domain.
% 

if MS_BM4D_dwt3_log == 1
    for n_levels = 3     % number of decomposition level

        params = common_params;
        
        X = @run_bm4d_mix_log;
        params.method.n_levels = n_levels;
        params.method.wname = 'sym4';
        
        params.noise_estimator = @(y) estimate_noise_dt2_max(y)*2.3;
        
        output_folder_name = sprintf('benchmark_bm4d_mix_log_dt2_%s_%d',...
            params.method.wname, n_levels);
        benchmark_X_on_dt2(output_folder_name,X,params)
        evaluate_metrics_dt2(output_folder_name,params);

    end
end



%% MS BM4D (dualtree3) *******
% Allowed filter banks are:
% ("nearsym5_7"), 'nearsym13_19', 'antonini', or 'legall'

if MS_BM4D_dualtree3 == 1
    for n_levels = 3     % number of decomposition level

        params = common_params;
        
        X = @run_bm4d_mix_dualtree3;
        params.method.n_levels = n_levels;
        params.method.filter_bank = 'nearsym5_7'; % antonini legall
        
        params.noise_estimator = @(y) estimate_noise_dt2_max(y)*2 - 1;
        
        output_folder_name = sprintf('benchmark_bm4d_mix_dualtree3_dt2_%s_%d',...
            params.method.filter_bank, n_levels);
        benchmark_X_on_dt2(output_folder_name,X,params)
        evaluate_metrics_dt2(output_folder_name,params);

    end
end



%% MS BM4D (dualtree3) in the logarithm domain
% Allowed filter banks are:
% ("nearsym5_7"), 'nearsym13_19', 'antonini', or 'legall'

if MS_BM4D_dualtree3_log == 1
    for n_levels = 3     % number of decomposition level

        params = common_params;
        
        X = @run_bm4d_mix_dualtree3_log;
        params.method.n_levels = n_levels;
        params.method.filter_bank = 'antonini'; 
        
        params.noise_estimator = @(y) estimate_noise_dt2_max(y)*2.8;
        
        output_folder_name = sprintf('benchmark_bm4d_mix_dualtree3_log_dt2_%s_%d',...
            params.method.filter_bank, n_levels);
        benchmark_X_on_dt2(output_folder_name,X,params)
        evaluate_metrics_dt2(output_folder_name,params);

    end
end



