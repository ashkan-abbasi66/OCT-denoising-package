% TODO
%   "evaluate_metrics_dt3"
%   

clear
clc

addpath(genpath('./Methods'));
addpath('./Metrics/utils');

dataset_name = 'dt3';

% ===========================================
% ===========================================
% Examples:
%   Denoise the whole volume at once:
%       start_frame = 1;
%       n_frames = 100;
%       batch_size = -1;
%   Denoise part of a volume:
%       start_frame = 50;
%       n_frames = 55;
%       batch_size = -1;
%   Denoise the whole volume in batches:
%       start_frame = 1;
%       n_frames = 100;
%       batch_size = 5;

params.start_frame = 1;
params.n_frames = 100; % >= 2 
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
BM4D_iid = 0; %%%
BM4D_iid_log = 0;
MS_BM4D_dwt3 = 0; 
MS_BM4D_dwt3_log = 0; %% 
MS_BM4D_dualtree3 = 0;
MS_BM4D_dualtree3_log = 1; %%%



%% Sparse K-SVD

if KSVDS == 1
    
    params = common_params;
    
    X = @run_ksvds;
    params.get_params = @get_params_ksvds;
    
    params.noise_estimator = @(y) estimate_noise_dt1_max(y); % no multiplication
    
    output_folder_name = sprintf('benchmark_ksvds_%s',dataset_name);
    benchmark_X_on_dt3(output_folder_name,X,params)
    %%%evaluate_metrics_dt3(output_folder_name,params);

end


%% Sparse K-SVD in the logarithm domain

if KSVDS_log == 1
    
    params = common_params;   
        
    X = @run_ksvds_log;
    params.get_params = @get_params_ksvds;
    
    params.noise_estimator = @(y) estimate_noise_dt1_min(y);
    
    output_folder_name = sprintf('benchmark_ksvds_log_%s',dataset_name);
    benchmark_X_on_dt3(output_folder_name,X,params)
    %%%evaluate_metrics_dt3(output_folder_name,params);

end



%% WMF

if WMF == 1

    params = common_params;

    X = @run_wmf;
    params.get_params = @get_params_wmf_dt1;
    
    output_folder_name = sprintf('benchmark_wmf_%s',dataset_name);
    benchmark_X_on_dt3(output_folder_name,X,params)
    params.preprocess = @(y) make_size_even(y);
    %%%evaluate_metrics_dt3(output_folder_name,params);

end



%% WMF in the logarithm domain

if WMF_log == 1

    X = @run_wmf_log;
    params.get_params = @get_params_wmf_dt1;

    output_folder_name = sprintf('benchmark_wmf_log_%s',dataset_name);
    benchmark_X_on_dt3(output_folder_name,X,params)
    params.preprocess = @(y) make_size_even(y);
    %%%evaluate_metrics_dt3(output_folder_name,params);

end




%% V-BM4D

if VBM4D == 1
    
    params = common_params;
    
    X = @run_vbm4d;
    
    output_folder_name = sprintf('benchmark_vbm4d_%s',dataset_name);  
    benchmark_X_on_dt3(output_folder_name,X,params)
    %%%evaluate_metrics_dt3(output_folder_name,params);

end



%% V-BM4D in the logarithm domain

if VBM4D_log == 1
    
    params = common_params;
    
    X = @run_vbm4d_log;
    
    output_folder_name = sprintf('benchmark_vbm4d_log_%s',dataset_name);
    benchmark_X_on_dt3(output_folder_name,X,params)
    %%%evaluate_metrics_dt3(output_folder_name,params);
end



%% Tensor DL

if TENSOR_DL == 1
    
    params = common_params;
    
    X = @run_tensor_dl;
    
    params.noise_estimator = @(y) estimate_noise_dt1_max(y);
    
    output_folder_name = sprintf('benchmark_tensor_dl_%s',dataset_name);
    benchmark_X_on_dt3(output_folder_name,X,params)
    d = 4;
    params.preprocess = @(y) crop_image(y,d);
    %%%evaluate_metrics_dt3(output_folder_name,params);
    
end



%% Tensor DL in the logarithm domain

if TENSOR_DL_log == 1
    
    params = common_params;  
    
    X = @run_tensor_dl_log;
    
    params.noise_estimator = @(y) estimate_noise_dt1_max(y);
    
    output_folder_name = sprintf('benchmark_tensor_dl_log_%s',dataset_name);
    benchmark_X_on_dt3(output_folder_name,X,params)
    d = 4;
    params.preprocess = @(y) crop_image(y,d);
    %%%evaluate_metrics_dt3(output_folder_name,params);

end



%% Original BM4D: BM4D with voxel-wise noise estimation

if BM4D_org == 1
    
    params = common_params;  
    
    X = @run_bm4d;
    
    output_folder_name = sprintf('benchmark_bm4d_%s',dataset_name);
    benchmark_X_on_dt3(output_folder_name,X,params)
    %%%evaluate_metrics_dt3(output_folder_name,params);
    
end



%% Original BM4D in the logarithm domain

if BM4D_org_log == 1
    
    params = common_params;
    
    X = @run_bm4d_log;
    
    output_folder_name = sprintf('benchmark_bm4d_log_%s',dataset_name);  
    benchmark_X_on_dt3(output_folder_name,X,params)
    %%%evaluate_metrics_dt3(output_folder_name,params);
    
end


%% BM4D with iid noise assumption

if BM4D_iid == 1
    
    params = common_params;
    
    X = @run_bm4d_iidnoise;
    
    params.noise_estimator = @(y) estimate_noise_dt1_max(y);
    
    output_folder_name = sprintf('benchmark_bm4d_iidnoise_%s',dataset_name);
    benchmark_X_on_dt3(output_folder_name,X,params)
    %%%evaluate_metrics_dt3(output_folder_name,params);
    
end


%% BM4D with iid noise assumption in the logarithm domain

if BM4D_iid_log == 1
    
    params = common_params;
    
    X = @run_bm4d_iidnoise_log;
    
    params.noise_estimator = @(y) estimate_noise_dt1_min(y);
    
    output_folder_name = sprintf('benchmark_bm4d_iidnoise_log_%s',dataset_name);
    benchmark_X_on_dt3(output_folder_name,X,params)
    %%%evaluate_metrics_dt3(output_folder_name,params);
    
end



%% MS BM4D (dwt3)

if MS_BM4D_dwt3 == 1
    for n_levels = 5     % number of decomposition level

        params = common_params;
        
        X = @run_bm4d_mix;
        params.method.n_levels = n_levels;
        params.method.wname = 'db7';
        
        params.noise_estimator = @(y) estimate_noise_dt1_max(y);
        
        output_folder_name = sprintf('benchmark_bm4d_mix_%s_%s_%d',...
            dataset_name, params.method.wname,n_levels);
        benchmark_X_on_dt3(output_folder_name,X,params)
    %%%evaluate_metrics_dt3(output_folder_name,params);

    end
end



%% MS BM4D (dwt3) in the logarithm domain

if MS_BM4D_dwt3_log == 1
    for n_levels = 5     % number of decomposition level

        params = common_params;
        
        X = @run_bm4d_mix_log;
        params.method.n_levels = n_levels;
        params.method.wname = 'db7';
        
        params.noise_estimator = @(y) estimate_noise_dt1_min(y);
        
        output_folder_name = sprintf('benchmark_bm4d_mix_log_%s_%s_%d',...
            dataset_name, params.method.wname,n_levels);
        benchmark_X_on_dt3(output_folder_name,X,params)
    %%%evaluate_metrics_dt3(output_folder_name,params);

    end
end



%% MS BM4D (dualtree3)

if MS_BM4D_dualtree3 == 1
    for n_levels = 3     % number of decomposition level

        params = common_params;
        
        X = @run_bm4d_mix_dualtree3;
        params.method.n_levels = n_levels;
        params.method.filter_bank = 'antonini';
        
        params.noise_estimator = @(y) estimate_noise_dt1_max(y);
        
        output_folder_name = sprintf('benchmark_bm4d_mix_dualtree3_%s_%s_%d',...
            dataset_name, params.method.filter_bank, n_levels);
        benchmark_X_on_dt3(output_folder_name,X,params)
    %%%evaluate_metrics_dt3(output_folder_name,params);

    end
end



%% MS BM4D (dualtree3) in the logarithm domain

if MS_BM4D_dualtree3_log == 1
    for n_levels = 3     % number of decomposition level

        params = common_params;
        
        X = @run_bm4d_mix_dualtree3_log_fast; %%%%% FAST VERSION %%%% LESS PADDING
        params.method.n_levels = n_levels;
        params.method.filter_bank = 'antonini';
        
        params.noise_estimator = @(y) estimate_noise_dt1_max(y); % estimate_noise_dt1_min
        
        output_folder_name = sprintf('benchmark_bm4d_mix_dualtree3_log_%s_%s_%d',...
            dataset_name, params.method.filter_bank, n_levels);
        
        benchmark_X_on_dt3(output_folder_name,X,params)
        %%%evaluate_metrics_dt1(output_folder_name,params);

    end
end



