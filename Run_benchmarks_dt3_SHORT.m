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

% params.start_frame = 1;
% params.n_frames = 100; % >= 2 
% params.batch_size = 10;
% ===========================================
% ===========================================




%% Select a denoising method

% To run a method, set its corresponding variable to 1
BM4D_iid = 0; %%%
BM4D_iid_log = 0;
MS_BM4D_dwt3_log = 0; %% 
MS_BM4D_dualtree3_log = 1; %%%

%% ===========================================

DO_DENOISE = false;

% used only for evaluation
params.file_type = 'frames'; % "frames" or "batchwise"
params.files_pattern = sprintf('*%s*.mat*', params.file_type);
params.frame_numbers = [15, 45, 55, 65, 85];

params.save_mat = true;

% Directly used by "benchmark_X_on_..." and "evaluate_metrics_..."
common_params = params;



%% BM4D with iid noise assumption

if BM4D_iid == 1
    
    params = common_params;
    
    X = @run_bm4d_iidnoise;
    
    params.noise_estimator = @(y) estimate_noise_dt1_max(y);
    
    output_folder_name = sprintf('benchmark_bm4d_iidnoise_%s',dataset_name);
    
    if DO_DENOISE
        benchmark_X_on_dt3(output_folder_name,X,params)
    else
        evaluate_metrics_dt3('benchmark_bm4d_iidnoise_dt3_SELECTED_IMAGES_FOR_EVAL',params);
    end
    
end


%% BM4D with iid noise assumption in the logarithm domain

if BM4D_iid_log == 1
    
    params = common_params;
    
    X = @run_bm4d_iidnoise_log;
    
    params.noise_estimator = @(y) estimate_noise_dt1_min(y);
    
    output_folder_name = sprintf('benchmark_bm4d_iidnoise_log_%s',dataset_name);
    
    if DO_DENOISE
        benchmark_X_on_dt3(output_folder_name,X,params)
    else
        evaluate_metrics_dt3('benchmark_bm4d_iidnoise_log_dt3_SELECTED_IMAGES_FOR_EVAL',params);
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

        if DO_DENOISE
            benchmark_X_on_dt3(output_folder_name,X,params)
        else
            evaluate_metrics_dt3('benchmark_bm4d_mix_log_dt3_db7_5_IMAGES',params);
        end

    end
end



%% MS BM4D (dualtree3) in the logarithm domain

if MS_BM4D_dualtree3_log == 1
    for n_levels = 3     % number of decomposition level

        params = common_params;
        
%         X = @run_bm4d_mix_dualtree3_log_fast; %%%%% FAST VERSION %%%% LESS PADDING
        X = @run_bm4d_mix_dualtree3_log; %%%%% FAST VERSION %%%% LESS PADDING
        params.method.n_levels = n_levels;
        params.method.filter_bank = 'nearsym5_7';
%         params.method.filter_bank = 'antonini';
        
        params.noise_estimator = @(y) estimate_noise_dt1_max(y); % estimate_noise_dt1_min
        
        output_folder_name = sprintf('benchmark_bm4d_mix_dualtree3_log_%s_%s_%d',...
            dataset_name, params.method.filter_bank, n_levels);
        
        if DO_DENOISE
            benchmark_X_on_dt3(output_folder_name,X,params)
        else
            evaluate_metrics_dt3('benchmark_bm4d_mix_dualtree3_log_dt3_nearsym5_7_3_IMAGES',params);
        end

    end
end

