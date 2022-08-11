%{

IMPORTANT NOTE ABOUT NOISE DISTRIBUTION IN THIS DATASET (dt2)

We have provided two functions for noise estimation on images from this
dataset:
  estimate_noise_dt2_max
  estimate_noise_dt2_min

Both of these methods CANNOT correctly estimate the noise level for this
dataset. Therefore, researchers usually increase the estimated noise
level by adding or multipling constants. They experimentally and viauslly
select an approporiate noise level. 

Each volume has 128 frames (slices). 

Assume that you want to use 8 frames to denoise a selected frame. 
E.g., to denoise frame number 021, a subset of frames with numbers between
21-4 to 21-4+8 is used.

We use the following NOISE_FACTORs:

To denoise frame number 021, we set the NOISE_FACTOR to 1.25.
To denoise frame number 069, set NOISE_FACTOR to 1.5.
To denoise frame number 109, set NOISE_FACTOR to 1.75.

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
params.test_indices = 1:13; % [3,5]
params.frame_numbers = ones(length(params.test_indices),1)*[21,69, 109]; % [10, 60] for all images
% params.frame_numbers = ones(length(params.test_indices),1)*[13,21,29,37,45,53,61,77];
% ----------------------------------------------------------

% params.valid_rows = 150:512 - 1; % 150:380 - 1 -- TEST 100:
params.valid_rows = 51:562;

dataset_name = 'dt2';

% ===========================================
% ===========================================
% Examples:
%   Denoise the whole volume at once: (stores selected frames)
%       start_frame = 1;
%       n_frames = 128;
%       batch_size = -1;
%   Denoise the whole volume in batches: (stores selected frames)
%       start_frame = 1;
%       n_frames = 128;
%       batch_size = 8;
%   Denoise part of a volume: (stores middle frame)
%       start_frame = 50;
%       n_frames = 55;

%       batch_size = -1;



%% ===========================================

DO_DENOISE = true;

% TEST ALL LOOP -- REMOVE THIS LOOP %%%%%%%%%%%
for fr = params.frame_numbers(1,:)
    
    % TO DENOISE SELECTED FRAMES BY CONSIDERING A SUBSET OF FRAMES AROUND
    % THEM
    batchsize = 4;
    params.start_frame = fr - batchsize/2;%1;
    params.n_frames = batchsize;%128; % >= 2 
    params.batch_size = -1;
    
%     params.start_frame = 1;
%     params.n_frames = 128;
%     params.batch_size = -1;
    
    noise_factors = {1.25};
    % Frame number 21 => noise_factor = 1
    
params.save_mat = true;

% Directly used by "benchmark_X_on_..." and "evaluate_metrics_..."
common_params = params;



%% Select a denoising method

% To run a method, set its corresponding variable to 1
BM4D_iid = 1; %%%%
BM4D_iid_log = 1; %%%%
MS_BM4D_dwt3_log = 1; %%%% 
MS_BM4D_dualtree3_log = 1; %%%%



%% BM4D with iid noise assumption -- %%%%


if BM4D_iid == 1
    for coeff = noise_factors
        
        params = common_params;

        X = @run_bm4d_iidnoise;

%         params.noise_estimator = @(y) estimate_noise_dt2_max(y)*2 - 1;
        params.noise_estimator = @(y) estimate_noise_dt2_max(y)*coeff{1};

        output_folder_name = sprintf('benchmark_bm4d_iidnoise_%s_times%g', ...
            dataset_name, coeff{1});
        
        if DO_DENOISE
            benchmark_X_on_dt2_new(output_folder_name,X,params)
        else
            evaluate_metrics_dt2(output_folder_name,params); 
        end
        
    end
end


%% BM4D with iid noise assumption in the logarithm domain -- %%%%


if BM4D_iid_log == 1
    for coeff = noise_factors
        
        params = common_params;

        X = @run_bm4d_iidnoise_log;

%         params.noise_estimator = @(y) estimate_noise_dt2_max(y)*2 - 1;
        params.noise_estimator = @(y) estimate_noise_dt2_max(y)*coeff{1};
        
        output_folder_name = sprintf('benchmark_bm4d_iidnoise_log_%s_times%g', ...
            dataset_name, coeff{1});
        
        if DO_DENOISE
            benchmark_X_on_dt2_new(output_folder_name,X,params)
        else
            evaluate_metrics_dt2(output_folder_name,params); 
        end
        
    end
end



%% MS BM4D (dwt3) in the logarithm domain -- %%%%
% for this dataset, it seems that the log version is not better than the
% spatial domain.
% 

wnames = {'sym4'};

if MS_BM4D_dwt3_log == 1
    for coeff = noise_factors
        for n_levels = 1     % number of decomposition level

            params = common_params;

            X = @run_bm4d_mix_log;
            params.method.n_levels = n_levels;
            params.method.wname = wnames{1};

%             params.noise_estimator = @(y) estimate_noise_dt2_max(y)*2.3;
            params.noise_estimator = @(y) estimate_noise_dt2_max(y)*coeff{1};

            output_folder_name = sprintf('benchmark_bm4d_mix_log_%s_%s_%d_times%g',...
                dataset_name, params.method.wname, n_levels, coeff{1});

        if DO_DENOISE
            benchmark_X_on_dt2_new(output_folder_name,X,params)
        else
            evaluate_metrics_dt2(output_folder_name,params); 
        end

        end
    end
end



%% MS BM4D (dualtree3) in the logarithm domain -- %%%%
% Allowed filter banks are:
% ("nearsym5_7"), 'nearsym13_19', 'antonini', or 'legall'

% noise_factors = {1, 1.25, 1.5};

if MS_BM4D_dualtree3_log == 1
    for coeff = noise_factors
        for n_levels = 2     % number of decomposition level

            params = common_params;

            X = @run_bm4d_mix_dualtree3_log; % run_bm4d_mix_dualtree3_log_fast
            params.method.n_levels = n_levels;
            params.method.filter_bank = 'nearsym5_7'; 

    %         params.noise_estimator = @(y) estimate_noise_dt2_max(y)*2.8;
%             params.noise_estimator = @(y) estimate_noise_dt2_max(y)*2 - 1;
            params.noise_estimator = @(y) estimate_noise_dt2_max(y)*coeff{1};

            output_folder_name = sprintf('benchmark_bm4d_mix_dualtree3_log_%s_%s_%d_times%g',...
                dataset_name, params.method.filter_bank, n_levels, coeff{1});
            
            if DO_DENOISE
                benchmark_X_on_dt2_new(output_folder_name,X,params)
            else
                evaluate_metrics_dt2(output_folder_name,params); 
            end

        end
    end
end


end
