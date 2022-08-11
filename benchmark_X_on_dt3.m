function benchmark_X_on_dt3(output_folder_name,X,params)
% INPUT:
%   output_folder_name: a folder with the given name is created in 
%     the Results folder.
% 
%   X: handle to a denoising function (e.g., @run_bm4d)
% 
%   params: parameters used by this function and the method
% 
% OUTPUT:
%   the outputs are saved in the Result folder
% 
% 
% 
%% Check and get inputs

check_required_inputs(output_folder_name,X,params);

start_frame = params.start_frame;
n_frames = params.n_frames;
batch_size = params.batch_size;


if isfield(params,'save_mat')
    save_mat = params.save_mat;
else
    save_mat = false;
end

if isfield(params,'noise_estimator')
    noise_estimator = params.noise_estimator;
else
    noise_estimator = -1;
end


%% 

dataset_path = './Datasets/dt3_Farsiu_Ophthalmology_2013/';
result_path = fullfile('./Results/',output_folder_name);
% command line outputs will be saved here:
diary_file_name = [output_folder_name,'.txt'];

if ~exist(result_path,'dir')
    mkdir(result_path);
end

noisy_volumes = dir(fullfile(dataset_path, '*.mat'));

N_volumes=length(noisy_volumes);% number of test images

running_times = zeros(N_volumes,1);

diary_file_path = fullfile(result_path,diary_file_name);
if exist(diary_file_path, 'file')
    delete(diary_file_path); 
end
diary(diary_file_path)



%%

for ii=1:N_volumes
    
    if N_volumes == 1
        filename = noisy_volumes.name;
    else
        filename = noisy_volumes(ii).name;
    end
    
    fprintf('\nDenoising volume #%d ... \n',ii)
    fprintf('   File name: %s \n',filename)
    
    load(fullfile(dataset_path,filename));
    % "images" contain a noisy volume with size 512x1000x100
    
    % **
    % Read noisy images
    % **
    
    %main_frame = 1; %%%%%%%
    
    left_frame = start_frame;
    right_frame = start_frame + n_frames - 1;
    
    noisy_imgs  = images(:,:,left_frame:right_frame);
    
    noisy_imgs = make_size_even(noisy_imgs);
    
    % **
    % Denoising method
    % **
    
    denoised_imgs = noisy_imgs;
    
    rt = 0; 
    if batch_size > 0
        for jj = 1: batch_size: size(noisy_imgs,3)
            last_index = jj + batch_size - 1;
            fprintf('\n\n----> Frame %d to %d is processing now ... \n', ...
                jj, last_index)
            noisy_batch = noisy_imgs(:,:, jj:last_index);
            [denoised_batch,run_time] = denoise_volume(noise_estimator, ...
                noisy_batch, params, X);
            rt = rt + run_time;
            
            denoised_batch = denoised_batch(:,:,1:batch_size);
            
            denoised_imgs(:,:, jj: last_index) = denoised_batch;          
        end
    else
        [denoised_imgs,run_time] = denoise_volume(noise_estimator, ...
            noisy_imgs, params, X);
        rt = run_time;
    end
    
    running_times(ii) = rt;
    
    if size(denoised_imgs,3) > 1
        denoised_imgs = denoised_imgs(:,:,1:n_frames);
    end
    
    % save the whole output
    if save_mat == true
        [~, fn, ~] = fileparts(filename);
        if batch_size>0
            mat_file_name = sprintf('%s_batchwise_bs%d.mat',fn, batch_size);
        else
            mat_file_name = sprintf('%s_frames_%d_%d.mat',fn, left_frame, right_frame);
        end
        save(fullfile(result_path,mat_file_name),...
            'denoised_imgs','run_time')
    end
    
end

fprintf('Average of running times = %g\n',mean(running_times))

diary off

end

function [denoised_imgs,run_time] = denoise_volume(noise_estimator, noisy_imgs, params, X)
    rt = 0;
    sigma_value = -1;
    if isa(noise_estimator,'function_handle')
        tic
        sigma_value = noise_estimator(noisy_imgs);
        rt = toc;
    end
    
    if isfield(params,'get_params')
        get_params = params.get_params;
    else
        get_params = -1;
    end
    
    if isa(get_params,'function_handle')
        params.method = get_params(noisy_imgs,sigma_value);
    else
        params.method.sigma_value = sigma_value;
    end
    
    if sigma_value ~= -1 
        fprintf('Applied noise level: %4.0f \n',sigma_value);
    end

    [denoised_imgs,run_time] = X(noisy_imgs,params.method);
    run_time = run_time + rt;
    
end


function check_required_inputs(output_folder_name,X,params)

    if isempty(output_folder_name)
        msg = 'Name of the output folder is required.';
        error(msg)
    end
    
    if ~isa(X,'function_handle')
        msg = 'A function handle to a denoising method is required.';
        error(msg)
    end
    
    if ~exist('params','var') || isempty(params)
        msg = 'A structure containing the input parameters is required.';
        error(msg)
    end

%     % indices of test volumes for denoising
%     if isfield(params,'test_indices')
%         test_indices = params.test_indices;
%         if ~all([test_indices >= 1, test_indices <= 18])
%             msg = ['This dataset has 18 volumes. \n',...
%                 '""test_indices"" is outside the valid range.'];
%             error(sprintf(msg))
%         end
%     else
%         error("""test_indices"" is required.")
%     end
%     
    if ~isfield(params,'n_frames')
        error("ERROR: ""n_frames"" is required.")
    end
    
    if ~all([params.n_frames >1, params.n_frames <=100])
        msg = ['ERROR: For this dataset, ""n_frames"" must be a '...
            'number in the range [2,100].'];
        error(msg)
    end
end