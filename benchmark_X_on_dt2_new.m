function benchmark_X_on_dt2_new(output_folder_name,X,params)
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

test_indices = params.test_indices;

% matrix of frame indices to be saved in separate TIF files.
frame_numbers = params.frame_numbers; 


if isfield(params,'valid_rows')
    valid_rows = params.valid_rows;
else
    valid_rows = -1;
end


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

if isfield(params,'get_params')
    get_params = params.get_params;
else
    get_params = -1;
end



%% 

dataset_path = './Datasets/dt2_topcon_oct1000_seg_normal/';
result_path = fullfile('./Results/',output_folder_name);
% command line outputs will be saved here:
diary_file_name = [output_folder_name,'.txt'];

if ~exist(result_path,'dir')
    mkdir(result_path);
end


N_images=length(test_indices);% number of test images

running_times = zeros(N_images,1);


diary_file_path = fullfile(result_path,diary_file_name);
if exist(diary_file_path, 'file')
    delete(diary_file_path); 
end
diary(diary_file_path)




%%
for ii= 1:N_images  
    
    vol_number = test_indices(ii);
    
    left_frame = start_frame;
    right_frame = start_frame + n_frames - 1;
    
    filename = sprintf('%0.2d.mat',vol_number);
    
    fprintf('\nDenoising volume #%d from %d to %d frames... \n',...
        vol_number, left_frame, right_frame)
    fprintf('   Filename: %s\n', filename)

    % **
    % Read noisy images
    % **
    
    load(fullfile(dataset_path,filename),'imn')
    
    
    if length(valid_rows) > 1
        noisy_imgs = imn(valid_rows, :,left_frame:right_frame);
    else
        noisy_imgs = imn(:, :,left_frame:right_frame);
    end
    
    noisy_imgs = make_size_even(noisy_imgs);
    
    % **
    % Denoising method
    % **  

    
%     denoised_imgs = noisy_imgs;
    
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
            mat_file_name = sprintf('%s_batchwise.mat',fn);
        else
            mat_file_name = sprintf('%s_frames_%d_%d.mat',fn, left_frame, right_frame);
        end
        save(fullfile(result_path,mat_file_name),...
            'denoised_imgs','run_time')
    end
    
    save_selected_frames(result_path, denoised_imgs, frame_numbers(ii,:), ...
    batch_size, left_frame, n_frames, vol_number);

    
end

fprintf('Average of SUM of the running times = %g\n',mean(running_times))
if batch_size >0
    num_batches = size(noisy_imgs,3) / batch_size;
    fprintf('Average of the running times for each BATCH = %g\n',...
        mean(running_times)/num_batches)
end

diary off

end



%% ========================================================================

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

    % indices of test volumes for denoising
    if isfield(params,'test_indices')
        test_indices = params.test_indices;
        if ~all([test_indices >= 1, test_indices <= 13])
            msg = ['This dataset has 13 volumes. \n',...
                '""test_indices"" is outside the valid range.'];
            error(sprintf(msg))
        end
    else
        error("""test_indices"" is required.")
    end
    
    if ~isfield(params,'n_frames')
        error("ERROR: ""n_frames"" is required.")
    end
    
    if ~all([params.n_frames >1, params.n_frames <=128])
        msg = ['ERROR: For this dataset, ""n_frames"" must be a '...
            'number in the range [2,128].'];
        error(msg)
    end
end



%% =======================================================================

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


%%


    
function save_selected_frames(result_path, denoised_imgs, frame_numbers, ...
    batch_size, left_frame, n_frames, vol_number)
    % 
    % save selected frames from the volume as separate TIF files
    %
    % Set "frame_numbers" to frame_numbers(ii,:)
    
    TOTAL_FRAMES = 128;
    
    
    if batch_size > 0 || n_frames == TOTAL_FRAMES
        inds = frame_numbers; % indices to be saved as separate TIF files
    else
        inds = floor(size(denoised_imgs, 3)/2);
    end
    
    for jj = 1: length(inds)
        if batch_size > 0 || n_frames == TOTAL_FRAMES
            % the whole volume is denoised in a batch-wise manner
            output_filename = sprintf('%0.2d_%0.3d.tif',vol_number,inds(jj));
            im_out = denoised_imgs(:,:,inds(jj));
        else
            % it assumes that you want to store the middle frame
            output_filename = sprintf('%0.2d_%0.3d.tif',...
                vol_number,left_frame + inds);
            im_out = denoised_imgs(:,:,inds + 1);
        end
        
        output_path = fullfile(result_path,output_filename);
        imwrite(uint8(im_out),output_path,'tif');
        
        fprintf('--- "%s" is saved.\n', output_filename)
    end
end

% function [frame_range,main_frame] = get_frame_range(frame_number,n_frames,max_frame_number)
% % 
% % Returns a valid range of frame numbers for 3D denoising
% % 
%     if any([frame_number < 1, frame_number > max_frame_number])
%         msg = ['ERROR: The loaded volume has ', ...
%             num2str(size(imn,3)) ' frames.\n',...
%             'Input frame range is not valid.'];
%         error(sprintf(msg))
%     end
%     
%     if rem(n_frames,2) ~= 1
%             error('ERROR: Number of frames must be odd')
%     end
%     nf = (n_frames - 1) / 2;
% 
%     fprintf('Input frame range is:\n')
%     frame_range = frame_number - nf : frame_number + nf;
%     disp(frame_range)
% 
%     % All frame numbers must be greater than or equal to 1
%     modified = false;    
%     if ~all(frame_range >= 1) 
%         % number of frame numbers violating the condition
%         num = nnz(frame_range < frame_number); 
% 
%         % use the range of valid frame numbers
%         tmp = frame_range(frame_range > frame_number);
%         tmp = tmp(end:-1:1);
% 
%         % modify the range of valid numbers
%         frame_range(frame_range < frame_number) = tmp(1:num);
%         
%         modified = true;
%     end
%     
%     % All frame numbers must be less than or equal to "max_frame_number"
%     if ~all(frame_range <= max_frame_number) 
%         % number of frame numbers violating the condition
%         num = nnz(frame_range > max_frame_number); 
%         
%         % use the range of valid frame numbers
%         tmp = frame_range(frame_range < max_frame_number);
%         tmp = tmp(end:-1:1);
%         
%         % modify the range of valid numbers
%         frame_range(frame_range > max_frame_number) = tmp(1:num);
%         
%         modified = true;
%     end
%     
%     if modified == true
%         msg = 'Input frame range is not valid. \nIt is modified as follows:';
%         warning(sprintf(msg))
%         disp(frame_range)
%     end
%     
%     main_frame = nf + 1;
% 
% end