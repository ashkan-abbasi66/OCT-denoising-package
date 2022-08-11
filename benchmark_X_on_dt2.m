function benchmark_X_on_dt2(output_folder_name,X,params)
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

test_indices = params.test_indices;
frame_numbers = params.frame_numbers;
n_frames = params.n_frames;

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
    for frame_number = frame_numbers(ii,:)
    
        img_number = test_indices(ii);

        % strnumber = num2str(img_number);

        fprintf('\nDenoising volume #%d ... \n',img_number)

        % **
        % Read noisy images
        % **

        load(fullfile(dataset_path,sprintf('%0.2d.mat',img_number)),'imn')

        max_frame_number = size(imn,3);
        [frame_range,main_frame] = get_frame_range(frame_number,n_frames,max_frame_number);

        if length(valid_rows) > 1
            noisy_imgs = imn(valid_rows, :,frame_range);
        else
            noisy_imgs = imn(:,:,frame_range);
        end



        % **
        % Denoising method
        % **

        output_filename = sprintf('%0.2d_%0.3d.tif',img_number,frame_number);
        output_path = fullfile(result_path,output_filename);


        % noise estimation and its runtime (rt)
        rt = 0; 
        sigma_value = -1;
        if isa(noise_estimator,'function_handle')
            tic
            sigma_value = noise_estimator(noisy_imgs);
            rt = toc;
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
        running_times(ii) = run_time;

        if size(denoised_imgs,3) > 1
            % The method outputs a volume
            denoised_imgs = denoised_imgs(:,:,1:n_frames);
            im_out = denoised_imgs(:,:,main_frame);
        else
            % The method outputs an image
            im_out = denoised_imgs;
        end

        % **
        % Save just one image(slice) from the output volume
        % **

        imwrite(uint8(im_out),output_path,'tif');

        % save the whole output
        if save_mat == true
            mat_file_name = sprintf('denoised_volume_%0.2d_%0.3d.mat',...
                img_number,frame_number);
            save(fullfile(result_path,mat_file_name),...
                'denoised_imgs','run_time')
        end
    end
end

fprintf('Average of running times = %g\n',mean(running_times))

diary off

end

% ========================================================================
% ========================================================================

function check_required_inputs(output_folder_name,X,params)

    if isempty(output_folder_name)
        msg = 'ERROR: Name of the output folder is required.';
        error(msg)
    end
    
    if ~isa(X,'function_handle')
        msg = 'ERROR: A function handle to a denoising method is required.';
        error(msg)
    end
    
    if ~exist('params','var') || isempty(params)
        msg = 'ERROR: A structure containing the input parameters is required.';
        error(msg)
    end

    % indices of test volumes for denoising
    if isfield(params,'test_indices')
        test_indices = params.test_indices;
        if ~all([test_indices >= 1, test_indices <= 13])
            msg = ['ERROR: This dataset has 13 volumes. \n',...
                '""test_indices"" is outside the valid range.'];
            error(sprintf(msg))
        end
    else
        error("""test_indices"" is required.")
    end

    if ~isfield(params,'frame_numbers')
        error("ERROR: ""frame_numbers"" is required.")
    end
    
    msg = ['ERROR: number of rows in ""params.frame_numbers"" must be ' ...
    'equal to the number of indices in ""params.test_indices""'];
    assert(size(params.frame_numbers,1) == length(params.test_indices),msg)
    
    if ~isfield(params,'n_frames')
        error("ERROR: ""n_frames"" is required.")
    end
    
end


function [frame_range,main_frame] = get_frame_range(frame_number,n_frames,max_frame_number)
% 
% Returns a valid range of frame numbers for 3D denoising
% 
    if any([frame_number < 1, frame_number > max_frame_number])
        msg = ['ERROR: The loaded volume has ', ...
            num2str(size(imn,3)) ' frames.\n',...
            'Input frame range is not valid.'];
        error(sprintf(msg))
    end
    
    if rem(n_frames,2) ~= 1
            error('ERROR: Number of frames must be odd')
    end
    nf = (n_frames - 1) / 2;

    fprintf('Input frame range is:\n')
    frame_range = frame_number - nf : frame_number + nf;
    disp(frame_range)

    % All frame numbers must be greater than or equal to 1
    modified = false;    
    if ~all(frame_range >= 1) 
        % number of frame numbers violating the condition
        num = nnz(frame_range < frame_number); 

        % use the range of valid frame numbers
        tmp = frame_range(frame_range > frame_number);
        tmp = tmp(end:-1:1);

        % modify the range of valid numbers
        frame_range(frame_range < frame_number) = tmp(1:num);
        
        modified = true;
    end
    
    % All frame numbers must be less than or equal to "max_frame_number"
    if ~all(frame_range <= max_frame_number) 
        % number of frame numbers violating the condition
        num = nnz(frame_range > max_frame_number); 
        
        % use the range of valid frame numbers
        tmp = frame_range(frame_range < max_frame_number);
        tmp = tmp(end:-1:1);
        
        % modify the range of valid numbers
        frame_range(frame_range > max_frame_number) = tmp(1:num);
        
        modified = true;
    end
    
    if modified == true
        msg = 'Input frame range is not valid. \nIt is modified as follows:';
        warning(sprintf(msg))
        disp(frame_range)
    end
    
    main_frame = nf + 1;

end