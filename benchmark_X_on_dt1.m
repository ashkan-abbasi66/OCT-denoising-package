function benchmark_X_on_dt1(output_folder_name,X,params)
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
n_frames = params.n_frames;

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

dataset_path = './Datasets/dt1_Bioptigen_SDOCT/';
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

for ii=1:N_images
    
    img_number = test_indices(ii);
    
    strnumber = num2str(img_number);
    
    fprintf('\nDenoising volume #%d ... \n',img_number)
    
    img_path = fullfile(dataset_path,strnumber);
      
    % **
    % Read noisy images
    % **
    
    test_img = double(imread(fullfile(img_path,'test.tif')));
    
    noisy_imgs = zeros(size(test_img,1),size(test_img,2),5);
    
    noisy_imgs(:,:,1) = test_img; % main_frame
    noisy_imgs(:,:,2) = double(imread(fullfile(img_path,'1.tif')));
    noisy_imgs(:,:,3) = double(imread(fullfile(img_path,'2.tif')));
    noisy_imgs(:,:,4) = double(imread(fullfile(img_path,'3.tif')));
    noisy_imgs(:,:,5) = double(imread(fullfile(img_path,'4.tif')));
    
    main_frame = 1;


    noisy_imgs  = noisy_imgs(:,:,1:n_frames);
    
    noisy_imgs = make_size_even(noisy_imgs);
    
    % **
    % Denoising method
    % **
    
    output_filename = sprintf('%0.2d.tif',img_number);
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
        denoised_imgs = denoised_imgs(:,:,1:n_frames);
        %im_out = H_weighted_averaging_fusion(denoised_imgs,400);
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
        mat_file_name = sprintf('denoised_volume_%0.2d.mat',img_number);
        save(fullfile(result_path,mat_file_name),...
            'denoised_imgs','run_time')
    end
    
end

fprintf('Average of running times = %g\n',mean(running_times))

diary off

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

    % indices of test volumes for denoising
    if isfield(params,'test_indices')
        test_indices = params.test_indices;
        if ~all([test_indices >= 1, test_indices <= 18])
            msg = ['This dataset has 18 volumes. \n',...
                '""test_indices"" is outside the valid range.'];
            error(sprintf(msg))
        end
    else
        error("""test_indices"" is required.")
    end
    
    if ~isfield(params,'n_frames')
        error("ERROR: ""n_frames"" is required.")
    end
    
    if ~all([params.n_frames >1, params.n_frames <=5])
        msg = ['ERROR: For this dataset, ""n_frames"" must be a '...
            'number in the range [2,5].'];
        error(msg)
    end
end