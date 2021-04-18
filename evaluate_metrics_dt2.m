function evaluate_metrics_dt2(output_folder_name,params)
% 
% INPUT:
%   OUTPUT_FOLDER_NAME: name of the folder where the denoised images are
%   stored.
%   TEST_INDICES: indices of the images.E.g., test_indices = 1:18
%   PREPROCESS: function handle to run over the image before further
%   processing.
%   
% OUTPUT:
%   computed metric values are stored in an Excel file into the output
%   folder. Also, a table is shown to the user.
% 
% USAGE:
%   benchmark_bm4d_dt1
%   evaluate_metrics_dt1("benchmark_bm4d_dt1",1)
% 
% 


test_indices = params.test_indices;
frame_numbers = params.frame_numbers;

assert(size(params.frame_numbers,1) == length(test_indices))

if isfield(params,'preprocess')
    preprocess = params.preprocess;
else
    preprocess = -1;
end


% **
% Settings
% ** 

addpath('./Metrics')

dataset_path = './Datasets/dt2_topcon_oct1000_seg_normal/';

% output images will be saved here:
result_path = fullfile('./Results/',output_folder_name);

rois_path = './Metrics/rois_for_dt2/';

background_indices = [1]; % indices of ROIs which are background regions

N_images = numel(frame_numbers); % Number of output images

column_names={'MSR','CNR','ENL','TP','EP'};
N_metrics = length(column_names); 


% pre-allocating memory for creating a table of measures.
metric_values=zeros(N_images,N_metrics);% #images by #quality measures & time
row_names=cell(N_images,1);


% **
% loop over all outputs/results
% **

kk = 1;

whandle = waitbar(0,'Please wait...');
for ii = 1:length(test_indices)
    waitbar(ii / N_images)
    
    img_number = test_indices(ii);
    
    for jj =  1:length(frame_numbers)
        
        frame_number = frame_numbers(ii,jj);
        
        
        % **
        % Read data
        % **

        % read noisy image
        load(fullfile(dataset_path,sprintf('%0.2d.mat',img_number)),'imn')
        imn = imn(150:512 - 1,:,frame_number);


        % read the output/denoised image - %%%%%%%%%%%%%%%%%%%%%
        output_filename = sprintf('%0.2d_%0.3d.tif',img_number,frame_number);
        try
            im_out_path = fullfile(result_path,output_filename);
            im_out = double(imread(im_out_path));
        catch
            msg = 'The file does not exist:';
            close(whandle)
            error('INPUT ERROR: \n %s\n%s',msg,im_out_path)
        end

        % load ROIs
        posfname=[rois_path sprintf('%0.2d_%0.3d',img_number,frame_number)];
        load(posfname,'pos');

        [roi,~] = get_roi_pos(im_out,pos);


        if preprocess ~= -1    
            im_out = preprocess(im_out);
        end

        % **
        % Compute Metrics
        % **

        % [PSNR,SSIM]=comp_psnr_ssim(Truth ,im_out);

        [MSR,CNR,ENL] = comp_MSR_CNR_ENL(roi,background_indices);

        TP = comp_TP(im_out,imn,pos,background_indices);
        EP = comp_EP(im_out,imn,pos,background_indices);

        metric_values(kk,:) = [MSR,CNR,ENL,TP,EP];


        row_names{kk} = output_filename;
        
        kk = kk + 1;
        
    end
    
end
close(whandle) 

% Add the average row
dim = 1;
metric_values(end+1,:) = mean(metric_values(1:N_images,:),dim);
row_names{end+1} = 'Avg.';


% **
% Save the metric results in a table and store it in a excel file
% **

assert(size(metric_values,2) == (length(column_names)),...
    "metric results does not compatible with the number of columns")


table_results = array2table(metric_values);
table_results.Properties.VariableNames = column_names;
table_results.Properties.RowNames = row_names;
table_results.Properties.DimensionNames{1} = 'ImageNumber';

excel_file_path = fullfile(result_path,strcat(output_folder_name,'.xlsx'));
writetable(table_results,excel_file_path,...
    'Sheet',1,'Range','B2',...
    'WriteRowNames',true,...
    'FileType','spreadsheet')

fprintf('\nThe numerical results are stored in:\n%s\n',excel_file_path);

% **
% draw a table
% **

fh=figure;
columns_formats=repmat({'short g'},1,numel(column_names));

ui_table_obj = ...
    uitable('Parent', fh,...
    'Data', metric_values,...
    'RowName',row_names,...
    'ColumnName', column_names,...
    'ColumnFormat', columns_formats,...
    'Units','normalized','Position',[0 0 1 1]);
ui_table_obj.FontSize=12;
fh.Name=[output_folder_name];% figure title