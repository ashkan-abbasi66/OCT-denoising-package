%% 

clear
clc

addpath('./utils/')

% Get an image to place ROIs over it
[img,img_fname] = get_data();

% Get ROI file
[file,pth,indx] = uigetfile( ...
{'*.mat'}, ...
   sprintf('ROI for %s',img_fname));
load(fullfile(pth,file))


%% Show ROIs over an image

opts.fname = get_file_name();
% opts.fname = 'im_out_3dsbsdi';
opts.fpath = './test_folder';
% opts.LineWidth = 5;
% opts.Color = 'yellow';
opts.LineWidth = 1%4;%4 - 6
opts.Color = 'red';

put_ROIs(img,pos,opts);


opts.mag_factor = 1;
show_mag_ROIs(img,pos,opts)



%% ============================================
%               AUXILIARY FUNCTIONS
% =============================================

% ----------------------------------------------



% ----------------------------------------------
function imgRGB = put_ROIs(input_img,pos,opts)
    
    % Prepare Inputs
    
     input_img = mat2gray(input_img);
%     input_img = (input_img)/255;
    
    if exist('opts','var')
        if isfield(opts,'fname') && isfield(opts,'fpath')
            fname = opts.fname;
            fpath = opts.fpath;
        else
            fname = [];
        end
        
        if isfield(opts,'Color')
            Color = opts.Color;
        else
            Color = 'red';
        end
        
        if isfield(opts,'LineWidth')
            LineWidth = opts.LineWidth;
        else
            LineWidth = 1;
        end
    end
    
    
    % Create an RGB image with ROIs over it.
    [r,g,~]=size(input_img);
    rgb=zeros(r,g,3); 
    rgb(:,:,1)=input_img;
    rgb(:,:,2)=rgb(:,:,1);
    rgb(:,:,3)=rgb(:,:,1);
    imgRGB=rgb; 

    % insert ROIs
    insertRect = @(im,pos) insertShape(im,...
        'Rectangle',pos,...
        'LineWidth',LineWidth,...
        'Color',Color);
    
    if ~iscell(pos)
        pos = mat2cell(pos,1);
    end
    
    for j = 1:length(pos)
        imgRGB = insertRect(imgRGB,pos{j});
    end
    
    
    % show and [save it]
    figure,imshow(imgRGB)
    
    if ~isempty(fname)
        
        if ~exist('fpath','dir')
            mkdir(fpath);
        end
        
        imwrite(imgRGB,fullfile(fpath,[fname,'.tif']));
    end
    
end


% ----------------------------------------------
function show_mag_ROIs(input_img,pos,opts)
% Show magnified ROIs

    input_img = mat2gray(input_img);
%     input_img = input_img/255;
    
    mag_factor = 2.5;
    if exist('opts','var')
        if isfield('opts','mag_factor')
            mag_factor = opts.mag_factor;
        end
    end
    
    if ~iscell(pos)
        pos = mat2cell(pos,1);
    end
    
    for j=1:length(pos)
        cropped_im=imcrop(input_img,pos{j});
        mag_img = imresize(cropped_im,mag_factor);
        figure,imshow(mag_img)
        title(sprintf('ROI #%d',j));
        
        if exist('opts','var') && ~isempty(opts)

            fname = opts.fname;
            fpath = opts.fpath;

            if ~exist('fpath','dir')
                mkdir(fpath);
            end

            imwrite(mag_img,fullfile(fpath,...
                sprintf([fname,'_%d.tif'],j)));
        end
    end
end


% -------------------------------------------------
function fname = get_file_name()
    prompt = {'output file name:'};
    dlgtitle = 'get file name';
    dims = [1 50];
    definput = {'im_test'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    fname = answer{1};
    if isempty(fname)
        fname = definput{1};
    end
end