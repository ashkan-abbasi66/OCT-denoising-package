function [input_img,filename] = get_data()
% Assumption: 
%   There are two input types:
%       .MAT file containing "imn" and "im"
%       .tif file
% 

    [file,pth,indx] = uigetfile( ...
    {'*.tif';'*.mat'}, ...
       'Input MAT-file/Image');


    [filepath,name,ext] = fileparts(file);

    % is input file a MAT-file?
    if strcmpi(ext,'.mat')
        
        load(fullfile(pth,file))
        
        if exist('imn','var')
            input_img = imn;
        elseif exist('im','var')
            input_img = im;
        else
            error('The MAT-file does not contain "imn" or "im".');
        end
        
        nFrames = size(input_img,3);
        if nFrames>1
            msg = sprintf('Select one image number (between 1 and %d): ', nFrames);
            slice_num = str2double(inputdlg(msg));
        end
        
        input_img = input_img(:,:,slice_num);
        
    else
        % the input file is an image.
        input_img = double(imread(fullfile(pth,file)));
        
    end
    
    filename = file;
end