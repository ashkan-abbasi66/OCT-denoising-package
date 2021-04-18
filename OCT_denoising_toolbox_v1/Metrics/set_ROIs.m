% 
% 
% For MAT-files, it is assumed that the image or images are stored
% in a variable named "im" or "imn".
% 

clear
clc

addpath('./utils/')

% Get an image to place ROIs over it
[img,img_fname] = get_data();
img = img(150:380,:,:);
[~,img_fname,~] = fileparts(img_fname);

[nb_rois, background_indices] = get_inputs();

maxI = 255;
figure;
imshow(img/maxI);

[~,pos,~] = get_n_roi(img,nb_rois);

mat_fname = [img_fname '.mat'];
save(mat_fname,'pos')

imshow_roi(img/maxI,pos)
title(sprintf('ROI positions are saved in \"%s\"',mat_fname),'Interpreter','none')



%%
% ==============================================================
%                        AUXILIARY FUNCTIONS                   %
% ==============================================================
function [nb_rois, background_indices] = get_inputs()
    prompt = {...
        'Number of ROIs:',...
        'Index (or indices) of background region(s):'};
    dlgtitle = 'ROI numbers and indices';
    dims = [1 50];
    definput = {'6','1'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);

    if ~isempty(answer)
        nb_rois = eval(answer{1}); % number of ROIs
        background_indices  = eval(answer{2});
    end
end