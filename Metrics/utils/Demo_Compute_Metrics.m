% 
% This script shows how to compute the following image quality assessment
% metrics:
% 
% PSNR, SSIM, MSR, CNR, ENL, TP, and EP.
% 
% 
% NOTE: (How to select regions of interest (ROIs))
% ROIs are generally divided into two categories: 
%   1) foreground region (desired region)
%   2) background region (undesired region). 
% Some metrics need one type of ROI while the other metrics need both types 
% of ROIs. Please, appropriately select ROIs and set the background ROIs
% indices ("background_indices") correctly in this script.
% 
% 

%% 

clear
clc


dataset_path = '../../Datasets/dt1_Bioptigen_SDOCT/3/';
result_path = '../../Results/benchmark_bm4d_iidnoise';

% noisy image
imn_path = fullfile(dataset_path,'test.tif');
imn = double(imread(imn_path));

% ground-truth image
im_path = fullfile(dataset_path,'average.tif');
im = double(imread(im_path));

% output/denoised image
im_out_path = fullfile(result_path,'03.tif');
im_out = double(imread(im_out_path));


maxI = 255;

% % figure,imshow(imn/maxI),title('noisy image')
% % figure,imshow(im/maxI),title('ground-truth image')
% % figure,imshow(im_out/maxI),title('denoised image')


%% Compute PSNR and SSIM metrics

[PSNR,SSIM] = comp_psnr_ssim(im,im_out,maxI)


%% Select some ROIs

input_img = imn;

N = 6; % number of ROIs that you want to select

figure;
imshow(input_img/maxI);

[~,pos,~] = get_n_roi(input_img,N);

% Uncomment to show the selected ROI over the output image
% show_roi(im_out/maxI,pos)

% Uncomment to store position array in a .MAT file
% save('demo_saved_position.mat','pos')



%% Compute MSR, CNR, and ENL metrics for the output image

% load('demo_saved_position.mat')

[roi,~] = get_roi_pos(im_out,pos);

background_indices = [1]; % indices of ROIs which are background regions

[MSR,CNR,ENL] = comp_MSR_CNR_ENL(roi,background_indices)



%% Texture Preservation and Edge Preservation metrics

ref_image = imn;

TP = comp_TP(im_out,ref_image,pos,background_indices)

EP = comp_EP(im_out,ref_image,pos,background_indices)


