clc
clear all

addpath('./utils')
addpath('./dependency/fkmeans')
addpath('./dependency/kmeans')
addpath('./dependency/poblano_toolbox')
addpath('./tensor_toolbox')

[imaVOL,scaninfo] = loadminc('t1_icbm.mnc');
input=double(imaVOL(:,:,1:60));
input_clean= normalized(input);
sigma=0.2;
msi_sz  =  size(input_clean);
noisy_input = input_clean + sigma * randn(msi_sz);  % add Gaussian noise


%%

sigma=0.3;
tic
[ clean_img, ~, ~, ~, ~ ] = TensorDL(noisy_input, sigma);
toc

figure;imshow(noisy_input(:,:,1),[])
figure;imshow(clean_img(:,:,1),[])