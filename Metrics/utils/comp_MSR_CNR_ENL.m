function [MSR,CNR,ENL]=comp_MSR_CNR_ENL(roi,background_indices)
% 
% Computes three non-referenced quality assessment metrics (MSR, CNR, and
% ENL).
% 
% 
% INPUT
%   roi: cell array of ROI vectors
%   background_indices: vector indicating which ROIs are background ROI.
%   e.g., background_indices = [1] means that the first region is the only
%   background region.
%   
% OUTPUT
%   MSR: mean-to-standard-deviation ratio [1]
%   CNR: Contrast to Noise Ratio [2]
%   ENL: Equivalent number of looks [2]
% 
% 
% References:
%   [1] Bao, Paul, and Lei Zhang. "Noise reduction for magnetic resonance
%   images via adaptive multiscale products thresholding." IEEE
%   transactions on medical imaging 22.9 (2003): 1089-1099.
% 
%   [2] Pizurica, Aleksandra, et al. "Multiresolution denoising for optical
%    coherence tomography: a review and evaluation." Current Medical
%    Imaging 4.4 (2008): 270-284.
% 
% 


N=numel(roi); % number of ROI vectors

foreground_indices=setdiff(1:N,background_indices);


% ===========
% Compute MSR
% MSR = average of mean to standard deviation ratios on foreground regions

f_mean=0; % mean of intensity values in the foreground ROIs
f_std=0;  % std of intensity values in the foreground ROIs
MSR=0;

for j=foreground_indices
    m=mean(roi{j});
    s=std(roi{j});
    MSR=MSR+m/s;
    f_mean=f_mean+m;
    f_std=f_std+s;
end


Nf=numel(foreground_indices);% number of foreground regions

MSR=MSR/Nf;
f_mean=f_mean/Nf;
f_std=f_std/Nf;


% ===========
% Compute CNR
% CNR = measures the contrast between foreground regions and background
% regions.

% background ROI vectors
broi=roi{background_indices};

CNR=abs(f_mean-mean(broi))/sqrt(0.5*(f_std^2+std(broi)^2));


% ===========
% Compute ENL
ENL = mean(broi)^2/std(broi)^2;