function params = get_params_wmf_dt1(~,~)
% Returns required parameters to run the WMF denoising method on dt1
%


params.sizeDivisibleBy = 32;

params.k = 1.1;% Default: 1.1
params.p = 1.5;% Default: 1.5
params.maxLevel = 5;
params.weightMode = 2; % 0/2/4
params.basis = 'dualTree'; % haar, dualTree

