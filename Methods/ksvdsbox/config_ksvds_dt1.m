function [params,run_time] = config_ksvds_dt1(noisy_imgs)

run_time = 0;

tic

% **
% Noise estimation for each frame
% **

Nframes = size(noisy_imgs,3);
sigma_vals = zeros(Nframes,1);

for j = 1:Nframes
    sigma_vals(j) = ceil(function_stdEst(noisy_imgs(:,:,j)));
end

sigma_value = max(sigma_vals);

run_time = run_time + toc;


% **
% Overcomple DCT dictionary
% **

bs = 8; % spatial patch size is [bs, bs]
nf = size(noisy_imgs,3); % number of frames in the input volume

params3d.blocksize = [bs bs nf];                % size of block to operate on
params3d.dictsize = 1000;                       % size of trained dictionary
params3d.basedict{1} = odctdict(bs,10);         % the separable base dictionary -
params3d.basedict{2} = odctdict(bs,10);         % basedict{i} is the base dictionary
params3d.basedict{3} = odctdict(nf,10);         % for the i-th dimension

params3d.initA = speye(params3d.dictsize);  % initial A matrix (identity)


% **
% other parameters
% **
params3d.x = noisy_imgs;        % volume to denoise
params3d.Tdict = 16;            % sparsity of each trained atom - default: 16
params3d.sigma = sigma_value;   % noise power
params3d.maxval = 255;          % maximum sample intensity value
params3d.trainnum = 80000;      % number of training blocks
params3d.iternum = 15;          % number of training iterations
params3d.stepsize = 2;          % distance between adjecent blocks to denoise
params3d.lambda = 0;            % lagrange multiplier
params3d.gain = 1.08;           % noise gain - default: 1.15 ==>  washes out details in the image and reduces its contrast


% Filling PARAMS
params.params3d = params3d;