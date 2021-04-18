function params_ksvds = get_params_ksvds(noisy_imgs,sigma_value)

% **
% Overcomple DCT dictionary
% **

bs = 8; % spatial patch size is [bs, bs]
nf = size(noisy_imgs,3); % number of frames in the input volume

params_ksvds.blocksize = [bs bs nf];                % size of block to operate on
params_ksvds.dictsize = 1000;                       % size of trained dictionary
params_ksvds.basedict{1} = odctdict(bs,10);         % the separable base dictionary -
params_ksvds.basedict{2} = odctdict(bs,10);         % basedict{i} is the base dictionary
params_ksvds.basedict{3} = odctdict(nf,10);         % for the i-th dimension

params_ksvds.initA = speye(params_ksvds.dictsize);  % initial A matrix (identity)


% **
% other parameters
% **
params_ksvds.x = noisy_imgs;        % volume to denoise
params_ksvds.Tdict = 16;            % sparsity of each trained atom - default: 16
params_ksvds.sigma = sigma_value;   % noise power
params_ksvds.maxval = 255;          % maximum sample intensity value
params_ksvds.trainnum = 80000;      % number of training blocks
params_ksvds.iternum = 5%15;          % number of training iterations
params_ksvds.stepsize = 2;          % distance between adjecent blocks to denoise
params_ksvds.lambda = 0;            % lagrange multiplier
params_ksvds.gain = 1.08;           % noise gain - default: 1.15 ==>  washes out details in the image and reduces its contrast
