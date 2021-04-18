function noisevar = evar(y)

%EVAR   Noise variance estimation.
%   Assuming that the deterministic function Y has additive Gaussian noise,
%   EVAR(Y) returns an estimated variance of this noise.
%
%   Note:
%   ----
%   A thin-plate smoothing spline model is used to smooth Y. It is assumed
%   that the model whose generalized cross-validation score is minimum can
%   provide the variance of the additive noise. A few tests showed that
%   EVAR works very well with "not too irregular" functions.
%
%   Examples:
%   --------
%   % 1D signal
%   n = 1e6; x = linspace(0,100,n);
%   y = cos(x/10)+(x/50);
%   var0 = 0.02; % noise variance
%   yn = y + sqrt(var0)*randn(size(y));
%   evar(yn) % estimated variance
%
%   % 2D function
%   [xp,yp] = deal(0:.02:1);
%   [x,y] = meshgrid(xp,yp);
%   f = exp(x+y) + sin((x-2*y)*3);
%   var0 = 0.04; % noise variance
%   fn = f + sqrt(var0)*randn(size(f));
%   evar(fn) % estimated variance
%
%   % 3D function
%   [x,y,z] = meshgrid(-2:.2:2,-2:.2:2,-2:.2:2);
%   f = x.*exp(-x.^2-y.^2-z.^2);
%   var0 = 0.5; % noise variance
%   fn = f + sqrt(var0)*randn(size(f));
%   evar(fn) % estimated variance
%
%   % Other examples
%   Click <a href="matlab:web('http://www.biomecardio.com/matlab/evar.html')">here</a> for more examples
%
%   Note:
%   ----
%   EVAR is only adapted to evenly-gridded 1-D to N-D data.
%
%   See also VAR, STD, SMOOTHN
%
%   -- Damien Garcia -- 2008/04, revised 2009/10

error(nargchk(1,1,nargin));

d = ndims(y);
siz = size(y);
S = zeros(siz);
for i = 1:d
    siz0 = ones(1,d);
    siz0(i) = siz(i);
    S = bsxfun(@plus,S,cos(pi*(reshape(1:siz(i),siz0)-1)/siz(i)));
end
S = 2*(d-S);

% N-D Discrete Cosine Transform of Y
if exist('dctn','file')
    DCTy = dctn(y);
else
    error('MATLAB:evar:MissingFunction',...
        ['DCTN is required. Download <a href="matlab:web(''',...
        'http://www.biomecardio.com/matlab/dctn.html'')">DCTN</a>.'])
end

fminbnd(@func,-38,38);

function score = func(L)
    % Generalized cross validation score
    M = 1-1./(1+10^L*S.^2);
    noisevar = mean(DCTy(:).^2.*M(:).^2);
    score = noisevar/mean(M(:))^2;
end

end

