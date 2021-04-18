function imgRes = waveletMultiFrame(img, varargin)
% WAVELETMULTIFRAME Wavelet based denoising of multiple frame data
% A denoising method based on wavelet decompositions, that works on
% multiple frame data, i.e. all images show the same content, but with
% varying noise. Assumptions made: Noise is uncorrelated between the
% frames.
%
% Parameters:
%   img:     Matrix of noisy frames (each slice of the 3D matrix
%   corresponds to a frame)
%   Additional parameters:
%        'weightMode':      
%          defines the features used for filtering, default: 4
%          0 - Correlation Weight only     
%          1 - Significance Weight only    
%          2 - Combined E                  
%          3 - Combined K                  
%          4 - Combined P                  
%       'maxLevel'  - Maximum decomposition level, default: 3
%       'k'     - Parameter k default: 1 (used in weightmode 1, 2, 3, 4)
%       'p'     - Parameter p default: 5 (used in all weightmodes)
%       'r'     - Parameter r default: 2 (used in weightmode 1, 2, 3, 4)
%       'windowSize' - Window Size (2*windowSize+1) for calculating correlation,
%                 default: 2
%       'theta' - Theta function used for significance weights 
%                 ('theta1' or 'theta2'), default: 'theta2'
%       'basis' - Wavelet basis, default: 'haar'
% 
% A calling example:
% imgResult = waveletMultiFrame(imgFrames, 'k', 1.1, 'p', 1.5, ...
%                         'maxLevel', 5, 'weightMode', 4, 'basis', 'haar');
%
% The description of the algorithm and the parameters can be found in:
% Markus A. Mayer, Anja Borsdorf, Martin Wagner, Joachim Hornegger,
% Christian Y. Mardin, and Ralf P. Tornow: "Wavelet denoising of multiframe
% optical coherence tomography data", Biomedical Optics Express, 2012
%
% Implementation by Martin Wagner and Markus Mayer, 
% Pattern Recognition Lab, University of Erlangen-Nuremberg.
% This version was revised and commented in January 2012.
%
% You may use this code as you want. I would be grateful if you would go to
% my homepage look for articles that you find worth citing in your next
% publication:
% http://www5.informatik.uni-erlangen.de/en/our-team/mayer-markus
% Thanks, Markus

%--------------------------------------------------------------------------
% Read Input Parameters

% Set Standard Values for optional parameters
paramK = 1;
paramP = 5;
paramR = 2;
windowSize = 2;
theta = 'theta2';
basis = 'haar';
maxLevel = 3;
weightMode = 5;
cutBack = 0;

if (~isempty(varargin) && iscell(varargin{1}))
    varargin = varargin{1};
end

% Read Optional Parameters
for k = 1:2:length(varargin)
    if (strcmp(varargin{k}, 'k'))
        paramK = varargin{k+1};
    elseif (strcmp(varargin{k}, 'p'))
        paramP = varargin{k+1};
    elseif (strcmp(varargin{k}, 'r'))
        paramR = varargin{k+1};
    elseif (strcmp(varargin{k}, 'windowSize'))
        windowSize = varargin{k+1};
    elseif (strcmp(varargin{k}, 'theta'))
        theta = varargin{k+1};
    elseif (strcmp(varargin{k}, 'basis'))
        basis = varargin{k+1};
    elseif (strcmp(varargin{k}, 'maxLevel'))
        maxLevel = varargin{k+1};
    elseif (strcmp(varargin{k}, 'weightMode'))
        weightMode =  varargin{k+1};
    end
end

if strcmp(basis, 'dualTree')
    [Faf, Fsf] = FSfarras;
    [af, sf] = dualfilt1;
end

%--------------------------------------------------------------------------
% Perform Wavelet Decomposition

numFrames = size(img, 3);

% Most wavelet transforms require a special size 
% e.g. mod2 == 0 for each level - this is just example code how one could
% work around such limitations. 
if (maxLevel == 5) && (size(img,1) == 496) 
    imgTemp = zeros(512, size(img,2), numFrames);
    imgTemp(1:496,:,:) = img;
    img = imgTemp;
    clear imgTemp;
    cutBack = 1;
end

for i = numFrames:-1:1
    if strcmp(basis, 'dualTree')
        % Dual tree transformation taken from:
        % WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
        % http://taco.poly.edu/WaveletSoftware/
        
        % Normalization
        x = double(img(:,:,i))/sqrt(2);
        
        % Tree 1
        [x1 wt{1}{1}] = afb2D(x, Faf{1});      % Stage 1
        Coeff(1).approx(:,:,i) = x1;
        for l = 2:maxLevel
            [x1 wt{l}{1}] = afb2D(x1, af{1});  % Remaining stages
            Coeff(l).approx(:,:,i) = x1;
        end
        
        
        % Tree 2
        [x2 wt{1}{2}] = afb2D(x, Faf{2});      % Stage 1
        Coeff(1).approx(:,:,i) = complex(Coeff(1).approx(:,:,i), x2);
        for l = 2:maxLevel
            [x2 wt{l}{2}] = afb2D(x2, af{2});  % Remaining stages
            Coeff(l).approx(:,:,i) = complex(Coeff(l).approx(:,:,i), x2);
        end
        
        % Sum and difference
        for l = 1:maxLevel
            for direction = 1:3
                A = wt{l}{1}{direction};
                B = wt{l}{2}{direction};
                wt{l}{1}{direction} = (A+B)/sqrt(2);
                wt{l}{2}{direction} = (A-B)/sqrt(2);
            end
        end
        
        for l = maxLevel:-1:1
            Coeff(l).horizontal(:,:,i) = complex(wt{l}{1}{1},wt{l}{2}{1});
            Coeff(l).vertical(:,:,i) = complex(wt{l}{1}{2},wt{l}{2}{2});
            Coeff(l).diagonal(:,:,i) = complex(wt{l}{1}{3},wt{l}{2}{3});
        end
    else % All other wavelet transformations
        [wt(i).approx, wt(i).horizontal, wt(i).vertical, wt(i).diagonal] = swt2(img(:,:,i),maxLevel, basis);
        for l = maxLevel:-1:1
            Coeff(l).approx(:,:,i)  = wt(i).approx(:,:,l);
            Coeff(l).horizontal(:,:,i) = wt(i).horizontal(:,:,l);
            Coeff(l).vertical(:,:,i) = wt(i).vertical(:,:,l);
            Coeff(l).diagonal(:,:,i) = wt(i).diagonal(:,:,l);
        end
    end
end

% Free Memory that is not needed anymore
clear wt;
clear img;

%--------------------------------------------------------------------------
% Compute weights

if (weightMode ~= 1)
    weightsCorr = getCorrelationWeight(Coeff, windowSize);
end

if (weightMode ~= 0)
    if (weightMode ~= 3) 
        weightsSignificance = getSignificanceWeight(Coeff, paramK, theta);
    else
        for l = 1:maxLevel
            weightsCorr(l).approx = (paramK .* (1 - mat2gray(weightsCorr(l).approx).^paramP) + 1);
        end
        weightsSignificance = getSignificanceWeightCorr(Coeff, weightsCorr , theta);
        clear weightsCorr;
    end
end


%--------------------------------------------------------------------------
% Combine weights if necessary

switch weightMode
    case 0
        for l = 1:maxLevel
            weightsCorr(l).approx = weightsCorr(l).approx .^ paramP;
        end
    case 1
        for l = 1:maxLevel
            weightsSignificance(l).horizontal = weightsSignificance(l).horizontal .^ paramR;
            weightsSignificance(l).vertical = weightsSignificance(l).vertical .^ paramR;
            weightsSignificance(l).diagonal = weightsSignificance(l).diagonal .^ paramR;
        end
    case 2
        for l = 1:maxLevel
            weightsSignificance(l).horizontal = (weightsCorr(l).approx .^ paramP) .* exp(-(1-weightsSignificance(l).horizontal.^2).^paramR);
            weightsSignificance(l).vertical = (weightsCorr(l).approx .^ paramP) .* exp(-(1-weightsSignificance(l).vertical.^2).^paramR);
            weightsSignificance(l).diagonal = (weightsCorr(l).approx .^ paramP) .* exp(-(1-weightsSignificance(l).diagonal.^2).^paramR);
        end
        clear weightsCorr;
    case 3
        for l = 1:maxLevel
            weightsSignificance(l).horizontal = weightsSignificance(l).horizontal .^ paramR;
            weightsSignificance(l).vertical = weightsSignificance(l).vertical .^ paramR;
            weightsSignificance(l).diagonal = weightsSignificance(l).diagonal .^ paramR;
        end
    case 4
        for l = 1:maxLevel
            weightsSignificance(l).horizontal = weightsCorr(l).approx  .^ (1 + paramP .* (1 - weightsSignificance(l).horizontal.^paramR));
            weightsSignificance(l).vertical = weightsCorr(l).approx  .^ (1 + paramP .* (1 - weightsSignificance(l).vertical.^paramR));
            weightsSignificance(l).diagonal = weightsCorr(l).approx  .^ (1 + paramP .* (1 - weightsSignificance(l).diagonal.^paramR));
        end
        clear weightsCorr;
end

%--------------------------------------------------------------------------
% Wavelet Weighting

if (weightMode == 0) 
    for l = maxLevel:-1:1
        Coeff(l).horizontal = Coeff(l).horizontal .* weightsCorr(l).approx;
        Coeff(l).vertical = Coeff(l).vertical .* weightsCorr(l).approx;
        Coeff(l).diagonal = Coeff(l).diagonal .* weightsCorr(l).approx;
    end
else
    for l = maxLevel:-1:1
        Coeff(l).horizontal = Coeff(l).horizontal .* weightsSignificance(l).horizontal;
        Coeff(l).vertical = Coeff(l).vertical .* weightsSignificance(l).vertical;
        Coeff(l).diagonal = Coeff(l).diagonal .* weightsSignificance(l).diagonal;
    end
end

%--------------------------------------------------------------------------
% Averaging

% Calculate mean approximation coefficients and delete unused
% coefficients.
for l = maxLevel:-1:1
    meanApprox{l} = mean(Coeff(l).approx, 3);
    Coeff(l).approx = [];
end

for l = maxLevel:-1:1
    meanCoeffHorizontal{l} = mean(Coeff(l).horizontal, 3);
    meanCoeffVertical{l} = mean(Coeff(l).vertical, 3);
    meanCoeffDiagonal{l} = mean(Coeff(l).diagonal, 3);
end


if strcmp(basis, 'dualTree')
    for l = maxLevel:-1:1
        wt{l}{1}{1} = real(meanCoeffHorizontal{l});
        wt{l}{2}{1} = imag(meanCoeffHorizontal{l});
        
        wt{l}{1}{2} = real(meanCoeffVertical{l});
        wt{l}{2}{2} = imag(meanCoeffVertical{l});
        
        wt{l}{1}{3} = real(meanCoeffDiagonal{l});
        wt{l}{2}{3} = imag(meanCoeffDiagonal{l});
    end
    
    wt{maxLevel+1}{1} = real(meanApprox{maxLevel});
    wt{maxLevel+1}{2} = imag(meanApprox{maxLevel});
    
    imgRes = idualtree2D(wt, maxLevel, Fsf, sf);
else
    for l = maxLevel:-1:1
        meanApproxTemp(:,:,l) = meanApprox{l};
        meanCoeffHorizontalTemp(:,:,l) = meanCoeffHorizontal{l};
        meanCoeffVerticalTemp(:,:,l) = meanCoeffVertical{l};
        meanCoeffDiagonalTemp(:,:,l) = meanCoeffDiagonal{l};
    end
    
    imgRes = iswt2(meanApproxTemp, meanCoeffHorizontalTemp, meanCoeffVerticalTemp, meanCoeffDiagonalTemp, basis);
end

% Cut image size workaround back
if cutBack
    imgRes = imgRes(1:496, :);
end

end
