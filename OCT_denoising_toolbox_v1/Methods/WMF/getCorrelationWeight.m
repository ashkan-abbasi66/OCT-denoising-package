function weights = getCorrelationWeight(Coeff, windowSize)
% GETCORRELATIONWEIGHT
% Computes the correlation weight for the wavelet multiframe denoising as
% described in:
%
% Markus A. Mayer, Anja Borsdorf, Martin Wagner, Joachim Hornegger,
% Christian Y. Mardin, and Ralf P. Tornow: "Wavelet denoising of multiframe
% optical coherence tomography data", Biomedical Optics Express, 2012
%
% Implementation by Martin Wagner and Markus Mayer,
% Pattern Recognition Lab, University of Erlangen-Nuremberg.
% This version was revised in January 2012.
%
%--------------------------------------------------------------------------

n = 2 * windowSize + 1;
numFrames = size(Coeff(1).approx,3);

if isreal(Coeff(1).approx)
    disp('Computing correlation weight for real coefficients.');
else
    disp('Computing correlation weight for complex coefficients.');
    for l = 1:length(Coeff)
        Coeff(l).approx = abs(Coeff(l).approx);
    end
end

for l = 1:length(Coeff)
    fprintf(1,['Level ' num2str(l) 'Counts' num2str(numFrames)]);
    for m = 1:numFrames
        fprintf(1,'%c','.');
        corrVals = zeros(size(Coeff(l).approx,1), size(Coeff(l).approx,2), numFrames-1);
        windowA = im2col(padarray(Coeff(l).approx(:,:,m), [windowSize windowSize], 'circular'), [n n], 'sliding');
        b = 0;
        for k = 1:numFrames
            if k == m
                b = 1;
                continue;
            end
            windowB = im2col(padarray(Coeff(l).approx(:,:,k), [windowSize windowSize], 'circular'), [n n], 'sliding');
            corrVals(:,:,k-b) = col2im(correlationCol(windowA, windowB), [n n], [size(Coeff(l).approx,1)+n-1 size(Coeff(l).approx,2)+n-1], 'sliding');
        end
        weights(l).approx(:,:,m) = 0.5 .* (median(corrVals, 3)+1);
        clear corrVals;
    end
    disp(' ');
end

end

%--------------------------------------------------------------------------
% Pearsons correlation coefficient of to ROIs of images of the same size.
% Some error handling is performed (division by 0 in the
% correlation computation)

function corr = correlationCol(A, B)
meanA = mean(A);
meanB = mean(B);

varA = var(A);
varB = var(B);

varA(varA < 0.01) = 0.01;
varB(varB < 0.01) = 0.01;

corr = (mean((A - meanA(ones(size(A,1),1), :)) .* (B - meanB(ones(size(B,1),1),:))) ./ sqrt(varA .* varB));

end
