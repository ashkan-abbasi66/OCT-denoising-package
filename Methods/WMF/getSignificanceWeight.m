function weight = getSignificanceWeight(Coeff, paramK, method)
% GETSIGNIFICANCEWEIGHT
% Computes the significance weight for the wavelet multiframe denoising as
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

if (strcmp(method, 'theta1'))
    methodHandle = @getTheta1;
elseif (strcmp(method, 'theta2'))
    methodHandle = @getTheta2;
else
    disp('Error: Unkown function theta!!');
end
numFrames = size(Coeff(1).horizontal, 3);

if isreal(Coeff(1).horizontal)
    disp('Computing significance weight for real coefficients.');
else
    disp('Computing significance weight for complex coefficients.');
    for l = 1:length(Coeff)
        Coeff(l).horizontal = abs(Coeff(l).horizontal);
        Coeff(l).vertical = abs(Coeff(l).vertical);
        Coeff(l).diagonal = abs(Coeff(l).diagonal);
    end
end

for l = 1:length(Coeff)
    fprintf(1, [' Level ' num2str(l)]);
    for m = 1:numFrames
        fprintf(1,'%c','.');
        weight(l).horizontal(:,:,m) = methodHandle(Coeff(l).horizontal(:,:,m), Coeff(l).horizontal(:,:,[1:(m-1) m+1:end]), paramK);
        weight(l).vertical(:,:,m) = methodHandle(Coeff(l).vertical(:,:,m), Coeff(l).vertical(:,:,[1:(m-1) m+1:end]), paramK);
        weight(l).diagonal(:,:,m) = methodHandle(Coeff(l).diagonal(:,:,m), Coeff(l).diagonal(:,:,[1:(m-1) m+1:end]), paramK);
    end
    disp(' ');
end
end

%--------------------------------------------------------------------------
% Alternative computation of the significance weight
function weight = getTheta1(img, imgOthers, paramK)
sig = zeros(size(img));
for m = 1:size(imgOthers,3)
    sig = sig + (img - imgOthers(:,:,m)).^2;
end
sig = sqrt(sig ./ (size(imgOthers,3)));
weight = abs(img)./(paramK.*sig);
weight (weight > 1) = 1;

end

%--------------------------------------------------------------------------
% Computation as proposed in the abovementioned paper
function weight = getTheta2(img, imgOthers, paramK)
sig = zeros(size(img));
weight = sig;
for m = 1:size(imgOthers,3)
    sig = sig + (img - imgOthers(:,:,m)).^2;
    tmp = imgOthers(:,:,m);
    tmp(abs(tmp) < 1) = 1;
    weight = weight + abs(1- img./tmp);
end
sig = sqrt(sig ./ (size(imgOthers,3)));
weight = 1-mat2gray(weight ./ (size(imgOthers,3)));
weight (abs(img) > paramK .* sig) = 1;

end