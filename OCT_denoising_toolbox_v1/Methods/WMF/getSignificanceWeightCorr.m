function weight = getSignificanceWeightCorr( Coeff, Corr, method )
% GETSIGNIFICANCEWEIGHTCORR
% Computes an alternative combination of the correlation and the
% significance weight for multiframe wavelet denoising as described in:
% (This particular combination (combined parameter k) is not mentioned in
% the paper due to length restrictions. It did perform worse compared to
% the compined p method in our evaluation.
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
icount = size(Coeff(1).horizontal, 3);

if isreal(Coeff(1).horizontal)
    disp('Computing significance weight with correlation for real coefficients.');
else
    disp('Computing significance weight with correlation for complex coefficients.');
    for l = 1:length(Coeff)
        Coeff(l).horizontal = abs(Coeff(l).horizontal);
        Coeff(l).vertical = abs(Coeff(l).vertical);
        Coeff(l).diagonal = abs(Coeff(l).diagonal);
    end
end

for l = 1:length(Coeff)
    fprintf(1, [' Level ' num2str(l)]);
    for m = 1:icount
        fprintf(1,'%c','.');
        
        weight(l).horizontal(:,:,m) = methodHandle(Coeff(l).horizontal(:,:,m), Coeff(l).horizontal(:,:,[1:(m-1) m+1:end]), Corr(l).approx(:,:,m));
        weight(l).vertical(:,:,m) = methodHandle(Coeff(l).vertical(:,:,m), Coeff(l).vertical(:,:,[1:(m-1) m+1:end]), Corr(l).approx(:,:,m));
        weight(l).diagonal(:,:,m) = methodHandle(Coeff(l).diagonal(:,:,m), Coeff(l).diagonal(:,:,[1:(m-1) m+1:end]), Corr(l).approx(:,:,m));
    end
    disp(' ');
end

end

%--------------------------------------------------------------------------
% Alternative computation of the significance weight
function weight = getTheta1(img, imgOthers, Corr)
sig = zeros(size(img));
for m = 1:size(imgOthers,3)
    sig = sig + (img - imgOthers(:,:,m)).^2;
end
sig = sqrt(sig ./ (size(imgOthers,3)));
weight = abs(img)./(Corr.*sig);
weight (weight > 1) = 1;

end

%--------------------------------------------------------------------------
% Computation as proposed in the abovementioned paper
function weight = getTheta2(img, imgOthers, Corr)
sig = zeros(size(img));
weight = sig;
for m = 1:size(imgOthers,3)
    sig = sig + (img - imgOthers(:,:,m)).^2;
    tmp = imgOthers(:,:,m);
    tmp(abs(tmp) < 1) = 1;
    weight = weight + abs(1- img./tmp);
end
sig = sqrt(sig ./ (size(imgOthers,3)));
weight = 1-mat2gray(weight ./ (size(imgOthers,3)-1));
weight (abs(img) > Corr .* sig) = 1;

end