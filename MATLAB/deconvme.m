function [O,chi2,J]=deconvme(I,PSF,BG,sigma,alpha,NUMIT)
%DECONVME Maximum entropy deconvolution
% The entropy follows Gull & Skilling's definition.
%
% I: observed data, with gaussian noise
% PSF: point spread function
% BG: background
% sigma: standard deviation of the gaussian noise in observed data
% alpha: relative weight between the goodness-of-fit and the entropy
%
% Return:
% O: estimated image of the object
% chi2: chi square of the fit
% J: chi2/2 - alpha*S, the minimization target
%
% Reference:
% [1]. E. Pantin and J.-L. Starck, Deconvolution of astronomical images using
% the multiscale maximum entropy method, Astronomy and Astrophysics Supp.
% Series, Vol. 118, p575-p585, 1996.
% [2]. M. Raydan and B. F. Svaiter, Relaxed Steepest Descent and Cauchy-Barzilai-Borwein
% Method, Computational Optimization and Applications, 21 (2002), pp. 155-167.
J=zeros(NUMIT,1);
chi2=J;
O=I;
sigma2=sigma^2;
gamma=(sigma2)*0.5;
IPSF=rot90(PSF,2);
if isempty(alpha)
    alpha=0.5*max(PSF(:))/sigma;
end
for k=1:NUMIT
    delta=I-imconv(O,PSF);
    chi2(k)=sum(delta(:).^2)/sigma2;
    lnterm=log(O./BG);
    J(k)=chi2(k)/2 - alpha*sum(sum(O-BG-O.*lnterm));
    G=0-imconv(delta/sigma2,IPSF)+alpha*(1./BG+lnterm-1);
    O=O-gamma*G;
    O=O.*double(O>=BG)+BG*double(O<BG);
    if k>1
        if abs(J(k)-J(k-1))<1e-3
            break
        end
    end
end
return