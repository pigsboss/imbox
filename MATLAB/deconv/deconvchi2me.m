function [O,chi2,J]=deconvchi2me(I,PSF,BG,sigma,alpha,NUMIT)
%DECONVCHI2ME Maximum entropy deconvolution
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
IPSF=rot90(PSF,2);
if isempty(alpha)
    alpha=0.5*max(PSF(:))/sigma;
end
rng('shuffle')
for k=1:NUMIT
    I_virt=imconv(O,PSF);
    delta=I-I_virt;
    chi2(k)=sum(delta(:).^2)/sigma2;
    lnterm=log(O./BG);
    J(k)=chi2(k)/2 - alpha*sum(sum(O-BG-O.*lnterm));
    G=-imconv(delta/sigma2,IPSF) + alpha*lnterm;
    Q=imconv(G,PSF);
    %gamma=(Q(:)'*I_virt(:) - I(:)'*Q(:))/sum(Q(:).^2);
    gamma=(BG*(Q(:)'*I_virt(:)-Q(:)'*I(:))+alpha*sigma2*(Q(:)'*O(:)-sum(G(:))))/...
        (BG*sum(Q(:).^2)+alpha*sigma2*sum(G(:).^2));
    O=O-rand(1)*2*gamma*G;
    O=O.*double(O>=BG)+BG*double(O<BG);
    if k>1
        if (abs(J(k)-J(k-1))/J(k-1))<1e-4
            disp('chi-square-max-entropy converges.');
            break
        end
    end
end
chi2=chi2(1:k);
J=J(1:k);
return