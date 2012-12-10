function [O,chi2]=deconvchi2(I,PSF,bg,sigma,NUMIT)
%DECONVCHI2 Using chi square fitting to solve deconvolution problem[1].
%
% The minimization of chi2 is implemented with steepest descent method with
% random step size[2].
% 
% I: observed data, with gaussian noise
% PSF: point spread function
% BG: background
% sigma: standard deviation of the gaussian noise in observed data
%
% Return:
% O: estimated image of the object
% chi2: chi square of the fit
%
% Reference:
% [1]. E. Pantin and J.-L. Starck, Deconvolution of astronomical images using
% the multiscale maximum entropy method, Astronomy and Astrophysics Supp.
% Series, Vol. 118, p575-p585, 1996.
% [2]. M. Raydan and B. F. Svaiter, Relaxed Steepest Descent and Cauchy-Barzilai-Borwein
% Method, Computational Optimization and Applications, 21 (2002), pp. 155-167.
rng('shuffle')
chi2=zeros(NUMIT,1);
O=I;
sigma2=sigma^2;
IPSF=rot90(PSF,2);
for k=1:NUMIT
    I_virt=imconv(O,PSF);
    delta=I-I_virt;
    chi2(k)=sum(delta(:).^2)/sigma2;
    G=-imconv(delta/sigma2,IPSF);
    Q=imconv(G,PSF);
    gamma=(Q(:)'*I_virt(:) - I(:)'*Q(:))/sum(Q(:).^2);
    O=O-rand(1)*2*gamma*G;
    O=O.*double(O>=bg)+bg*double(O<bg);
    if k>1
        if abs(chi2(k)-chi2(k-1))<1e-3
            disp('chi-square converges.');
            break
        end
    end
end
return