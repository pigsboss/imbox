function [O,chi2]=deconvchi2(I,PSF,varargin)
%DECONVCHI2 Using chi square fitting to solve deconvolution problem[1].
%
% The minimization of chi2 is implemented with steepest descent method with
% random step size[2].
% 
%Syntax:
% deconvchi2(I, PSF)
% deconvchi2(I, PSF, NUMIT)
% deconvchi2(I, PSF, BG, NUMIT)
% deconvchi2(I, PSF, BG, sigma, NUMIT)
% deconvchi2(I, PSF, LL, UL, sigma, NUMIT)
% deconvchi2(I, PSF, InitO, LL, UL, sigma, NUMIT)
%
% I: observed data, with gaussian noise
% PSF: point spread function
% BG (LL): background (lower limit)
% sigma: standard deviation of the gaussian noise in observed data
% UL: upper limit
% InitO: initial estimate
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

NUMIT=20;
sigma=[];
bg=0;
ul=[];
InitO=[];
switch nargin
    case 2
    case 3
        NUMIT=varargin{nargin-2};
    case 4
        bg=varargin{1};
        NUMIT=varargin{nargin-2};
    case 5
        bg=varargin{1};
        sigma=varargin{nargin-3};
        NUMIT=varargin{nargin-2};
    case 6
        bg=varargin{1};
        ul=varargin{2};
        sigma=varargin{nargin-3};
        NUMIT=varargin{nargin-2};
    case 7
        InitO=varargin{1};
        bg=varargin{2};
        ul=varargin{3};
        sigma=varargin{nargin-3};
        NUMIT=varargin{nargin-2};
    otherwise
        error('Syntax error.')
end
if isempty(InitO)
    InitO=I;
end
if isempty(sigma)
    sigma=sigmaclipping(I-sdwt2(I,1),3);
    disp(sigma)
end
rng('shuffle')
chi2=zeros(NUMIT,1);
O=InitO;
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
    if ~isempty(ul)
        O=O.*double(O<=ul)+ul*double(O>ul);
    end
    if k>1
        if (abs(chi2(k)-chi2(k-1))/chi2(k-1))<1e-8
            disp('chi-square converges.');
            break
        end
    end
end
chi2=chi2(1:k);
return