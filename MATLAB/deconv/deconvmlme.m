function [O,J]=deconvmlme(I,PSF,varargin)
%DECONVMLME Maximum likelihood deconvolution with maximum entropy constraint.
%The iteration stops once number of iterations exceeds NUMIT or correction
%becomes smaller than EPS.
%
%Syntax:
%  deconvml(I, PSF, NUMIT)
%  deconvml(I, PSF, LL, NUMIT)
%  deconvml(I, PSF, LL, UL, NUMIT)
%  deconvml(I, PSF, InitO, LL, UL, NUMIT)
%  deconvml(I, PSF, EPS)
%  deconvml(I, PSF, LL, EPS)
%  deconvml(I, PSF, LL, UL, EPS)
%  deconvml(I, PSF, InitO, LL, UL, EPS)
%  deconvml(I, PSF, InitO, LL, UL, alpha, EPS, NUMIT)
%
%Input arguments:
%  I - observed data with Poisson noise
%  PSF - point spread function
%  LL - lower limit
%  UL - upper limit
%  InitO - initial estimate of image
%  alpha - maximum entropy constraint factor
%  NUMIT - maximum number of iterations ( >1 )
%  EPS - maximum precision requirement  ( <1 )
%
%Return:
%  O

%% parsing input arguments
NUMIT=100;
EPS=1e-3;
LL=[];
UL=[];
InitO=[];
alpha=[];
switch nargin
    case 2
    case 3
        if varargin{nargin-2}>1
            NUMIT=varargin{nargin-2};
        else
            EPS=varargin{nargin-2};
        end
    case 4
        LL=varargin{1};
        if varargin{nargin-2}>1
            NUMIT=varargin{nargin-2};
        else
            EPS=varargin{nargin-2};
        end
    case 5
        LL=varargin{1};
        UL=varargin{2};
        if varargin{nargin-2}>1
            NUMIT=varargin{nargin-2};
        else
            EPS=varargin{nargin-2};
        end
    case 6
        InitO=varargin{1};
        LL=varargin{2};
        UL=varargin{3};
        if varargin{nargin-2}>1
            NUMIT=varargin{nargin-2};
        else
            EPS=varargin{nargin-2};
        end
    case 8
        InitO=varargin{1};
        LL=varargin{2};
        UL=varargin{3};
        alpha=varargin{4};
        EPS=varargin{5};
        NUMIT=varargin{6};
    otherwise
        error('Syntax error.')
end
if isempty(alpha)
    alpha=1e-3;
end
if isempty(LL)
    disp('Please drag a rectangle on the image for background estimation.')
    [~,LL]=imflux(I,[]);
end
if isempty(InitO)
    InitO=I;
end
if isscalar(LL)
    LL=ones(size(I))*LL;
end
if isscalar(UL)
    UL=ones(size(I))*UL;
end
J=zeros(NUMIT,1);
O=InitO;
OTF=psf2otf(PSF);
rng('shuffle')
u=rand(NUMIT,1);
for k=1:NUMIT
    R=real(ifft2(OTF.*fft2(O)));
    R(R==0)=eps;
    phi=abs(real(ifft2(fft2(I./R + eps).*conj(OTF)))...
        - alpha*log(O)).^(0.5*u(k));
    J(k)=max(abs(phi(:)-1));
    O=O.*phi;
    O=max(O,LL);
    if ~isempty(UL),    O=min(O,UL);end
    if J(k)<=EPS
        disp('Precision requirement has been satisfied.')
        break
    end
    if k>2
        if J(k)>=J(k-1)
            disp('Convergence has been reached.')
%             break
        end
    end
end
O=O*sum(I(:))/sum(O(:));
J=J(1:k);
return