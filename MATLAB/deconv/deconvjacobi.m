function O=deconvjacobi(I,PSF,varargin)
%DECONVJACOBI Deconvolution using Jacobi method.
%Syntax:
%  deconvjacobi(I, PSF)
%  deconvjacobi(I, PSF, LL)
%  deconvjacobi(I, PSF, LL, UL)
%  deconvjacobi(I, PSF, InitO, LL, UL)
%  deconvjacobi(I, PSF, InitO, LL, UL, L)
%  deconvjacobi(I, PSF, InitO, LL, UL, L, omega)
%  deconvjacobi(I, PSF, InitO, LL, UL, L, omega, NUMIT)
%
%Input arguments:
%  I - observed data
%  PSF - point spread function
%  LL - lower limit of true image
%  UL - upper limit of true image
%  InitO - initial estimate of true image
%  L - correlation transform order
%  omega - relaxation factor
%    (0,1) - under relaxation
%    =1 (default) - relaxation (original jacobi method)
%    >1    - over relaxation
%  NUMIT - max number of iterations
%
%Returns:
%  O - estimate of true image
%
%Reference:
%  [1] T.-P. Li and M. Wu, A Direct Restoration Method for Spectral and Image
%    Analysis, Astrophysics and Space Science, Vol. 206, 1993.

%% default parameters
LL=0;
UL=[];
L=0;
InitO=[];
omega=1;
NUMIT=100;
path('./lib',path)

%% parsing input arguments
switch nargin
    case 2
    case 3
        LL=varargin{1};
    case 4
        LL=varargin{1};
        UL=varargin{2};
    case 5
        InitO=varargin{1};
        LL=varargin{2};
        UL=varargin{3};
    case 6
        InitO=varargin{1};
        LL=varargin{2};
        UL=varargin{3};
        L=varargin{4};
    case 7
        InitO=varargin{1};
        LL=varargin{2};
        UL=varargin{3};
        L=varargin{4};
        omega=varargin{5};
    case 8
        InitO=varargin{1};
        LL=varargin{2};
        UL=varargin{3};
        L=varargin{4};
        omega=varargin{5};
        NUMIT=varargin{6};
    otherwise
        error('Syntax error.')
end
if isempty(LL)
    LL=0;
end
if isempty(omega)
    omega=1;
end
if isempty(L)
    L=0;
end
if isempty(InitO)
    InitO=I;
end
if isempty(NUMIT)
    NUMIT=100;
end
if isscalar(LL)
    LL=LL*ones(size(I));
end
if ~isempty(UL)
    if isscalar(UL)
        UL=UL*ones(size(I));
    end
end
%% L-folder correlation transform
for l=1:L
    I=imconv(I,rot90(PSF,2));
    PSF=imconv(PSF,rot90(PSF,2));
end
PSFRmd=ifftshift(PSF);
PSFDiag=PSFRmd(1,1);
PSFRmd(1,1)=0;
PSFRmd=fftshift(PSFRmd);
NTF=sum(I(:)-LL(:)); % net total flux
rng('shuffle')
O=InitO;
for k=1:NUMIT
    u=rand(1);
    sigma=omega*u*(I-imconv(O,PSFRmd))/PSFDiag;
    O=(1-omega*u)*O+sigma;
    O=O.*double(O>LL)+LL.*double(O<LL);
    if ~isempty(UL)
%         O=O.*double(O<UL)+UL.*double(O>UL);
        O=((O-LL)/max(O(:)-LL(:))).*UL+LL;
    end
    O=(O-LL)/sum(O(:)-LL(:))*NTF+LL;
    if max(abs(sigma(:)))/max(abs(O(:)))<=1e-3
        disp(['Convergence has been reached in ',int2str(k),'-th step.'])
%         break
    end
end
return