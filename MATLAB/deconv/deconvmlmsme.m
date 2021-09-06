function [O,J]=deconvmlmsme(I,PSF,varargin)
%DECONVMLMSME Maximum likelihood deconvolution with multiscale maximum
%entropy constraint.
%
%Syntax:
%  deconvmlmsme(I, PSF)
%  deconvmlmsme(I, PSF, NUMIT)
%  deconvmlmsme(I, PSF, InitO, LL, MaxScale, NUMIT)
%  deconvmlmsme(I, PSF, InitO, LL, UL, sigma, MaxScale, alpha, NUMIT)
%
%Input arguments:
%  I - observed data, with Poisson noise
%  PSF - point spread function
%  NUMIT - maximum number of iterations
%  InitO - initial estimate of image
%  LL - lower limit
%  UL - upper limit
%  sigma - standard deviation of background noise
%  MaxScale - maximum wavelet scale
%  alpha - parameter balancing the entropy and the log-likelihood
%
%Returns:
%
%References:
%

%% set path
path('./lib',path)
path('./lib/msme',path)
path('./lib/wavelet',path)
%% set default parameters
NUMIT=10;
LL=[];
UL=[];
MaxScale=[];
alpha=[];
sigma=[];
InitO=[];
% sigma=imconv(sigma,fspecial('gauss',size(I),5));
%% parsing input arguments
switch nargin
    case 2
    case 3
        NUMIT=varargin{1};
    case 6
        InitO=varargin{1};
        LL=varargin{2};
        MaxScale=varargin{3};
        NUMIT=varargin{4};
    case 9
        InitO=varargin{1};
        LL=varargin{2};
        UL=varargin{3};
        sigma=varargin{4};
        MaxScale=varargin{5};
        alpha=varargin{6};
        NUMIT=varargin{7};
    otherwise
        error('Syntax error.')
end
if isempty(LL)
    [~,LL]=imflux(I,[]);
end
if isempty(sigma)
    sigma=sqrt(I);
end
if isscalar(sigma)
    sigma=sigma*ones(size(I));
end
if isempty(MaxScale)
    MaxScale=lastpow2(min(size(I)));
end
if isempty(alpha)
    alpha=0.5*max(PSF(:))./sigma;
end
if isscalar(alpha)
    alpha=alpha*max(PSF(:))./sigma;
end
if isempty(NUMIT)
    NUMIT=10;
end
if isempty(InitO)
    InitO=I;
end
if isscalar(LL)
    LL=LL*ones(size(I));
end
if isscalar(UL)
    UL=UL*ones(size(I));
end
%% prepare workspace
[d1,d2]=size(I);
J=zeros(NUMIT,1);
kPadSz=2^(MaxScale+1);
BSF=bsplinefilter(size(PSF));
BSF=BSF(:,:,1:MaxScale);
DDF=zeros(size(I));DDF(1)=1;DDF=fftshift(DDF);
WTF=DDF-BSF(:,:,1);
IWTF=zeros(size(BSF));
IWTF(:,:,1)=rot90(WTF,2);
for k=2:MaxScale
    WTF=BSF(:,:,k-1)-BSF(:,:,k);
    IWTF(:,:,k)=rot90(WTF,2);
end
[~,w]=dwt2(I,MaxScale);
sigmams=mssigma(MaxScale);
msigma=zeros(size(w));
for k=1:MaxScale
    msigma(:,:,k)=sigma*sigmams(k);
end
M=mssupht(w,msigma,BSF);
A=1-M;
I=padarray(I,[kPadSz kPadSz],'symmetric','both');
sigma=padarray(sigma,[kPadSz kPadSz],'symmetric','both');
msigma=padarray(msigma,[kPadSz kPadSz],'symmetric','both');
alpha=padarray(alpha,[kPadSz kPadSz],'symmetric','both');
PSF=padarray(PSF,[kPadSz kPadSz],0,'both');
A=padarray(A,[kPadSz kPadSz],1,'both');
IWTF=padarray(IWTF,[kPadSz kPadSz],0,'both');
O=InitO;
O=padarray(O,[kPadSz kPadSz],'symmetric','both');
LL=padarray(LL,[kPadSz kPadSz],'symmetric','both');
if ~isempty(UL)
    UL=padarray(UL,[kPadSz kPadSz],'symmetric','both');
end
km=1e-2;
m=km*msigma;
IPSF=rot90(PSF,2);
for k=1:NUMIT
    Gmsme=zeros(size(O));
    [~,WO]=atrous2(O,MaxScale);
    for l=1:MaxScale
        flags=double(abs(WO(:,:,l))>1e-8);
%         WO(:,:,l)=WO(:,:,l).*flags+abs(1-flags);
%         lnterm=log(abs(WO(:,:,l))).*flags;
        WO(:,:,l)=WO(:,:,l).*flags+m(:,:,l).*abs(1-flags);
        lnterm=log(abs(WO(:,:,l)./m(:,:,l))).*flags;
        if sum(isnan(lnterm(:)))>0
            error('NaN fault detected.')
        end
        Gmsme=Gmsme+imconv(A(:,:,l).*sign(WO(:,:,l)).*lnterm,...
            IWTF(:,:,l)).*msigma(:,:,l);
    end
    phi=imconv(I./imconv(O,PSF),IPSF)-alpha.*Gmsme./sigma;
    J(k)=max(max(abs(phi((1+kPadSz):(d1+kPadSz), (1+kPadSz):(d2+kPadSz))-1)));
    if k>2
        if J(k)>=J(k-1)
            disp('Convergence has been reached.')
            break
        end
    end
    O=O.*phi;
    O=O.*double(O>LL)+LL.*double(O<LL);
    if ~isempty(UL)
        O=O.*double(O<UL)+UL.*double(O>UL);
    end
end
J=J(1:k);
O=O((1+kPadSz):(d1+kPadSz), (1+kPadSz):(d2+kPadSz));
return

function M=mssupht(w,sigma,h)
%MSSUPHT multi-resolution support function, with 3-sigma hard threshold
% w - multiresolution wavelet data
% sigma - multiresolution noise levels of the wavelet data
%
[H,W,J]=size(w);
M=zeros(size(w));
K=sqrt(2*log10(H*W));
% g=fspecial('gauss',[H W],2);
for k=1:J
    M(:,:,k)=double(w(:,:,k) >= (K*sigma(:,:,k)));
    M(:,:,k)=imconv(M(:,:,k),h(:,:,k));
end
return
