function [J,chi2,MinTar]=deconvmsme(varargin)
%DECONVMSME multiscale maximum entropy deconvolution
%
%Syntax:
% deconvmsme(I, PSF)
% deconvmsme(I, PSF, NUMIT)
% deconvmsme(I, PSF, NUMIT, LL)
% deconvmsme(I, PSF, NUMIT, LL, UL)
% deconvmsme(I, PSF, NUMIT, LL, UL, sigma)
% deconvmsme(I, PSF, NUMIT, LL, UL, sigma, MaxScale)
% deconvmsme(I, PSF, NUMIT, LL, UL, sigma, MaxScale, alpha)
%
%Input arguments:
% I - observed data, with Gaussian noise
% PSF - point spread function
% LL (optional, default: []) - background constraint
% NUMIT (optional, default: 20) - maximum number of iterations
% MaxScale (optional, default: log2(size_of_I)-2) - maximum scale of wavelet transforms
% sigma (optional, default: semi-automatically determined)- standard deviation of noise in observed data I
% alpha (optional, default: automatically determined) - parameter balancing the entropy and the chi-square terms
% 
%Reference:
% [1]. E. Pantin and J.-L. Starck, Deconvolution of astronomical images using
% the multiscale maximum entropy method, Astronomy and Astrophysics Supp.
% Series, Vol. 118, p575-p585, 1996.
% [2]. T.-P. Li and M. Wu, A Direct Restoration Method for Spectral and 
% Image Analysis, Astrophysics and Space Science. 1993, p91-p102.
% [3]. M. Raydan and B.F. Svaiter, Relaxed Steepest Descent and 
% Cauchy-Barzilai-Borwein Method, Computational Optimization and Applications,
% Vol. 21, 2002, p155-p167.
%

% parsing input arguments:
[J, PSF, NUMIT, LL, UL, sigma, MaxScale, alpha,...
    sizeI, classI, sigma2] = ...
    parse_inputs(varargin{:});

MinTar = zeros(NUMIT,1); % initialize minimizing target
chi2 = MinTar; % initialize chi-square term
kPadSz = 2.^((1:MaxScale) + 1); % kernel-padding size for wavelet transform
BSF = bsplinefilter(sizeI); % b-spline filter (scaling function)
BSF = BSF(:,:,1:(MaxScale+1));
[~,w] = sdwt2(J{1},MaxScale);
msigma = sigma*mssigma(MaxScale); % multi-resolution sigma
M = mssupht(w,msigma,BSF); % multiresolution support
A = 1-M;
km = eps;
m = km*msigma;
OTF = psf2otf(PSF);
WTFconj = cell(MaxScale,1);
for k = 1:MaxScale
    WTFconj{k} = conj(psf2otf(padarray(BSF(:,:,k)-BSF(:,:,k+1),...
        [kPadSz(k),kPadSz(k)],0,'both')));
end
clear BSF DDF w M km PSF

% initialize random step size
rng('shuffle')
% u = rand(NUMIT,1);
MinTarUL = sqrt(numel(J{1})) + (numel(J{1})).^(1/4);
MinTarLL = sqrt(numel(J{1})) - (numel(J{1})).^(1/4);
TotalFlux = sum(ianscombe(J{1}(:)));
global VERBOSE
if VERBOSE >= 1
  hbar = waitbar(0,'Multiscale ME deconvolving...');
  tMSME = tic;
end
for k = 1:NUMIT
    Reblurred = real(ifft2(fft2(J{2}).*OTF));
    delta = J{1} - Reblurred;
    chi2(k) = sum(delta(:).^2)/sigma2;
    Gchi2 = -real(ifft2(fft2(delta).*conj(OTF)))/sigma2;% gradient
    Gmsme = zeros(sizeI);
    Sms = 0;
    [~,WO] = sdwt2(J{2});
    for l = 1:MaxScale
        flags = double(abs(WO(:,:,l)) > eps);
        WO(:,:,l) = WO(:,:,l).*flags+m(l)*abs(1-flags);
        lnterm = log(abs(WO(:,:,l)/m(l))).*flags;
        if sum(isnan(lnterm(:)))>0
            error('Found NaN!')
        end
        Sms = Sms+sum(sum(A(:,:,l).*(abs(WO(:,:,l))-m(l)-abs(WO(:,:,l)).*lnterm)*msigma(l)));
        WWO = padarray(A(:,:,l).*sign(WO(:,:,l)).*lnterm, [kPadSz(l), kPadSz(l)], 'symmetric', 'both');
        WWO = real(ifft2(fft2(WWO).*WTFconj{l}))*msigma(l);
        Gmsme = Gmsme+WWO((1+kPadSz(l)):(sizeI(1)+kPadSz(l)),...
            (1+kPadSz(l)):(sizeI(2)+kPadSz(l)));
    end
    G = Gchi2+Gmsme*alpha/sigma;
    MinTar(k) = chi2(k)-Sms*alpha/sigma;
    Q = real(ifft2(fft2(G).*OTF));% step size
    gamma = (Q(:)'*Reblurred(:) - J{1}(:)'*Q(:))/sum(Q(:).^2);
    J{2} = J{2}-gamma*G; % steepest descent
%     J{2} = max(J{2}-2*u(k)*gamma*G, LL); % steepest descent
    if ~isempty(LL), J{2} = max(J{2}, LL);end
    if ~isempty(UL), J{2} = min(J{2}, UL);end
    % total flux constraint
%     J{2} = ianscombe(J{2});
%     J{2} = anscombe((TotalFlux/sum(J{2}(:)))*J{2});
    if k>1 && MinTar(k)<=MinTarUL
        if (abs(MinTar(k)-MinTar(k-1))/MinTar(k-1))<1e-3 || MinTar(k)<MinTarLL
            disp('Converges.');
            break
        end
    end
    if VERBOSE >= 1
      tElapsed = toc(tMSME);
      tRemain = ((NUMIT-k)/k)*tElapsed;
      if mod(k,2)==0
        waitbar(k/NUMIT,hbar,['Multiscale ME deconvolving, ',int2str(round(tRemain)),' seconds left...']);
      end
    end
end
if exist('hbar','var')==1,  close(hbar);end
chi2 = chi2(1:k);
MinTar = MinTar(1:k);
if strcmp(classI, 'notcell')
    J=J{2};
end
return

function M=mssupht(w,sigma,h)
%MSSUPHT multi-resolution support function, with 3-sigma hard threshold
% w - multiresolution wavelet data
% sigma - multiresolution noise levels of the wavelet data
%
[H,W,J]=size(w);
if length(sigma)~=J
    error('Maximum scales of the wavelet data and its noise level do not match.')
end
M=zeros(size(w));
K=sqrt(2*log10(H*W));
% g=fspecial('gauss',[H W],2);
for k=1:J
    M(:,:,k)=double(w(:,:,k) >= (K*sigma(k)));
    M(:,:,k)=imconv(M(:,:,k),h(:,:,k));
end
return

function [J, PSF, NUMIT, LL, UL, sigma, MaxScale, alpha,...
    sizeI, classI, sigma2]=...
    parse_inputs(varargin)
NUMIT=[];NUMIT_d=20;
LL=[];
UL=[];
MaxScale=[];
alpha=[];
sigma=[];
narginchk(2,9)
J = varargin{1};
PSF = varargin{2};
switch nargin
    case 3
        NUMIT=varargin{3};
    case 4
        NUMIT=varargin{3};
        LL=varargin{4};
    case 5
        NUMIT=varargin{3};
        LL=varargin{4};
        UL=varargin{5};
    case 6
        NUMIT=varargin{3};
        LL=varargin{4};
        UL=varargin{5};
        sigma=varargin{6};
    case 7
        NUMIT=varargin{3};
        LL=varargin{4};
        UL=varargin{5};
        sigma=varargin{6};
        MaxScale=varargin{7};
    case 8
        NUMIT=varargin{3};
        LL=varargin{4};
        UL=varargin{5};
        sigma=varargin{6};
        MaxScale=varargin{7};
        alpha=varargin{8};
end
if ischar(J),         J=imread(J);J=mean(double(J),3);end
if ischar(PSF),       PSF=imread(PSF);PSF=mean(double(PSF),3);end
if ischar(LL),        LL=str2double(LL);end
if ischar(sigma),     sigma=str2double(sigma);end
if ischar(MaxScale),  MaxScale=str2double(MaxScale);end
if ischar(NUMIT),     NUMIT=str2double(NUMIT);end
if isempty(NUMIT),    NUMIT = NUMIT_d;end
NUMIT = round(NUMIT);
% check J
if ~iscell(J)
    classI='notcell';
    I=J;
    J={};
    J{1}=I;
    J{2}=I;
    clear I
else
    classI='cell';
end
sizeI = size(J{1});
% check PSF
sizePSF = size(PSF);
if any(sizeI ~= sizePSF)
  error('PSF and I must be the same size.')
end

if isempty(LL) && isempty(sigma)
    disp('Please select a region to estimate the background and the noise level.')
    [~,LL,sigma] = imflux(J{1},[]);
    disp(['The background level estimate is ',num2str(LL)])
    disp(['The noise level estimate is ',num2str(sigma)])
end
% if isempty(LL)
%     imagesc(J{1});drawnow;disp('Please select a region to estimate the background level.')
%     [~,LL] = imflux;disp(['The background level estimate is ',num2str(LL)])
% end
if isempty(sigma)
    imagesc(J{1});drawnow;disp('Please select a region to estimate the noise level.')
    I2=imcrop;
    sigma=sigmaclipping(I2);disp(['The noise level estimate is ',num2str(sigma)])
end
if isempty(MaxScale), MaxScale=lastpow2(min(size(J{1})));end
if isempty(alpha),    alpha=0.5*max(PSF(:))/sigma;end
if isscalar(LL),      LL=LL*ones(sizePSF);end
if isscalar(UL),      UL=UL*ones(sizePSF);end
sigma2 = sigma^2;
MaxScale = min(MaxScale,lastpow2(min(size(J{1}))));
return

% function [O, chi2, MinTar] = corechi2msme(O,I,OTF,sigma,sigma2,msigma,...
%     m,A,WTFconj,kPadSz,MaxScale,alpha,u,sizeI)
% Reblurred = real(ifft2(fft2(O).*OTF));
% delta = I - Reblurred;
% chi2 = sum(delta(:).^2)/sigma2;
% Gchi2 = -real(ifft2(fft2(delta).*conj(OTF)))/sigma2;% gradient
% Gmsme = zeros(sizeI);
% Sms = 0;
% [~,WO] = sdwt2(O);
% for l = 1:MaxScale
%     flags = double(abs(WO(:,:,l)) > eps);
%     WO(:,:,l) = WO(:,:,l).*flags+m(l)*abs(1-flags);
%     lnterm = log(abs(WO(:,:,l)/m(l))).*flags;
% %     WO(:,:,l) = (abs(WO(:,:,l))+eps).*sign(WO(:,:,l));
% %     lnterm = log(abs(WO(:,:,l)/m(l)));
%     if sum(isnan(lnterm(:)))>0
%         error('Found NaN!')
%     end
%     Sms = Sms+sum(sum(A(:,:,l).*(abs(WO(:,:,l))-m(l)-abs(WO(:,:,l)).*lnterm)*msigma(l)));
%     WWO = padarray(A(:,:,l).*sign(WO(:,:,l)).*lnterm, [kPadSz(l), kPadSz(l)], 'symmetric', 'both');
%     WWO = real(ifft2(fft2(WWO).*WTFconj{l}))*msigma(l);
%     Gmsme = Gmsme+WWO((1+kPadSz(l)):(sizeI(1)+kPadSz(l)),...
%         (1+kPadSz(l)):(sizeI(2)+kPadSz(l)));
% end
% G = Gchi2+Gmsme*alpha/sigma;
% MinTar = chi2-Sms*alpha/sigma;
% Q = real(ifft2(fft2(G).*OTF));% step size
% gamma = (Q(:)'*Reblurred(:) - I(:)'*Q(:))/sum(Q(:).^2);
% O = O-2*u*gamma*G; % steepest descent
% return
