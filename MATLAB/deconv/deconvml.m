function J=deconvml(varargin)
%DECONVML Maximum likelihood deconvolution.
%The iteration stops once number of iterations exceeds NUMIT or correction
%becomes smaller than EPS.
%
%Syntax:
%   DECONVML(I, PSF)
%   DECONVML(I, PSF, NUMIT)
%   DECONVML(I, PSF, NUMIT, LL)
%   DECONVML(I, PSF, NUMIT, LL, UL)
%   DECONVML(I, PSF, NUMIT, LL, UL, WEIGHT)
%   DECONVML(I, PSF, NUMIT, LL, UL, WEIGHT, READOUT)
%   DECONVML(I, PSF, NUMIT, LL, UL, WEIGHT, READOUT, DAMPN)
%   DECONVML(I, PSF, NUMIT, LL, UL, WEIGHT, READOUT, DAMPN, DAMPT)
%   DECONVML(I, PSF, NUMIT, LL, UL, WEIGHT, READOUT, DAMPN, DAMPT, ACCORDER)
%
%Input arguments:
%   I        is the observed data. I can be either a 2-D double matrix, a 3-D double matrix,
%   or a cell. If it is cell its elements are defined as:
%   I{1} is the observed data (2-D matrix or 3-D matrix).
%   I{2}, I{3}, ..., I{L-1} are the last estimate and estimates before the last
%   one.
%   I{L} is an auxiliary array generated in this program.
%
%   PSF      is the point spread function. PSF and I must be the same size.
%   If both PSF and I are 3-D matrix (or I is cell contains 3-D matrix) the
%   program will run in synthesis mode.
%
%   NUMIT    (optional) is the number of iterations. Default is 20.
%
%   LL       (optional) is the lower limit. Default is 0.
%
%   UL       (optional) is the upper limit. Default is null.
%
%   WEIGHT   (optional) is assigned to each pixel to reflect its recording
%   quality in the detector. A good pixel is assigned weight 1 while a bad pixel
%   is excluded by assigning it 0 weight. This matrix can be determined by
%   flat-field correction. Default is a unity matrix with the same size as the
%   observed data.
%
%   READOUT  (optional) is an array (or a value) corresponding to the
%   additive noise (e.g., background, foreground noise) and the variance 
%   of the read-out camera noise. READOUT has to be in the units of the
%   image. Default is 0.
%
%   DAMPN    (optional) is the damping parameter N. Default is 0.
%   DAMPT    (optional) is the damping threshold. Default is 0 (no damping).
%   AccOrder (optional) is the acceleration order. 
%
%ABOUT SYNTHESIS MODE
%   If the observed data is a collection from different detectors (with
%   different PSF) and all the observations correspond to the same object, the
%   program will run in synthesis mode. In this mode both the observed data I
%   and the PSF are 3-D matrix. The third dimension of each of them corresponds
%   to different observations.
%
%Return:
%   J        is either a 2-D array (if the input data I is a 2-D array) or a
%   cell (if the input data I is a cell).
%
%References:
%
%   [1] D.S.C. Biggs and M. Andrews, Acceleration of iterative image restoration
%   algorithms, Applied Optics, Vol. 36, No. 8, 1997.
%   [2] T.-P. Li and M. Wu, A Direct Restoration Method for Spectral and Image
%   Analysis, Astrophysics and Space Science, Vol. 206, 1993.
%   [3] R.J. Hanisch, R.L. White, and R.L. Gilliland. in "Deconvolution of Images 
%   and Spectra", Ed. P.A. Jansson, 2nd ed., Academic Press, CA, 1997.
%

% parse input arguments:
[J, PSF, NUMIT, LL, UL, WEIGHT, READOUT, ...
    DAMPN, DAMPT22, AccOrder, classI, ~, lenPSF, sizeJ, lenJ] = parse_inputs(varargin{:});

% convert PSF (point spread function)  into OTF (optical transfer function):
if lenPSF > 1
  OTF = PSF;
  for k = 1:lenPSF
    OTF(:,:,k) = psf2otf(PSF(:,:,k));
  end
else
  OTF = psf2otf(PSF);
end


% weighted observed data
% wI = max(WEIGHT.*(READOUT + J{1}),0);
wI = max(J{1},0);

scale = mean(real(ifft2(conj(OTF).*fft2(WEIGHT))),3) + sqrt(eps);
% imagesc(scale)
% clear WEIGHT;

% number of Taylor terms involved in the acceleration
AccTerms = AccOrder+1;

% generating taylor series coefficients matrix
TCM=zeros(AccOrder+1);
for k=0:AccOrder
    TCM(k+1,(k+1):(AccOrder+1)) = ...
        (1./factorial(0:(AccOrder-k)))*(1/factorial(k))*(-1)^(k);
end

curStep = AccTerms*any(J{lenJ}(:)~=0);
lambda = 0;

global VERBOSE
if VERBOSE >= 2
  hbar = waitbar(0,'RL deconvolving...');
  tML = tic;
end
for k = curStep + (1:NUMIT)
    % compute step size
    if k > AccTerms,
        lambda = (J{lenJ}(:,1).'*J{lenJ}(:,2))/...
            (J{lenJ}(:,2).'*J{lenJ}(:,2) +eps);
        lambda = max(min(lambda,1),0);% stability enforcement
        lambda = lambda^(1/AccOrder);
    end
    
    % make finit difference approximation
    lambdan = lambda.^(0:AccOrder)';
    TCML = TCM*lambdan;
    Y = zeros(sizeJ);
    for n = 0:AccOrder
        Y = Y+sum(TCML(n+1,:))*J{2+n};
    end
    
    % lower limit constraint
    Y = max(Y, LL);
    
    % upper limit constraint
    if ~isempty(UL), Y = min(Y, UL);end
    
    % compute damped RL iterative coefficient
    CC = damplucy(Y,OTF,DAMPN,DAMPT22,wI,READOUT,lenPSF);
    
    % update J
    for n=1:AccOrder
        J{AccOrder+3-n} = J{AccOrder+2-n};
    end
    
    % lower limit constraint
    J{2} = max(Y.*mean(real(ifft2(conj(OTF).*CC)),3)./scale, LL);
    
    % upper limit constraint
    if ~isempty(UL), J{2} = min(J{2}, UL);end
    
    J{AccOrder+3} = [J{2}(:)-Y(:) J{AccOrder+3}(:,1)];
    if VERBOSE >= 2
      tElapsed = toc(tML);
      tRemain = ((NUMIT-k)/k)*tElapsed;
      if mod(k,20)==0
        waitbar(k/NUMIT,hbar,['RL deconvolving, ',int2str(round(tRemain)),' seconds left...']);
      end
    end
end
if exist('hbar','var')==1,  close(hbar);end
if strcmp(classI, 'notcell')
    J=J{2};
end
return

function f=damplucy(Y,OTF,DAMPN,DAMPT22,wI,READOUT,lenPSF)
if lenPSF > 1
  ReBlurred = ...
    real(ifft2(fft2(padarray(Y,[0 0 lenPSF-1],'replicate','post')).*OTF));
else
  ReBlurred = real(ifft2(fft2(Y).*OTF));
end
ReBlurred = ReBlurred + READOUT;
ReBlurred(ReBlurred == 0) = eps;
AnEstim = wI./ReBlurred + eps;
if DAMPT22 == 0,% No Damping
  ImRatio = AnEstim;
else % Damping of the image relative to DAMPAR22 = (N*sigma)^2
  g = (wI.*log(AnEstim)+ ReBlurred - wI)./DAMPT22;
  g = min(g,1);
  G = (g.^(DAMPN-1)).*(DAMPN-(DAMPN-1)*g);
  ImRatio = 1 + G.*(AnEstim - 1);
end;
f = fft2(ImRatio);
return

function [J, PSF, NUMIT, LL, UL, WEIGHT, READOUT, ...
    DAMPN, DAMPT22, AccOrder, classI, sizePSF, lenPSF, sizeJ, lenJ]=...
    parse_inputs(varargin)
NUMIT=[];NUMIT_d=20;
LL=[];LL_d=0;
UL=[];
WEIGHT=[];WEIGHT_d=1;
READOUT=[];READOUT_d=0;
DAMPN=[];DAMPN_d=0;
DAMPT=[];DAMPT_d=0;
AccOrder=[];AccOrder_d=1;
narginchk(2,11);
J=varargin{1};
PSF=varargin{2};
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
        WEIGHT=varargin{6};
    case 7
        NUMIT=varargin{3};
        LL=varargin{4};
        UL=varargin{5};
        WEIGHT=varargin{6};
        READOUT=varargin{7};
    case 8
        NUMIT=varargin{3};
        LL=varargin{4};
        UL=varargin{5};
        WEIGHT=varargin{6};
        READOUT=varargin{7};
        DAMPN=varargin{8};
    case 9
        NUMIT=varargin{3};
        LL=varargin{4};
        UL=varargin{5};
        WEIGHT=varargin{6};
        READOUT=varargin{7};
        DAMPN=varargin{8};
        DAMPT=varargin{9};
    case 10
        NUMIT=varargin{3};
        LL=varargin{4};
        UL=varargin{5};
        WEIGHT=varargin{6};
        READOUT=varargin{7};
        DAMPN=varargin{8};
        DAMPT=varargin{9};
        AccOrder=varargin{10};
end
% check J
if ~iscell(J)
    classI='notcell';
    I=J;
    J={};
    J{1}=I;
    J{2}=mean(I,3);
    clear I
else
    classI='cell';
end
sizeI = size(J{1});
% check PSF
sizePSF=size(PSF);
if any(sizeI ~= sizePSF)
  error('PSF and I must be the same size.')
end
switch length(sizePSF)
  case 2
    sizeJ = sizePSF;
    lenPSF = 1;
  case 3
    disp('Synthesis mode activated.')
    sizeJ = sizePSF(1:2);
    lenPSF = sizePSF(3);
  otherwise
    error('Unsupported PSF.')
end
if isempty(AccOrder),   AccOrder=AccOrder_d;end
if isempty(DAMPT),      DAMPT=DAMPT_d;end
if isempty(DAMPN),      DAMPN=DAMPN_d;end
if isempty(LL),         LL=LL_d;end
if isempty(WEIGHT),     WEIGHT=WEIGHT_d;end
if isempty(READOUT),    READOUT=READOUT_d;end
if isempty(NUMIT),      NUMIT=NUMIT_d;end
if isscalar(WEIGHT),    WEIGHT=WEIGHT*ones(sizePSF);end
if isscalar(READOUT),   READOUT=READOUT*ones(sizePSF);end
if isscalar(LL),        LL=LL*ones(sizeJ);end
if isscalar(UL),        UL=UL*ones(sizeJ);end
if isscalar(DAMPT) && DAMPT~=0, DAMPT=DAMPT*ones(sizePSF);end
DAMPT22=(DAMPT.^2)/2;
AccOrder=max(round(AccOrder),0);
lenJ=3+AccOrder;

switch length(J)
    case lenJ
    case 2
        for k=3:(lenJ-1)
            J{k}=0;
        end
        J{lenJ}(prod(sizeJ), 2)=0;
    otherwise
        error('Wrong length of I.')
end
return
