function varargout = imclean(varargin)
%IMCLEAN Improved Iterative Maximum likelihood CLEAN algorithm.
global VERBOSE
[J,PSF,NUMIT,NUMML,NUMMSME,LL,UL,...
  WEIGHT,READOUT,DAMPN,DAMPT,AccOrder,...
  sigma,MaxScale,alpha,...
  SmoothScale,~,VERBOSE] = ...
  parse_inputs(varargin{:});
if VERBOSE > 0
  hf = figure('Name','IMCLEAN');
  subplot(2,3,1);imagesc(J{1});axis image;axis xy;
  title('Observed data');xlabel('x pixels');ylabel('y pixels');drawnow
  subplot(2,3,2);imagesc(PSF);axis image;axis xy;
  title('PSF');xlabel('x pixels');ylabel('y pixels');drawnow
end
for k = 1:NUMIT
  %   J = DECONVML(I, PSF, NUMIT, LL, UL, WEIGHT, READOUT, DAMPN, DAMPT, ACCORDER)
  JML = deconvml(J{1},PSF,NUMML,J{3},UL,WEIGHT,READOUT,DAMPN,DAMPT,AccOrder);
  [Btab,Bimg] = blobex(JML,SmoothScale);
  if VERBOSE > 0
      printbtab(Btab(1:3,:),'sources');
  end
  Stab = blobfilter(Btab);
  Simg = tab2map(Stab,sizeI);
  J{2} = Simg;
  %   J = deconvmsme(J,PSF,NUMIT,LL,UL,sigma,MaxScale,alpha)
  JMSME = deconvml(anscombe(J{1}-J{2}),PSF,NUMMSME,...
      anscombe(LL),anscombe(UL),sigma,MaxScale,alpha);
  J{3} = ianscombe(JMSME);
  Res = J{1} - imconv(J{2}+J{3},PSF);
  if VERBOSE > 0
    figure(hf);
    subplot(2,3,3);imagesc(gsmooth(Bimg,2));axis image;axis xy;
    title('Resolved blobs');xlabel('x pixels');ylabel('y pixels');drawnow
    subplot(2,3,4);imagesc(gsmooth(J{2},2));axis image;axis xy;
    title('Resolved point sources');xlabel('x pixels');ylabel('y pixels');drawnow
    subplot(2,3,5);imagesc(J{3});axis image;axis xy;
    title('Resolved background');xlabel('x pixels');ylabel('y pixels');drawnow
    subplot(2,3,6);imagesc(Res);axis image;axis xy;
    title('Residual map');xlabel('x pixels');ylabel('y pixels');drawnow
  end
end
% if exist('hf','var')==1, close(hf);end
varargout{1} = J;
varargout{2} = Res;
return

function [J,PSF,NUMIT,NUMML,NUMMSME,LL,UL,...
    WEIGHT,READOUT,DAMPN,DAMPT,AccOrder,...
    sigma,MaxScale,alpha,...
    SmoothScale,lenPSF,VERBOSE] = ...
    parse_inputs(varargin)
NUMIT=[];NUMIT_d=20;
NUMML=[];NUMML_d=50; % number of ML (RL) iterations per IMCLEAN iteration.
NUMMSME=[];NUMMSME_d=20; % number of MSME iterations per IMCLEAN iteration.
LL=[];LL_d=0; % initial background constraint
UL=[];
SmoothScale=[]; SmoothScale_d=10;
VERBOSE = []; VERBOSE_d = 1; VERBOSE_max = 2;
WEIGHT=[];WEIGHT_d=1;
READOUT=[];READOUT_d=0;
DAMPN=[];DAMPN_d=0;
DAMPT=[];DAMPT_d=0;
AccOrder=[];AccOrder_d=1;
sigma = [];
% sigma_d = 'Auto';
alpha = [];
MaxScale = [];
narginchk(2,100);
J=cell(1,3); % J{1}, J{2} and J{3} are observed data, resolved point sources and resolved background.
J{1}=varargin{1};
PSF=varargin{2};
sizeI = size(J{1}); % size of observed data matrix
sizePSF = size(PSF); % size of PSF matrix
sizeJ = sizeI(1:2); % size of image matrix
if nargin>=3 && isnumeric(varargin{3}), NUMIT = varargin{3};end
k = 3;
while k < nargin
  if isscalar(varargin{k}) || ischar(varargin{k})
    switch upper(varargin{k})
      case {'BACKGROUND','LL','BG'}
        k = k+1;
        LL = varargin{k};
      case {'UL','UPPERLIMIT'}
        k = k+1;
        UL = varargin{k};
      case {'GAIN'}
        k = k+1;
        WEIGHT = varargin{k};
      case {'WEIGHT'}
        k = k+1;
        READOUT = varargin{k};
      case {'ACCELERATEORDER','ACCORDER'}
        k = k+1;
        AccOrder = varargin{k};
      case {'DAMPTHRESHOLD','DAMPT','DT'}
        k = k+1;
        DAMPT = varargin{k};
      case {'DAMPORDER','DAMPN','DN'}
        k = k+1;
        DAMPN = varargin{k};
      case {'SIGMA'}
        k = k+1;
        sigma = varargin{k};
      case {'ALPHA'}
        k = k+1;
        alpha = varargin{k};
      case {'MAXIMUMSCALE','MAXSCALE','MS'}
        k = k+1;
        MaxScale = varargin{k};
      case {'MODEL','INITIAL'}
        k = k+1;
        J{2} = varargin{k};
      case {'NUMML','ML'}
        k = k+1;
        NUMML = varargin{k};
      case {'NUMMSME','MSME','ME'}
        k = k+1;
        NUMMSME = varargin{k};
      case {'SMOOTHSCALE','SMOOTH'}
        k = k+1;
        SmoothScale = varargin{k};
      case {'VERBOSE','V'}
        k = k+1;
        VERBOSE = varargin{k};
    end
    k = k+1;
  else
    error('Syntax error. Please run "help imclean" to refer to the document.')
  end
end

if isempty(VERBOSE),    VERBOSE=VERBOSE_d;end

% ML parameters:
if isempty(AccOrder),   AccOrder=AccOrder_d;end
if isempty(DAMPT),      DAMPT=DAMPT_d;end
if isempty(DAMPN),      DAMPN=DAMPN_d;end
if isempty(LL),         LL=LL_d;end
if isempty(WEIGHT),     WEIGHT=WEIGHT_d;end
if isempty(READOUT),    READOUT=READOUT_d;end
if isempty(NUMIT),      NUMIT=NUMIT_d;end
if isempty(NUMML),      NUMML=NUMML_d;end
if isempty(SmoothScale),SmoothScale=SmoothScale_d;end
if isempty(NUMMSME),    NUMMSME=NUMMSME_d;end
if isscalar(WEIGHT),    WEIGHT=WEIGHT*ones(sizeI);end
if isscalar(READOUT),   READOUT=READOUT*ones(sizeI);end
if isscalar(UL),        UL=UL*ones(sizeJ);end
if isscalar(DAMPT) && DAMPT~=0, DAMPT=DAMPT*ones(sizeI);end
AccOrder=max(round(AccOrder),0);
VERBOSE = max(min(floor(VERBOSE),VERBOSE_max),0);

if any(sizeI < sizePSF)
    error('Size of PSF must not exceed size of I.')
end

if length(sizeI) ~= length(sizePSF)
    error('PSF and I matrices must be the same dimension.')
end

if VERBOSE >= 1
%   disp( '++++++++++++++++++++++++++++++++++')
  disp(['Size of observed data (by pixels): ',int2str(sizeI(1)),' x ',int2str(sizeI(2))]);
  disp(['Size of PSF (by pixels):           ',int2str(sizePSF(1)),' x ',int2str(sizePSF(2))]);
  if length(sizePSF) > 2
    disp('Multiple observations found.')
  end
%   disp( '++++++++++++++++++++++++++++++++++')
end

[~,~,lenPSF] = size(PSF);
if any(sizeI ~= sizePSF)
  PSFext = zeros(sizeI);
  for k = 1:lenPSF
    sizePSF = size(PSF);
    PSF(:,:,k) = ifftshift(PSF(:,:,k));
    PSFtmp = [PSF(1:round(sizePSF(1)/2),:,k);...
        zeros(sizeI(1)-sizePSF(1),sizePSF(2));...
        PSF((1+round(sizePSF(1)/2)):sizePSF(1),:,k)];
    sizePSF = size(PSFtmp);
    PSFext(:,:,k) = [PSFtmp(:,1:round(sizePSF(2)/2)),...
        zeros(sizePSF(1),sizeI(2)-sizePSF(2)),...
        PSFtmp(:,(1+round(sizePSF(2)/2)):sizePSF(2))];
%     sizePSF = size(PSF);
    PSFext(:,:,k) = fftshift(PSFext(:,:,k));
  end
  PSF = PSFext;
  clear PSFext PSFtmp
end

if ischar(LL)
    LL = bgrd(mean(J{1},3),mean(PSF,3),LL);
elseif isscalar(LL)
    LL = LL*ones(sizeJ);
end

% check J
if isempty(J{2})
  J{2} = zeros(sizeI);
end
if isempty(J{3})
  J{3} = LL;
end

if any(LL < 0)
    error('Background must be non-negative.')
end

% MSME parameters:
Ipoiss = anscombe(mean(J{1},3));
sigma_d = 'POISS';
PSF = mean(PSF,3);
lenPSF = 1;
if isempty(sigma),    sigma = sigma_d;end
if ischar(sigma)
switch upper(sigma)
  case {'POISSON','POISS','P'}
    sigma = 1;
  case {'AUTO','A'}
    sigma = sigmaclipping(Ipoiss - sdwt2(Ipoiss,1));
  case {'MANUAL','MAN','M'}
    figure;imagesc(Ipoiss);drawnow;
    disp('Please select a region to estimate the noise level.')
    I2=imcrop;
    sigma=sigmaclipping(I2 - sdwt2(I2,1));
    disp(['The noise level estimate is ',num2str(sigma)])
  otherwise
    error('Unsupported parameter SIGMA. Please run "help dd" to refer to the document.')
end
end
if sigma==0,          alpha=0;end
if isempty(alpha),    alpha = 0.5*max(PSF(:))/sigma;end
if isempty(MaxScale), MaxScale = lastpow2(min(sizeJ));end

% show input data infomation
if VERBOSE >= 1
  if lenPSF > 1
    for k = 1:lenPSF
      showImStat(['Observed data ',int2str(k),' statistics'],J{1}(:,:,k),LL)
    end
  else
    showImStat('Observed data statistics',J{1},LL)
  end
end
return

function Stab = blobfilter(Btab)
%BLOBFILTER Blob filter rejects blobs which are not significant according to the
%input null distribution and returns blobs which are more likely to be sources
%than those rejected.

Stab = Btab;
return

function I = tab2map(T,sizeI)
[NUMROW,lenT] = size(T);
if NUMROW < 3
    error('Input table must contain 3 rows at least.')
end
f = T(1,:);
x = min(max(round(T(2,:)),1),sizeI(2));
y = min(max(round(T(3,:)),1),sizeI(1));
I = zeros(sizeI);
for n = 1:lenT
    I(y(n),x(n))=f(n);
end
return
