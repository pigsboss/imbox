function varargout = dd(varargin)
%DD Direct demodulation method.
%
%Syntax:
%   dd(I,PSF)
%   dd(I,PSF,NUMIT)
%   dd(I,PSF,NUMIT,'ParamName1',ParamValue1,'ParamName2',ParamValue2,...)
%   [J, Res] = dd(...)
%
%Input arguments:
%   I                is the blurred image. It can be either a matrix or a cell.
%                    If I is a cell the first element of the cell is the blurred
%                    image (observed data) and the second element of the cell is
%                    the initial estimate of the true image.
%   PSF              is the point spread function. Its size must be the same as
%                    size of I.
%   NUMIT (optional) is the number of iterations (default is 20).
%
%Parameters:
%   'BACKGROUND' (alias 'bg', 'LowerLimit', 'll') is the background level used
%   as a constraint. The value can be one of the following:
%     1. a scalar specifies a constant background level explicitly, or
%     2. a matrix of the same size as I, which specifies a variable background
%     level, or
%     3. a string identifier 'manual', 'spline', 'intrinsic', 'gaussian' or
%     'clean' estimate the background during runtime.
%   Default is 0.
%
%   'UPPERLIMIT' (alias 'ul') is the upper limit constraint. The value can be
%   either a scalar (specifies a constant upper limit) or a matrix of the same
%   size as I.
%   Default is positive infinity.
%
%   'GAIN' is assigned to each pixel to reflect its recording
%   quality in the detector. A good pixel is assigned weight 1 while a bad pixel
%   is excluded by assigning it 0 weight. This matrix can be determined by
%   flat-field correction. Default is a unity matrix with the same size as the
%   observed data.
%
%   'FOLD' is number of folds of correlation transform applied to
%   both sides of modulation equation. Default is 0.
%
%   'WEIGHT' is an array (or a value) corresponding to the
%   additive noise (e.g., background, foreground noise) and the variance
%   of the read-out camera noise. This parameter has to be in the units of the
%   image and has to be either a scalar or a matrix of the same size as I.
%   Default is 0.
%
%   'AccelerateOrder' (alias 'AccOrder') is order of finite difference
%   approximation used to accelerate RL iterations. 0 is RL iteration without
%   acceleration. 1 is the default. Orders higher than 1 are available only for
%   tests. Don't use order higher than 6.
%
%   'DampThreshold' (alias 'dt') is the threshold of damped RL. Pixels beyond
%   threshold are taken as sharp structures of the image thus are not damped.
%   The rest of the image are damped during the RL deconvolution. This function
%   is useful to suppress the noise magnification of RL deconvolution.
%   Default is 0.
%
%   'DampOrder' (alias 'dn') is the order of damping effect. Orders higher than
%   10 have no significant difference on the result.
%   Default is 10.
%
%   'SIGMA' is the standard deviation of the additive noise in the observed data.
%   While 'Readout' parameter can be assigned to each pixel individually, this
%   parameter measures the deviation of the ensemble of all the pixels. This
%   parameter can be specified by:
%     1. a scalar explicitly, or
%     2. a string identifier 'manual', 'poisson', or 'auto'. Then the 'Sigma'
%     parameter will be determined during runtime according the the specified
%     identifier.
%   Default is 'auto'.
%
%   'ALPHA' is the parameter balances the data fidelity and the maximum entropy
%   constraint.
%   This parameter is determined automatically by default. If the 'Sigma'
%   parameter is specified properly there is no need to specify the 'Alpha'
%   parameter manually.
%
%   'MaximumScale' (alias 'MS', 'MaxScale') is the maximum scale of maximum
%   entropy constraint. The maximum entropy is applied to multiscale wavelet
%   coeffients of the image instead of the image itself.
%   Default is the maximum scale supported by the input, namely, if the size of I
%   is M x N the maximum scale is ceil(log2(min(M, N))).
%
%   'MaximumEntropy' (alias 'ME', 'MaxEnt') is the switch turns on and off the
%   multiscale maximum entropy constraint.
%   Default is false.
%
%   'VERBOSE' Verbose level of the runtime.
%     0 - quiet
%     1 - normal
%     2 - debug
%
%   'MODEL' (alias 'INITIAL') is the initial estimate of the true image.
%
%   All the parameters as well as all the identifiers are case-insensitive.
%
%ABOUT SYNTHESIS MODE
%   If the observed data is a collection from different detectors (with
%   different PSF) and all the observations correspond to the same object, the
%   program will run in synthesis mode. In this mode both the observed data I
%   and the PSF are 3-D matrix. The third dimension of each of them corresponds
%   to different observations.
%
%Returns:
%   J   (optional) is the reconstructed image.
%
%   Res (optional) is the residual image. The residual image is the difference of
%   the observed data and the reblurred image, i.e., the convolution of the
%   reconstructed image and the PSF.
%

global VERBOSE
[J,PSF,NUMIT,~,LL,UL,WEIGHT,READOUT,DAMPN,DAMPT,AccOrder,...
  sigma,MaxScale,alpha,...
  MaxEntOn,lenPSF,VERBOSE] = parse_inputs(varargin{:});
if MaxEntOn
%   J{1} = mean(J{1},3);
%   PSF = mean(PSF,3);
  J = deconvmsme(J,PSF,NUMIT,LL,UL,sigma,MaxScale,alpha);
  J{1} = ianscombe(J{1});
  J{2} = ianscombe(J{2});
else
  J = deconvml(J,PSF,NUMIT,LL,UL,WEIGHT,READOUT,DAMPN,DAMPT,AccOrder);
end
Res = J{1} - imconv(J{2},PSF);
if VERBOSE >= 1
  showImStat('Residual statistics',mean(Res,3),LL)
end
if nargout > 0
  varargout{1} = J{2};
  varargout{2} = Res;
else
  figure;imagesc(J{2});title('Reconstruction');axis xy
  if lenPSF > 1
    for k = 1:lenPSF
      figure;imagesc(J{1}(:,:,k));title(['Observed data ',int2str(k)]);axis xy
      figure;imagesc(PSF(:,:,k));title(['PSF ',int2str(k)]);axis xy
      figure;imagesc(Res(:,:,k));title(['Residual ',int2str(k)]);axis xy
    end
  else
    figure;imagesc(J{1});title('Observed data');axis xy
    figure;imagesc(PSF);title('PSF');axis xy
    figure;imagesc(mean(Res,3));title('Residual');axis xy
  end
end
return

function [datan,psfn] = nfoldcorr(data,psf,fold)
[~,~,lenPSF] = size(psf);
for k = 1:lenPSF
  datan = data(:,:,k);
  psfn = psf(:,:,k);
  sizePSF = size(psfn);
  for l = 1:fold
    datan = imconv(datan,psfn);
    psfn = imconv(psfn,psfn);
    kszPSF = kernelSize(psfn);
    if any(kszPSF==sizePSF)
      disp([int2str(l),' fold correlated PSF coefficients exceed boundary.'])
    end
  end
end
return

function [J,PSF,NUMIT,fold,LL,UL,WEIGHT,READOUT,DAMPN,DAMPT,AccOrder,...
  sigma,MaxScale,alpha,...
  MaxEntOn,lenPSF,VERBOSE] = ...
  parse_inputs(varargin)
NUMIT=[];NUMIT_d=20;
LL=[];LL_d=0;
UL=[];
fold=[];fold_d=0;
dsample=[];dsample_d=0;
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
MaxEntOn = []; MaxEntOn_d = false;
narginchk(2,100);
J=cell(1,2);
J{1}=varargin{1};
PSF=varargin{2};
sizeI = size(J{1}); % size of observed data matrix
sizePSF = size(PSF); % size of PSF matrix
sizeJ = sizeI(1:2); % size of image matrix
if nargin>=3 && isnumeric(varargin{3}), NUMIT = varargin{3};end
k = 3;
while k < nargin
%   disp(varargin{k})
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
      case {'MAXIMUMENTROPY','MAXENT','ME'}
        k = k+1;
        MaxEntOn = logical(varargin{k});
      case {'MODEL','INITIAL'}
        k = k+1;
        J{2} = varargin{k};
      case {'VERBOSE','V'}
        k = k+1;
        VERBOSE = varargin{k};
      case {'FOLD','F'}
        k = k+1;
        fold = varargin{k};
      case {'DOWNSAMPLE','DS'}
        k = k+1;
        dsample = varargin{k};
    end
    k = k+1;
  else
    error('Syntax error. Please run "help dd" to refer to the document.')
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
if isempty(fold),       fold=fold_d;end
if isempty(dsample),    dsample=dsample_d;end
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
  disp(['Size of observed data: ',int2str(sizeI(1)),' x ',int2str(sizeI(2))]);
  disp(['Size of PSF:           ',int2str(sizePSF(1)),' x ',int2str(sizePSF(2))]);
  if length(sizePSF) > 2
    disp('Multiple observations found.')
  end
%   disp( '++++++++++++++++++++++++++++++++++')
end

[~,~,lenPSF] = size(PSF);
PSF = padpsf(PSF,sizeI);
PSF_unfold = PSF;

if fold>0
    [J{1},PSF] = nfoldcorr(J{1},PSF_unfold,fold);
end

% check LL
if ischar(LL)
    LL = bgrd(mean(J{1},3),mean(PSF,3),LL);
    if VERBOSE >= 1
      disp(['Estimated mean background level: ',num2str(mean(LL(:)))])
    end
elseif isscalar(LL)
    LL = LL*ones(sizeJ);
end

if dsample>0
    J{1} = imdsample(J{1},dsample);
    PSF = imdsample(PSF,dsample);
end

sizeI = size(J{1}); % size of observed data matrix
sizeJ = sizeI(1:2); % size of image matrix

% check J
if isempty(J{2})
  J{2} = mean(J{1},3);
else
    if fold>0
        [J{2},~] = nfoldcorr(J{2},PSF_unfold,fold);
    end
    if dsample>0
        J{2} = imdsample(J{2},dsample);
    end
end

% check other aux matrix
X={LL,UL,WEIGHT,READOUT};
for k = 1:length(X)
    if ~isempty(X{k})
        if fold>0
            [X{k},~] = nfoldcorr(X{k},PSF_unfold,fold);
        end
        if dsample > 0
            X{k} = imdsample(X{k},dsample);
        end
    end
end
[LL,UL,WEIGHT,READOUT] = X{:};
clear X

if any(LL < 0)
    error('Background must be non-negative.')
end

% MSME parameters:
if isempty(MaxEntOn),   MaxEntOn=MaxEntOn_d;end
if MaxEntOn
  J{1} = anscombe(mean(J{1},3));
  J{2} = anscombe(J{2});
  LL = anscombe(LL);
  UL = anscombe(UL);
  sigma_d = 'POISS';
  PSF = mean(PSF,3);
  lenPSF = 1;
  if isempty(sigma),    sigma = sigma_d;end
  if ischar(sigma)
    switch upper(sigma)
      case {'POISSON','POISS','P'}
        sigma = 1;
      case {'AUTO','A'}
        sigma = sigmaclipping(J{1} - sdwt2(J{1},1));
      case {'MANUAL','MAN','M'}
        figure;imagesc(J{1});drawnow;
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
end

if VERBOSE >= 1
  if lenPSF > 1
    for k = 1:lenPSF
      showImStat(['Observed data ',int2str(k),' statistics'],J{1}(:,:,k),LL)
    end
  else
    showImStat('Observed data statistics',J{1},LL)
  end
end
if VERBOSE >= 2
  %TODO:  show all the parameters
  figure;imagesc(LL);title('Background');axis xy
%   figure;imagesc(UL);title('Background');axis xy
end
return
