function [sigmamean,sigmastd,sigmatab] = msgsigma(NUMIT,scale,MODE)
%MSGSIGMA sigma of multi-scale Laplacians of Gaussians of noissy image.
narginchk(0,3)
NUMIT_d = 32;
scale_d = (1:20);
N = 256;
MODE_d = 'generate';
rng('shuffle')
switch nargin
  case 0
    NUMIT = NUMIT_d;
    scale = scale_d;
    MODE = MODE_d;
  case 1
    if ischar(NUMIT)
      MODE = NUMIT;
      NUMIT = NUMIT_d;
    else
      MODE = MODE_d;
    end
    scale = scale_d;
  case 2
    if isempty(NUMIT)
      MODE = 'collect';
    else
      MODE = 'generate';
    end
  case 3
    if isempty(NUMIT)
      NUMIT = NUMIT_d;
    end
    if isempty(scale)
      scale = scale_d;
    end
    if isempty(MODE)
      MODE = MODE_d;
    end
end
scale = scale.^2;
switch MODE
  case {'generate','gen','run'}
    disp('generating MC...')
    NS = length(scale);
    sigmas = zeros(NUMIT,NS);
    ofile = ['msgsigma.',int2str(NUMIT),'x',int2str(NS),'.',...
      int2str(feature('getpid'))];
    parfor l = 1:NUMIT
    f = normrnd(100,1,N);
    L = laplams(f,scale);
      for k = 1:NS
        sigmas(l,k) = std(reshape(L(:,:,k),1,[]));
      end
    end
    sigmamean = mean(sigmas,1);
    sigmastd = std(sigmas,0,1);
    save([ofile,'.',datetimefilename,'.mat'])
  case {'collect'}
    disp('collecting MC results...')
    files = dir('msgsigma*.mat');
    numfiles = length(files);
    NS = length(scale);
    sigmamean = zeros(1,NS);
    sigmastd = zeros(1,NS);
    sigma2mean = zeros(1,NS);
    totalmc = 0;
    sigmas = cell(numfiles,1);
    for k = 1:numfiles
      mat = load(files(k).name);
      sigmas{k} = mat.sigmas;
      [mcinmat,~] = size(mat.sigmas);
      sigmainmat = mean(mat.sigmas,1);
      sigmamean = (sigmamean*totalmc +...
        interp1(mat.scale,sigmainmat,scale)*mcinmat) / (totalmc + mcinmat);
      sigma2mean = (sigma2mean*totalmc +...
        interp1(mat.scale,sigmainmat,scale).^2*mcinmat) / (totalmc + mcinmat);
      sigmastd = sqrt(sigma2mean - sigmamean.^2);
      totalmc = totalmc + mcinmat;
      linefitopt = linefit(log(sqrt(scale)),log(sigmamean));
      disp(['MC number: ',int2str(totalmc),...
        '; a = ',num2str(linefitopt(1)),...
        '; b = ',num2str(linefitopt(2)),...
        '; r = ',num2str(linefitopt(3)),...
        '; rms = ',num2str(linefitopt(4))])
    end
    sigmatab = cell2mat(sigmas);
end
return