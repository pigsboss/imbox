function varargout = inlmeans(varargin)
%INLMEANS observation-independent non-local Means Method
%
%SYNTAX
%  inlmeans(img,parameters,...)
%  img is image filename or image matrix.
%
%PARAMETERS
%  radius, (alias: r)        is the radius of Gaussian filter weighing the
%  Euclidean distance.
%
%  alpha, (alias: a)         is the standard deviation of Gaussian filter
%  weighing the variance of corresponding pixels in two neighborhood.
%
%  h                         is the standard deviation of Gaussian filter
%  weighing similarities of all other pixels between the given pixel.
%
%  repeats, (alias: reps, m) is the number of denoising repeats.
%
%  mc                        (boolean) specifies if MC method is used.
%
%  rotate, (alias: rot, t)   (boolean) specifies if rotation is used.
%
%  anime, (alias: movie, v)  (string) is the filename prefix of optional video
%  clips.
%
%  output, (alias: o)        (string) is the filename prefix of optional png
%  pictures.
%
%  colormap, (alias: cmap)   is colormap used for convert double-precision float
%  matrices to video frames.
%
%  pngcolormap, (alias: pngcmap)  is colormap used for convert double-precision 
%  float matrices to png pictures.
%
%  rpf                       is number of repeats per video frame.
%
%  fps                       is the frame rate of output (optional) video clip.
%  
%References:
% [1] Buades, et al., On image denoising methods, Int J Comput Vis, 76: 123 -
%     139, 2008.
% [2] Buades and Morel, A non-local algorithm for image denoising, in IEEE
%     international conference on computer vision and pattern recognition, 2005.
[img,rimg,pradius,~,h2,reps,wt,rot,ofile,anime,rpf,sizeI,mc,fps,...
  cmap,pngcmap] = parse_inputs(varargin{:});
rng('shuffle')
NumPix = numel(img);
N1 = sizeI(1);
simg = img; % sum image
z = ones(sizeI); % weights
beta = zeros(ceil(reps/rpf),1);
if ~isempty(ofile)
  ofile = [ofile,'.',int2str(feature('getpid'))];
end
if ~isempty(anime)
  anime = [anime,'.',int2str(feature('getpid')),'.',datetimefilename];
  dovid=VideoWriter([anime,'.denoise.avi'],'Uncompressed AVI');
  dovid.FrameRate = fps;
  open(dovid)
  zovid=VideoWriter([anime,'.simap.avi'],'Uncompressed AVI');
  zovid.FrameRate = fps;
  open(zovid)
end
f = rpf;
fn = 1;
[lencmap,~] = size(cmap);
for  l = 1:reps
  if mc == true
    n1 = mod(ceil(rand(1)*NumPix)-1,NumPix)+1;
  else
    n1 = l;
  end
  r1 = mod(n1,N1);
  c1 = (n1-r1)/N1 + 1;
  p1 = pxnbh(rimg,[r1,c1],pradius);
  s = zeros(sizeI);
  s = max(s,exp(-1*awedist(p1,rimg,wt).^2 / h2));
  if rot > 0
    for rotn = 1:rot
      p1r = pxnbh(rimg,[r1,c1],pradius,pi*2*rotn/(rot+1));
      s = max(s,exp(-1*awedist(p1r,rimg,wt).^2 / h2));
    end
  end
  s(n1) = 0;
  simg(n1) = simg(n1) + sum(img(:).*s(:));
  z(n1) = z(n1) + sum(s(:));
  simg = simg + img(n1)*s;
  z = z + s;
  dimg = simg./z; % denoised image
  if f == rpf
    if ~isempty(anime)
      frame = im2frame(uint8(lencmap*...
        (dimg - min(dimg(:)))/(max(dimg(:)) - min(dimg(:)))),cmap);
      writeVideo(dovid,frame);
      frame = im2frame(uint8(lencmap*(log10(z))/max(log(z(:)))),cmap);
      writeVideo(zovid,frame);
    end
    if ~isempty(ofile)
      save([ofile,'.snapshot.mat']);
      im2png(dimg,pngcmap,[ofile,'.snapshot.denoise.png']);
      im2png(log10(z),pngcmap,[ofile,'.snapshot.simap.png']);
    end
    beta(fn) = std(dimg(:) - img(:));
    fn = fn+1;
    f = 0;
  end
  f = f+1;
end
if ~isempty(anime)
  close(dovid);
  close(zovid);
end
if ~isempty(ofile)
  save([ofile,'.',datetimefilename,'.mat']);
  im2png(dimg,pngcmap,[ofile,'.',datetimefilename,'.denoise.png']);
  im2png(log10(z),pngcmap,[ofile,'.',datetimefilename,'.simap.png']);
  figure('Name','Method noise curve'),plot((1:length(beta))*rpf,beta)
  xlabel('repeats'),ylabel('\sigma'),title('Method noise \sigma')
  print('-depsc2',[ofile,'.',datetimefilename,'.beta.eps']);
end
if nargout == 0
  figure('Name','Input image');
  imagesc(img);colormap('jet');axis image;drawnow
  figure('Name','Denoised image');
  imagesc(dimg);colormap('jet');axis image;drawnow
else
  varargout{1} = dimg;
  varargout{2} = simg;
  varargout{3} = z;
  varargout{4} = beta;
end
return

function im2png(im,pngcmap,pngname)
  im = im-min(im(:));
  im = im/max(im(:));
  [lencmap,~] = size(pngcmap);
  if lencmap > 256
    disp('colormap has more than 256 entries is not supported.')
    imwrite(im,pngname,'png');
  else
    imwrite(im*lencmap,pngcmap,pngname,'png');
  end
return

function [img,rimg,pradius,alpha,h2,reps,wt,rot,ofile,anime,rpf,sizeI,...
  mc,fps,cmap,pngcmap] =...
  parse_inputs(varargin)
narginchk(1,100);
img = varargin{1};
% if iscell(img)
%     dict = img{2};
%     img = img{1};
% end
if ischar(img)
  disp(['Read input image from file: ',img])
  img = mean(double(imread(img)),3);
end
sizeI = size(img);
disp(['size of image: ',int2str(sizeI(1)),' x ',int2str(sizeI(2))])
pradius = [];pradius_d = 3;
alpha = [];alpha_d = 2;
wt = [];
h = [];
k = 2;
reps = []; reps_d = 100;
rot = []; rot_d = 0;
mc = []; mc_d = true;
rimg = []; rimg_d = img;
anime = '';
ofile = '';
rpf = []; rpf_d = 1;
fps = []; fps_d = 5;
cmap = []; cmap_d = jet(256);
pngcmap = []; pngcmap_d = gray(256);
while k < nargin
  switch lower(varargin{k})
    case {'wt','weight'}
      k = k+1;
      wt = varargin{k};
    case {'rimg','reference_image'}
      k = k+1;
      rimg = varargin{k};
    case {'r','radius'}
      k = k+1;
      pradius = varargin{k};
    case {'fps'}
      k = k+1;
      fps = varargin{k};
    case {'h'}
      k = k+1;
      h = varargin{k};
    case {'a','alpha'}
      k = k+1;
      alpha = varargin{k};
    case {'m','reps','repeats'}
      k = k+1;
      reps = varargin{k};
    case {'t','rot','rotate'}
      k = k+1;
      rot = varargin{k};
    case {'v','movie','anime'}
      k = k+1;
      anime = varargin{k};
    case {'c','mc','montecarlo'}
      k = k+1;
      mc = logical(varargin{k});
    case {'o','output'}
      k = k+1;
      ofile = varargin{k};
    case {'f','rpf','repsperframe'}
      k = k+1;
      rpf = varargin{k};
    case {'cmap','colormap'}
      k = k+1;
      cmap = varargin{k};
    case {'pngcmap','pngcolormap'}
      k = k+1;
      pngcmap = varargin{k};
    otherwise
      error('Unsupported parameters.')
  end
  k = k+1;
end
if isempty(pradius) && isempty(wt), pradius = pradius_d;end
if isempty(alpha) && isempty(wt), alpha = alpha_d;end
if isempty(pradius) && isempty(alpha) && ~isempty(wt)
  disp('custom weight function.')
  pradius = ceil(max(kernelSize(wt))/2);
  wt = pxnbh(wt,round(size(wt)/2),pradius);
end
if isempty(fps),     fps = fps_d;end
if isempty(cmap),    cmap = cmap_d;end
if isempty(pngcmap), pngcmap = pngcmap_d;end
if ischar(cmap),     cmap = stdcmap(cmap);end
if ischar(pngcmap),  pngcmap = stdcmap(pngcmap);end
if isempty(rimg),    rimg = rimg_d;end
if isempty(h)
  nlvl_img = sigmaclipping(img - sdwt2(img,1),3);
  nlvl_rimg = sigmaclipping(rimg - sdwt2(rimg,1),3);
  disp(['Estimated noise level: ', num2str(nlvl_img)]);
  h = sqrt(nlvl_img^2+nlvl_rimg^2);
end
if isempty(reps),  reps = reps_d;end
if isempty(rot),   rot = rot_d;end
if isempty(mc),    mc = mc_d;end
if isempty(rpf),   rpf = rpf_d;end
psize = [2*pradius+1 2*pradius+1];
if isempty(wt),wt = fspecial('gauss',psize,alpha);end
disp(['size of patch: ',int2str(psize(1)),' x ',int2str(psize(2))])
disp(['standard deviation of Gaussian weighted Euclidean distance: ',num2str(alpha)])
disp(['standard deviation of Gaussian weighted similarity: ',num2str(h)])
if ~isempty(anime)
  disp(['Animation filename: ', anime])
end
if ~isempty(ofile)
  disp(['Output filename prefix: ', ofile])
end
if rot > 0
  disp('Patch will be rotated.')
else
  disp('Patch will not be rotated.')
end
if mc == true
  disp('Image will be denoised by monte-carlo method.')
else
  disp('Image will be denoised from the first pixel to the last.')
end
h2 = h^2;
return

function cmap = stdcmap(cmapname)
cmapname = lower(cmapname);
switch cmapname
  case {'gray','jet','hsv','hot','cool','spring','autumn','summer','winter',...
      'bone','copper','pink','lines'}
    cmapfun = str2func(cmapname);
    cmap = cmapfun(256);
  otherwise
    error('Unsupported colormap.')
end
return
