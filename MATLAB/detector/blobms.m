function varargout = blobms(varargin)
%BLOBMS Multi-scale blob detector
%
%Syntax:
%   blobms
%   blobms(I)
%   blobms(I,scale)
%   blobms(I,scale,threshold)
%   blobms(I,scale,threshold,lineOpts)
%   [B, L, I, xgv, ygv, scale, lineOpts] = blobms(...)
%
%Input arguments:
%   I     (optional) is input 2D image. If absence, current image CData
%   while be taken as input.
%   scale (optional) is scale(s) of blob to detect.
%   threshold (optional) reject blobs that are less bright than threshold.
%
%Return:
%   B (optional) is the sizes of all detected blobs.
%   L (optional) is the Multi-scale Laplacian of the input image.
%
%Reference:
%   [1] http://en.wikipedia.org/wiki/Blob_detection
%   [2] http://en.wikipedia.org/wiki/Finite_difference_coefficients
%   [3] http://en.wikipedia.org/wiki/Gaussian_filter
MSGSIGMA_A = -0.94613;
MSGSIGMA_B = -1.0215;
IMTHRES = 0;
[I, xgv, ygv, scale, ~, sigma, threshold, lineopts] =...
  parse_inputs(varargin{:});
if isnumeric(threshold)
  if isscalar(threshold)
    IMTHRES = threshold;
    MSGTHRES = 0*scale;
  elseif numel(threshold) == numel(scale)
    MSGTHRES = threshold;
  else
    error('If threshold is a vector, its size must be the same as scale vector.')
  end
elseif ischar(threshold)
  switch lower(threshold)
    case {'auto'}
      if ~isempty(sigma)
        MSGTHRES = 3 * sigma * exp(MSGSIGMA_B)*scale.^MSGSIGMA_A;
      else
        MSGTHRES = zeros(size(scale));
      end
    case {'off'}
      MSGTHRES = zeros(size(scale));
    otherwise
      error('Unsupported threshold.')
  end
else
  error('Unsupported threshold.')
end
FFTTHRES = fftthreshold(size(I));
L = -1*laplams(I,scale);
M = maxima3(L);
B = zeros(size(I));
Ovl = zeros(size(I));
for k=2:(length(scale)-1)
  M(:,:,k) = M(:,:,k) & (L(:,:,k)>max(FFTTHRES,MSGTHRES(k))) & (I>IMTHRES);
  B = B + double(M(:,:,k))*scale(k);
  Ovl = Ovl + double(M(:,:,k));
end
[y,x] = find(Ovl>=2);
for k = 1:length(x)
  [~,idx] = max(squeeze(L(y(k),x(k),:)));
  B(y(k),x(k)) = scale(idx);
end
if isempty(lineopts)
  lineopts={'color','r'};
end

% nargoutchk(0,5)
switch nargout
  case 0
    figure;imagesc(xgv, ygv, I);colormap('gray');axis image
    [X,Y] = meshgrid(xgv,ygv);
    hold all
    idx = find(B);
    dx = abs(mean(diff(xgv)));
    dy = abs(mean(diff(ygv)));
    S = 2*pi*B(idx);
    ellidrawvars = [{[X(idx),Y(idx),S(:)*dx,S(:)*dy,zeros(size(S(:)))]'}, lineopts{:}];
    ellidraw(ellidrawvars{:});
    drawnow
  otherwise
    varargout{1} = B;
    varargout{2} = L;
    varargout{3} = I;
    varargout{4} = xgv;
    varargout{5} = ygv;
    varargout{6} = scale;
    varargout{7} = varargin;
end
return

function [I, xgv, ygv, scale,...
  canvas, sigma, threshold, lineopts, varargin] = parse_inputs(varargin)
I = [];
scale = [];
canvas=[];
narginchk(0,100);
xgv = [];
ygv = [];
NS = 20;
sigma = [];
threshold = [];
threshold_d = 'auto';
lineopts = [];
switch nargin
  case 0
  case 1
    I = varargin{1};
  case 2
    I = varargin{1};
    scale = varargin{2};
  otherwise
    I = varargin{1};
    scale = varargin{2};
    lineopts = getoptval(varargin,{'lineopts'});
    threshold = getoptval(varargin,{'threshold'});
    sigma = getoptval(varargin,{'sigma'});
end
if isempty(threshold), threshold = threshold_d;end
if isempty(I),      [I, xgv, ygv, canvas] = getcanvas;end
if isempty(I),      error('No image data specified.');end
sizeI = size(I);
if isempty(xgv),    xgv = 1:sizeI(2);end
if isempty(ygv),    ygv = 1:sizeI(1);end
if isempty(scale)
    scmax = min(sizeI)/NS;
    scale = ((1:NS)/NS*scmax);
end
if ischar(scale),   scale = str2double(scale);end
if length(scale)==2
    scmin = min(scale);
    scmax = max(scale);
    scale = (scmin:((scmax-scmin)/NS):scmax);
end
if length(scale)==1
    scmax = scale;
    scale = ((1:NS)/NS*scmax);
end
if any(scale<=0)
    error('scale must be positive.')
end
scale = (max(scale,0));
return

function [I, xgv, ygv, canvas] = getcanvas
%GETCANVAS get both CData (the image) and the handle of the current canvas object.
I = [];
canvas = [];
hf=get(0,'CurrentFigure');
if ~isempty(hf)
	ha=get(hf,'CurrentAxes');
	if ~isempty(ha)
		canvas=findobj(ha,'Type','image');
		if ~isempty(canvas)
			I=get(canvas,'CData');
      xgv = get(canvas,'XData');
      ygv = get(canvas,'YData');
      if numel(xgv)*numel(ygv) ~= numel(I)
        xgv = xgv(1):xgv(2);
        ygv = ygv(1):ygv(2);
      end
		end
	end
end
return
