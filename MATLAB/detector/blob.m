function varargout = blob(varargin)
%BLOB blob detector
%
%Syntax:
%   blob
%   blob(I)
%   blob(I,scale)
%   blob(I,scale,threshold)
%   blob(I,scale,threshold,scatterparam1,scatterparam2,...)
%   B = blob(...)
%   [B, L] = blob(...)
%
%Input arguments:
%   I     (optional) is input 2D image. If absence, current image CData
%   while be taken as input.
%   scale (optional) is scale of blob to detect.
%   threshold (optional) reject blobs that are less bright than threshold.
%
%Return:
%   B (optional) is a 2D logical array with the same size of input image,
%   which labels all detected blobs in the input image.
%   L (optional) is the Laplacian of the input image.
%
%Reference:
%   [1] http://en.wikipedia.org/wiki/Blob_detection
%   [2] http://en.wikipedia.org/wiki/Finite_difference_coefficients
%   [3] http://en.wikipedia.org/wiki/Gaussian_filter
[I, xgv, ygv, scale, ~, threshold, scatterparams] =...
  parse_inputs(varargin{:});
% L = dlt2(sdwt2(I,scale));
L = dlt2(gsmooth(I,scale))/max(max(fspecial('gauss',size(I),scale)));
Lp = 0.9*dlt2(gsmooth(I,scale*1.1))/...
  max(max(fspecial('gauss',size(I),scale*1.1)));
Ln = 0.9*dlt2(gsmooth(I,scale*0.9))/...
  max(max(fspecial('gauss',size(I),scale*0.9)));
% B = maxima(-L);
B = maxima3(-L,-Lp,-Ln);
if ~isempty(threshold)
  B = B & (I>=threshold);
end
if isempty(scatterparams)
  scatterparams={'r'};
end
nargoutchk(0,2)
switch nargout
    case 0
        figure;imagesc(xgv, ygv, I);
        hold all
        [yB,xB] = find(B);
        dx = abs(mean(diff(xgv)));
        dy = abs(mean(diff(ygv)));
        S = scale*pi/dx/dy;
        scattervars = [{xgv(xB)},{ygv(yB)},S,scatterparams{:}];
        scatter(scattervars{:});
        drawnow
    case 1
        varargout{1} = logical(B);
    case 2
        varargout{1} = logical(B);
        varargout{2} = L;
end
return

function [I, xgv, ygv, scale,...
  canvas, threshold, scatterparams] = parse_inputs(varargin)
I = [];
scale = [];scale_d = 0;
canvas=[];
threshold = [];
scatterparams = {};
narginchk(0,100);
xgv = [];
ygv = [];
switch nargin
  case 0
  case 1
    I = varargin{1};
  case 2
    I = varargin{1};
    scale = varargin{2};
  case 3
    I = varargin{1};
    scale = varargin{2};
    threshold = varargin{3};
  otherwise
    I = varargin{1};
    scale = varargin{2};
    threshold = varargin{3};
    scatterparams = {varargin(4:nargin)};
end
if isempty(I),      [I, xgv, ygv, canvas] = getcanvas;end
if isempty(I),      error('No image data specified.');end
sizeI = size(I);
if isempty(xgv),    xgv = 1:sizeI(2);end
if isempty(ygv),    ygv = 1:sizeI(1);end
if isempty(scale),  scale = scale_d;end
if ischar(scale),   scale = str2double(scale);end
scale = (max(scale,0));
return

function M = maxima(L)
% sigma = std(L(:));
% k = 0;
% t = k*sigma;
sizeL = size(L);
nElem = prod(sizeL);
nOps = 0; % number of fft operations
T = 10; % threshold
for k = 1:ndims(L)
    nffts = nElem/sizeL(k);
    nOps  = nOps + sizeL(k)*log2(sizeL(k))*nffts;
end
M = (L >= symmshift(L,[-1 -1])).*...
    (L >= symmshift(L,[-1  0])).*...
    (L >= symmshift(L,[-1  1])).*...
    (L >= symmshift(L,[ 0 -1])).*...
    (L >= symmshift(L,[ 0  1])).*...
    (L >= symmshift(L,[ 1 -1])).*...
    (L >= symmshift(L,[ 1  0])).*...
    (L >= symmshift(L,[ 1  1]));
M = M.*(L>eps*nOps*T);
return

function M = maxima3(L,Lp,Ln)
% sigma = std(L(:));
% k = 0;
% t = k*sigma;
sizeL = size(L);
nElem = prod(sizeL);
nOps = 0; % number of fft operations
T = 10; % threshold
for k = 1:ndims(L)
    nffts = nElem/sizeL(k);
    nOps  = nOps + sizeL(k)*log2(sizeL(k))*nffts;
end
M = (L >= symmshift(L,[-1 -1])).*...
    (L >= symmshift(L,[-1  0])).*...
    (L >= symmshift(L,[-1  1])).*...
    (L >= symmshift(L,[ 0 -1])).*...
    (L >= symmshift(L,[ 0  1])).*...
    (L >= symmshift(L,[ 1 -1])).*...
    (L >= symmshift(L,[ 1  0])).*...
    (L >= symmshift(L,[ 1  1])).*...
    (L >= symmshift(Lp,[-1 -1])).*...
    (L >= symmshift(Lp,[-1  0])).*...
    (L >= symmshift(Lp,[-1  1])).*...
    (L >= symmshift(Lp,[ 0 -1])).*...
    (L >= symmshift(Lp,[ 0  0])).*...
    (L >= symmshift(Lp,[ 0  1])).*...
    (L >= symmshift(Lp,[ 1 -1])).*...
    (L >= symmshift(Lp,[ 1  0])).*...
    (L >= symmshift(Lp,[ 1  1])).*...
    (L >= symmshift(Ln,[-1 -1])).*...
    (L >= symmshift(Ln,[-1  0])).*...
    (L >= symmshift(Ln,[-1  1])).*...
    (L >= symmshift(Ln,[ 0 -1])).*...
    (L >= symmshift(Ln,[ 0  0])).*...
    (L >= symmshift(Ln,[ 0  1])).*...
    (L >= symmshift(Ln,[ 1 -1])).*...
    (L >= symmshift(Ln,[ 1  0])).*...
    (L >= symmshift(Ln,[ 1  1]));
M = M.*(L>eps*nOps*T);
return

function [I, xgv, ygv, canvas] = getcanvas
%GETCANVAS get both CData (the image) and the handle of the current
%canvas object.
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
