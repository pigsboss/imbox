function varargout = ellipse(varargin)
%ELLIPSE detect blobs in given image and use elliptical areas to fit them.
%
%Syntax:
%   ellipse
%   ellipse(I)
%   ellipse(I,scale)
%   ellipse(I,scale,threshold)
%   ells = ellipse(...)
%   [xc,yc,xr,yr,phi] = ellipse(...)
%
%Input arguments:
%
%Returns
[I, xgv, ygv, scale, threshold] = parse_inputs(varargin{:});
rx = xgv(1);
ry = ygv(1);
dx = mean(diff(xgv));
dy = mean(diff(ygv));
[B, L] = blob(I,scale,threshold);
[yc, xc] = find(B);
sizeL = size(L);
[xL, yL] = meshgrid(1:sizeL(2), 1:sizeL(1));
NB = length(yc);
ellxc = zeros(1,NB);ellyc = ellxc;ellxr = ellxc;ellyr = ellxc;ellphi = ellxc;elleof = ellxc; % error of fitting
R = 2*round(sqrt(numel(L)/NB));
for k = 1:NB
    [xb, yb] = border((xc(k)),(yc(k)),xL,yL,-L,R);
    [ellxc(k),ellyc(k),ellxr(k),ellyr(k),ellphi(k),elleof(k)] = ellifit(xb,yb);
    ellxr(k) = min(ellxr(k), R);
    ellyr(k) = min(ellyr(k), R);
end
nargoutchk(0,6)
ell = [(ellxc-1)*dx+rx;...
  (ellyc-1)*dy+ry;...
  ellxr*abs(dx);...
  ellyr*abs(dy);...
  ellphi*sign(dx)*sign(dy);...
  elleof];
switch nargout
    case 0
        figure;imagesc(xgv,ygv,I);axis image;colormap gray;drawnow
        hold all
        for k = 1:NB
%             if  ell(6,k)<0.25
%             if ell(3,k)>=2^(scale-1) && ell(3,k)<=2^(scale+2) && ell(6,k)<0.25
                ellidraw(ell(:,k),'Color','r');
%             end
        end
    case 1
        varargout{1} = ell;
    otherwise
        varargout = cell(nargout,1);
        for n = 1:nargout
            varargout{n} = ell(n,:);
        end
end
return

function [xb, yb] = border(xc,yc,xL,yL,L,R)
%BORDER border of dark blob
%
%Input arguments
%   (xc, yc) is center of the current blob
%   xL and yL are the coordinate matrices of Laplacian L
%   L is the Laplacian
%   R is the maximum radius of the current blob
%
%criteria of border:
%   (1) sign of L flip, or
%   (2) monotonicity flip.

% number of vertices on the polygon approximation of the border
N = 32;
% interpolation
NI = 40;
% direction angles
phi = 2*pi*(0:(N-1))/N;
xR = R*cos(phi);
yR = R*sin(phi);
% each column of XR is a radius from the center to one of the vertices.
XR = ((0:(NI-1))/NI)'*xR + xc;
YR = ((0:(NI-1))/NI)'*yR + yc;
LR = interp2(xL,yL,L,XR,YR,'linear',-1);
xb = zeros(1,N);
yb = zeros(1,N);
for k = 1:N
    s = find(LR(:,k)<0,1);
    if isempty(s),  s = NI;end
    d = find(diff(LR(:,k))>0,1)+1;
    if isempty(d),      d = NI;end
    xb(k) = XR(min(s,d),k);
    yb(k) = YR(min(s,d),k);
end
return

function [I, xgv, ygv, scale,...
  threshold, canvas] = parse_inputs(varargin)
I = [];
scale = [];scale_d = 0;
threshold = [];threshold_d = 0;
canvas=[];
xgv = [];
ygv = [];
narginchk(0,3);
switch nargin
    case 1
        I = varargin{1};
    case 2
        I = varargin{1};
        scale = varargin{2};
    case 3
        I = varargin{1};
        scale = varargin{2};
        threshold = varargin{3};
end
if isempty(I),      [I, xgv, ygv, canvas] = getcanvas;end
if isempty(I),      error('No image data specified.');end
sizeI = size(I);
if isempty(xgv),    xgv = 1:sizeI(2);end
if isempty(ygv),    ygv = 1:sizeI(1);end
if isempty(scale),  scale = scale_d;end
if isempty(threshold), threshold = threshold_d;end
if ischar(scale),   scale = str2double(scale);end
scale = round(max(scale,0));
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
