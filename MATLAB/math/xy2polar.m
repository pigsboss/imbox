function [fp,rho,phi,x,y] = xy2polar(varargin)
%XY2POLAR Transform function f(x,y) on cartesian coordinate grid to
%f(rho,phi) on polar coordinate grid.
%SYNTAX:
%  [fp,rho,phi] = xy2polar(fxy,x,y,rho,phi,method,extrapval)
%  [fp,rho,phi] = xy2polar(fxy,x,y,method,extrapval,resample_rate)
%  [fp,rho,phi] = xy2polar(fxy,x,y,method,extrapval)
%  [fp,rho,phi] = xy2polar(fxy,x,y,method)
%  [fp,rho,phi] = xy2polar(fxy,x,y)
%  [fp,rho,phi] = xy2polar(fxy,method)
%  [fp,rho,phi] = xy2polar(fxy)
%INPUTS:
%  fxy     is the 2-D function defined on cartesian coordinate grid.
%  x and y are coordinates of cartesian grid.
%  method  is a string specifying the method of sampling. It can be one of
%          the following: 'tri', means triangle; 'nearest', 'linear',
%          'cubic', and 'spline'.
%  extrapval is value used for extrapolation.
%  resample_rate is resample frequency in proportion to the sample
%  frequency of fxy. >1 means up sample while <1 means down sample.
%          Because the polar grid is equidistant in radial and angular
%          coordinates, sampling on this grid is irregular. To make sure
%          the sampling is complete for the worst case, resample_rate
%          should greater than 4.

[fxy,x,y,rho,phi,method,extrapval] = parse_inputs(varargin{:});
method = lower(method);
X=rho.*cos(phi);
Y=rho.*sin(phi);
switch method
  case {'tri'}
    F=TriScatteredInterp(x(:),y(:),fxy(:));
    fp=F(X,Y);
    fp(isnan(fp)) = extrapval;
  case {'nearest'}
    [dy,dx] = size(fxy);
    ux = mean(diff(x(1,:)));
    uy = mean(diff(y(:,1)));
    xref = x(1,1);
    yref = y(1,1);
    xi = round((X(:)-xref)/ux) + 1;
    yi = round((Y(:)-yref)/uy) + 1;
    p2xy = ((xi-1)*dy+yi).*double(xi>=1).*double(xi<=dx).*...
      double(yi>=1).*double(yi<=dy);
    fp = extrapval*ones(size(X));
    for k = 1:numel(X)
      if p2xy(k)>0
        fp(k) = fxy(p2xy(k));
      end
    end
  otherwise
    fp = interp2(x,y,fxy,X,Y,method,extrapval);
end
return

function [fxy,x,y,rho,phi,method,extrapval] =...
  parse_inputs(varargin)
narginchk(1,7)
xgv = [];
ygv = [];
rho = [];
phi = [];
extrapval = []; extrapval_d = 0;
method = ''; method_d = 'nearest';
resample = []; resample_d = 4;
fxy = varargin{1};
[dy,dx]=size(fxy);
switch nargin
  case 2
    if ischar(varargin{2})
      method = varargin{2};
    else
      error('syntax error.')
    end
  case 3
    xgv = varargin{2};
    ygv = varargin{3};
  case 4
    xgv = varargin{2};
    ygv = varargin{3};
    method = varargin{4};
  case 5
    xgv = varargin{2};
    ygv = varargin{3};
    method = varargin{4};
    extrapval = varargin{5};
  case 6
    xgv = varargin{2};
    ygv = varargin{3};
    method = varargin{4};
    extrapval = varargin{5};
    resample = varargin{6};
  case 7
    xgv = varargin{2};
    ygv = varargin{3};
    rho = varargin{4};
    phi = varargin{5};
    method = varargin{6};
    extrapval = varargin{7};
end
if isempty(xgv),       xgv = (1:dx) - (1+dx)/2;end
if isempty(ygv),       ygv = (1:dy) - (1+dy)/2;end
if isempty(method),    method = method_d;end
if isempty(extrapval), extrapval = extrapval_d;end
if isempty(resample),  resample = resample_d;end
if isscalar(xgv),      xgv = (1:dx) - xgv;end
if isscalar(ygv),      ygv = (1:dy) - ygv;end
if ~strcmpi(method,'tri')
  if ~ismonotonic(xgv) || ~ismonotonic(ygv)
    error('x and y must be monotonic.');
  end
  if ~isvector(xgv)
    x = xgv;
    y = ygv;
  else
    [x,y] = meshgrid(xgv,ygv);
  end
else
  if isvector(xgv)
    [x,y] = meshgrid(xgv,ygv);
  else
    x = xgv;
    y = ygv;
  end
end
if isempty(rho)
  rho=sqrt(x(:).^2+y(:).^2);
  phi=atan2(y(:),x(:));
  rmax=max(rho);
  rmin=min(rho);
  pmax=max(phi);
  pmin=min(phi);
  Nxy=max(dx,dy);
  Nrho = Nxy / 2 * resample;
  Nphi = Nxy * resample;
  rgv=(((1:Nrho)-1)/(Nrho-1))*(rmax-rmin)+rmin;
  if any(x(:)>=0)
    pgv=(((1:Nphi)-1)/(Nphi-1))*(pmax-pmin)+pmin;
  else
    pmin = mod(pmin,2*pi);
    pmax = mod(pmax,2*pi);
    pgv=(((1:Nphi)-1)/(Nphi-1))*(pmax-pmin)+pmin;
  end
  [rho,phi]=meshgrid(rgv,pgv);
end
return
