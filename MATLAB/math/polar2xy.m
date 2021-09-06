function [fxy,x,y,rep] = polar2xy(varargin)
%POLAR2XY Transform f(rho,phi) on polar coordinate grid to f(x,y) on
%cartesian coordinate grid.
%SYNTAX:
%  [fxy,x,y] = polar2xy(fp,rho,phi,x,y,method,extrapval)
%  [fxy,x,y] = polar2xy(fp,rho,phi,x,y,extrapval)
%  [fxy,x,y] = polar2xy(fp,rho,phi,x,y)
%  [fxy,x,y] = polar2xy(fp,rho,phi,Nxy,extrapval)
%  [fxy,x,y] = polar2xy(fp,rho,phi,Nxy)
%  [fxy,x,y] = polar2xy(fp,rho,phi)
%  [fxy,x,y] = polar2xy(fp,Nxy,method)
%  [fxy,x,y] = polar2xy(fp,method)
%  [fxy,x,y] = polar2xy(fp,Nxy)
%  [fxy,x,y] = polar2xy(fp)
%INPUTS:
%  fp          is the 2-D function defined on polar grid.
%  rho and phi are radial and angular coordinates.
%  x and y     are x and y coordinates.
%  Nxy         is sample frequency along each axis of cartesian coordinate
%              grid.
%  extrapval   is extrapolate value for undefined data on the cartesian
%              coordinate grid.

[fp,~,~,rho,phi,Nxy,x,y,method,extrapval] = parse_inputs(varargin{:});
fxy = zeros(Nxy);
rep = zeros(Nxy);
px = rho.*cos(phi);
py = rho.*sin(phi);
if ismonotonic(x) && ismonotonic(y) && ~strcmpi(method,'tri')
  xgv = x(1,:);
  ygv = y(:,1);
  xref = xgv(1);
  yref = ygv(1);
  ux = mean(diff(xgv));
  uy = mean(diff(ygv));
  xi = round(((px(:)-xref)/ux))+1;
  yi = round(((py(:)-yref)/uy))+1;
  p2xy = ((xi-1)*Nxy + yi).*double(xi>=1).*double(xi<=(Nxy)).*...
    double(yi>=1).*double(yi<=(Nxy));
  for k = 1:numel(fp)
    if p2xy(k)>0
      fxy(p2xy(k)) = fxy(p2xy(k))+fp(k);
      rep(p2xy(k)) = rep(p2xy(k))+1;
    end
  end
  fxy(rep>0) = fxy(rep>0)./rep(rep>0);
  fxy(rep==0) = extrapval;
else
  [XI,YI] = meshgrid(1:Nxy);
  p2xi = min(max(round(xy2polar(XI,x,y,rho,phi,'tri',[])),1),Nxy);
  p2yi = min(max(round(xy2polar(YI,x,y,rho,phi,'tri',[])),1),Nxy);
  for k = 1:numel(fp)
    fxy(p2yi(k),p2xi(k)) = fxy(p2yi(k),p2xi(k))+fp(k);
    rep(p2yi(k),p2xi(k)) = rep(p2yi(k),p2xi(k))+1;
  end
  fxy(rep>0) = fxy(rep>0)./rep(rep>0);
  fxy(rep==0) = extrapval;
end
return

function [fp,drho,dphi,rho,phi,Nxy,x,y,method,extrapval] = parse_inputs(varargin)
narginchk(1,7)
fp = varargin{1};
[dphi,drho] = size(fp);
rho = [];
phi = [];
Nxy = [];
x = [];
y = [];
extrapval = []; extrapval_d = 0;
method = ''; method_d = 'average';
resample = 4;
switch nargin
  case 2
    if isscalar(varargin{2})
      Nxy = varargin{2};
    elseif ischar(varargin{2})
      method = varargin{2};
    end
  case 3
    if isscalar(varargin{2}) && ischar(varargin{3})
      Nxy = varargin{2};
      method = varargin{3};
    else
      rho = varargin{2};
      phi = varargin{3};
    end
  case 4
    rho = varargin{2};
    phi = varargin{3};
    Nxy = varargin{4};
  case 5
    rho = varargin{2};
    phi = varargin{3};
    if isinteger(varargin{4}) && isscalar(varargin{5})
      Nxy = varargin{4};
      extrapval = varargin{5};
    else
      x = varargin{4};
      y = varargin{5};
    end
  case 6
    rho = varargin{2};
    phi = varargin{3};
    x = varargin{4};
    y = varargin{5};
    extrapval = varargin{6};
  case 7
    rho = varargin{2};
    phi = varargin{3};
    x = varargin{4};
    y = varargin{5};
    method = varargin{6};
    extrapval = varargin{7};
end
if isempty(extrapval), extrapval = extrapval_d;end
if isempty(method),    method = method_d;end
if isempty(Nxy),         Nxy = max(drho*2,dphi)/resample;end
if isempty(rho) && isempty(x)
  xgv = (1:Nxy) - (1+Nxy)/2;
  ygv = xgv;
  [x,y] = meshgrid(xgv,ygv);
  crho = sqrt(x.^2+y.^2);
  cphi = atan2(y,x);
  rmax = max(crho(:));
  rmin = min(crho(:));
  pmax = max(cphi(:));
  pmin = min(cphi(:));
  rgv = ((1:drho)-1)/(drho-1)*(rmax-rmin)+rmin;
  pgv = ((1:dphi)-1)/(dphi-1)*(pmax-pmin)+pmin;
  [rho,phi] = meshgrid(rgv,pgv);
end
px = rho.*cos(phi);
py = rho.*sin(phi);
if isempty(x)
  xmin = min(px(:));
  xmax = max(px(:));
  ymin = min(py(:));
  ymax = max(py(:));
  xgv = ((1:Nxy)-1) / (Nxy-1) * (xmax-xmin) + xmin;
  ygv = ((1:Nxy)-1) / (Nxy-1) * (ymax-ymin) + ymin;
  [x,y] = meshgrid(xgv,ygv);
end
return