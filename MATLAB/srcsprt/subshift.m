function S = subshift(I,v,varargin)
%SUBSHIFT Sub-pixel shift of 1D or 2D image.
%INPUT
% I             is the 1D or 2D image to be shifted.
% v             is the shifting vector.
% interp_method is a string identifier of the method used for interpolation.
%RETURN
% S is the shifted image.
narginchk(2,100);
if nargin < 3
  varargin = {'linear',0};
end
if isvector(I)
  d = size(I(:));
  XI = 1:d;
  XS = XI - v;
  if ~isempty(varargin)
    S = interp1(XI,I,XS,varargin{:});
  else
    S = interp1(XI,I,XS);
  end
else
  [dy,dx] = size(I);
  xgv = 1:dx;
  ygv = 1:dy;
  [XI,YI] = meshgrid(xgv,ygv);
  [XS,YS] = meshgrid(xgv-v(2),ygv-v(1));
  if ~isempty(varargin)
    S = interp2(XI,YI,I,XS,YS,varargin{:});
  else
    S = interp2(XI,YI,I,XS,YS);
  end
end
return