function varargout = linefit(x,y,wgt)
%LINEFIT Linear fitting
%SYNTAX
% [a b r] = linefit(x,y,wgt)
%
%INPUT
% x and y are data points.
% wgt     is weight of each point.
%RETURN
% a and b are linear parameters s.t. y = a*x + b.
% r       is the linear correlation coefficient.

if nargin<3
  wgt = ones(size(x));
end
xw = x(:).*wgt(:);
yw = y(:).*wgt(:);
xmean = sum(xw(:));
ymean = sum(yw(:));
wmean = sum(wgt(:));
x2mean = sum(xw(:).^2);
y2mean = sum(yw(:).^2);
xymean = sum(xw(:).*yw(:));
A = [x2mean xmean; xmean wmean];
B = [xymean; ymean];
c = linsolve(A,B);
dx2 = x2mean - xmean^2/wmean;
dy2 = y2mean - ymean^2/wmean;
dxy = xymean - xmean*ymean/wmean;
r = dxy/sqrt(dx2*dy2);
yerr = y(:) - x(:)*c(1) - c(2);
rms = std(yerr(:));
switch nargout
  case 0
    figure('Name','Data and line fit')
    scatter(x(:),y(:),'rx')
    hold all
    plot(x(:),c(1)*x(:)+c(2),'b')
  case 1
    varargout{1} = [c;r;rms];
  otherwise
    varargout{1} = c(1);
    varargout{2} = c(2);
    varargout{3} = r;
    varargout{4} = rms;
end
return
