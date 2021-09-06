function [r, idx] = findroot(x,y,PRECISION)
%FINDROOT find all roots of equation y(x) = 0
narginchk(1,3)
PRECISION_d = 1e-3;
switch nargin
  case 1
    y = x;
    x = 1:length(y);
    PRECISION = PRECISION_d;
  case 2
    PRECISION = PRECISION_d;
end
x = x(:);
if std(y)>eps
  y = normalize(y);
else
  warning('imbox:findroot:yIsConstant',['y(x) is constant function.',...
    'Please try to increase sampling frequency.'])
end
y = round(y(:)/PRECISION)*PRECISION;
N = length(y);
rc1 = x(y==0); % root of class 1
rc2idx = find(((y(1:(N-1))>0) & (y(2:N)<0)) |...
  ((y(1:(N-1))<0) & (y(2:N)>0)));
rc2 = (x(rc2idx).*abs(y(rc2idx)) + x(rc2idx+1).*abs(y(rc2idx+1))) ./...
  abs(y(rc2idx) - y(rc2idx+1));
r = sort([rc1(:);rc2(:)]);
dx = mean(diff(x));
idx = (r-x(1))/dx + 1;
if isempty(r)
    disp('No root found with current precision.')
end
return