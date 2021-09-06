function I = cat2map(varargin)
%CAT2MAP Category to map
switch nargin
  case 2
    T = varargin{1};
    sizeI = varargin{2};
    xgv = 1:sizeI(2);
    ygv = 1:sizeI(1);
  case 3
    T = varargin{1};
    if isvector(varargin{2})
      xgv = varargin{2};
      ygv = varargin{3};
    else
      xgv = varargin{2}(1,:);
      ygv = varargin{3}(:,1);
    end
    sizeI = [length(ygv),length(xgv)];
end
[NUMROW,lenT] = size(T);
if NUMROW < 3
    error('Input table must contain 3 rows at least.')
end
f = T(1,:);
x = T(2,:);
y = T(3,:);
I = zeros(sizeI);
for n = 1:lenT
  if x(n) >= min(xgv) && x(n) <= max(xgv) &&...
      y(n) >= min(ygv) && y(n) <= max(ygv)
    [~,c] = min(abs(xgv - x(n)));
    [~,r] = min(abs(ygv - y(n)));
    I(r,c)=I(r,c)+f(n);
  end
end
return
