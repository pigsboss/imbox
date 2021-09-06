function y=zeroshift(x,shift)
%ZEROSHIFT Works just like circshift, except circshift
% pads bottom (right) elements to the top (left), vice versa, while
% zeroshift pads 0 elements.
sz = size(x);
ndims = length(sz);
d1 = sz(1);
switch ndims
  case 1
    d2 = 1;
    d3 = 1;
  case 2
    d2 = sz(2);
    d3 = 1;
  case 3
    d2 = sz(2);
    d3 = sz(3);
end
typex = str2func(class(x));
y=typex(zeros(d1,d2,d3));
sft = zeros(1,3);
for k = 1:numel(shift)
  sft(k) = shift(k);
end
if sft(1)>=0
  xgv1 = 1:(d1-sft(1));
  ygv1 = (1+sft(1)):d1;
else
  xgv1 = (1-sft(1)):d1;
  ygv1 = 1:(d1+sft(1));
end
if sft(2)>=0
  xgv2 = 1:(d2-sft(2));
  ygv2 = (1+sft(2)):d2;
else
  xgv2 = (1-sft(2)):d2;
  ygv2 = 1:(d2+sft(2));
end
if sft(3)>=0
  xgv3 = 1:(d3-sft(3));
  ygv3 = (1+sft(3)):d3;
else
  xgv3 = (1-sft(3)):d3;
  ygv3 = 1:(d3+sft(3));
end
y(ygv1,ygv2,ygv3) = x(xgv1,xgv2,xgv3);
return
