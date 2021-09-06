function y=zeroshift2(x,shift)
%ZEROSHIFT2 Works just like circshift on 2-D matrix, except circshift
% pads bottom (right) elements to the top (left), vice versa, while
% zeroshift2 pads 0 elements.
[d1,d2] = size(x);
y=zeros(d1,d2);
if shift(1)>=0
  if shift(2)>=0
    y((1+shift(1)):d1,(1+shift(2)):d2)=...
      x(1:(d1-shift(1)),1:(d2-shift(2)));
  else
    y((1+shift(1)):d1,1:(d2+shift(2)))=...
      x(1:(d1-shift(1)),(1-shift(2)):d2);
  end
else
  if shift(2)>=0
    y(1:(d1+shift(1)),(1+shift(2)):d2)=...
      x((1-shift(1)):d1,1:(d2-shift(2)));
  else
    y(1:(d1+shift(1)),1:(d2+shift(2)))=...
      x((1-shift(1)):d1,(1-shift(2)):d2);
  end
end
return
