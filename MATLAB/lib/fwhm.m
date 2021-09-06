function [W,WX,WY] = fwhm(X)
%FWHM Full width at half maximum
[x,y] = find(X>0.5*max(X(:)));
WX = max(x)-min(x);
WY = max(y)-min(y);
W = max(WX,WY);
return
