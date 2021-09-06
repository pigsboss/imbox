function [c,w]=atrous2(x,J)
%ATROUS A Trous wavelet transform.
% The convolution is calculated through FFT.
% The wavelet function is a spline of degree 3.
%
% x - original function, must be N x N
% J - maximum scale
%
% Returns:
% c - maximum scaling function (redundant)
% w - wavelet from scale 1 to J
%
% Reference:
% J.-L. Starck and F. Murtagh, Astronomical Image and Data Analysis, 2nd
% Edition, Springer, p31-p32
c=x;
[H, W] = size(x);
if J > 0
    w=zeros(H, W, J);
    for k = 0:J-1
        hc = zeros(H,1);
        hc(1)=3/8;
        hc(mod(2^k,H)+1)=1/4;
        hc(mod(2^(k+1),H)+1)=1/16;
        hc(mod(H-2^k,H)+1)=1/4;
        hc(mod(H-2^(k+1),H)+1)=1/16;
        hr = zeros(1,W);
        hr(1)=3/8;
        hr(mod(2^k,W)+1)=1/4;
        hr(mod(2^(k+1),W)+1)=1/16;
        hr(mod(W-2^k,W)+1)=1/4;
        hr(mod(W-2^(k+1),W)+1)=1/16;
        h=hc*hr;
        cnext=fftshift(imconv(h,c));
        w(:,:,k+1)=c-cnext;
        c=cnext;
    end
else
    w=[];
end
return
