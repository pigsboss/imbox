function [c,w]=dwt2(x,J)
%DWT2 2-D discrete wavelet transform
% A Trous method is used. The convolution is calculated through FFT.
%
% x - original function
% J - maximum scale
%
% Returns:
% c - maximum scaling function (redundant)
% w - wavelet from scale 1 to J
%
[d1,d2]=size(x);

if nargin==1
    J=lastpow2(min(d1,d2));
end
if J==0
    c=x;
    w=zeros(size(c));
    return
end
% kernel padding
kPadSz=2^J;
x=padarray(x,[kPadSz kPadSz],'symmetric','both');
% square padding
N=2^nextpow2(max(d1,d2)+2*kPadSz);
sPadSz=N*ones(1,2)-[d1 d2];
x=padarray(x,sPadSz,'symmetric','post');
% a trous wavelet transform
[c,w]=atrous2(x,J);
c=c((1+kPadSz):(d1+kPadSz),(1+kPadSz):(d2+kPadSz));
w=w((1+kPadSz):(d1+kPadSz),(1+kPadSz):(d2+kPadSz),:);
return

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
[N,~]=size(x);
if J > 0
    w=zeros(N,N,J);
    for k=0:J-1
        h=zeros(N,1);
        h(1)=3/8;
        h(mod(2^k,N)+1)=1/4;
        h(mod(2^(k+1),N)+1)=1/16;
        h(mod(N-2^k,N)+1)=1/4;
        h(mod(N-2^(k+1),N)+1)=1/16;
        h=h*h';
        cnext=fftshift(imconv(h,c));
        w(:,:,k+1)=c-cnext;
        c=cnext;
    end
else
    w=[];
end
return