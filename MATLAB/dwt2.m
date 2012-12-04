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
    J=floor(log(min(d1,d2))/log(2))-1;
end

if J==0
    c=x;
    w=zeros(size(c));
    return
end

% configuration
szEdge=2^(J); % size of purfled mirror, in pixels.

% mirroring boundary
%
% purfling with mirror edge
xpurf=mirrorpurfle(x,szEdge);
% find ceiling 2^n
D=max(d1+2*szEdge,d2+2*szEdge);
N=2^ceil(log(D)/log(2));
% zero-padding
X=zeros(N);
X(1:(d1+2*szEdge),1:(d2+2*szEdge))=xpurf;

% a trous wavelet transform
c=X;
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
        cnext=real(ifft2(fft2(c).*fft2(h)));
        w(:,:,k+1)=c-cnext;
        c=cnext;
    end
end
c=c((1+szEdge):(d1+szEdge),(1+szEdge):(d2+szEdge));
w=w((1+szEdge):(d1+szEdge),(1+szEdge):(d2+szEdge),:);
return
