function [c,w]=dwt2(varargin)
%DWT2 2-D discrete wavelet transform
% A Trous method is used. The convolution is calculated through FFT.
%
% Syntax:
%   [c,w]=dwt2(x)
%   [c,w]=dwt2(x,J)
%   [c,w]=dwt2(x,J,sz)
%
% x - original function, NxN matrix
% J (optional, default: lastpow2(N)-2) - maximum scale
% sz (optional, default: 0) - padding size
%
% Returns:
% c - maximum scaling function (redundant)
% w - wavelet from scale 1 to J
%
%% parsing input arguments:
J=[];
kPadSz=[];
switch nargin
    case 1
        x=varargin{1};
    case 2
        x=varargin{1};
        J=varargin{2};
    case 3
        x=varargin{1};
        J=varargin{2};
        kPadSz=varargin{3};
    otherwise
        error('Syntax error')
end

%% setup initial parameters:
[d1,d2,ch]=size(x);
if isempty(J) % the max_scale is determined automatically
    J=lastpow2(min(d1,d2));
end
if J<=0
    c=x;
    w=zeros(size(c));
    return
end
if isempty(kPadSz) % the padding size is determined automatically
    kPadSz=2^(J+1);
end

% kernel padding
x=padarray(x,[kPadSz kPadSz],'symmetric','both');

% square padding
% N=2^nextpow2(max(d1,d2)+2*kPadSz)-2*kPadSz;
% sPadSz=N*ones(1,2)-[d1 d2];
% x=padarray(x,sPadSz,'symmetric','post');

% a trous wavelet transform
if ch>1
    sz=size(x);
    c=zeros(sz);
    w=zeros(sz(1),sz(2),J,sz(3));
    for k=1:ch
        [c(:,:,k),w(:,:,:,k)]=atrous2(x(:,:,k),J);
    end
    c=c((1+kPadSz):(d1+kPadSz),(1+kPadSz):(d2+kPadSz),:);
    w=w((1+kPadSz):(d1+kPadSz),(1+kPadSz):(d2+kPadSz),:,:);
else
    [c,w]=atrous2(x,J);
    c=c((1+kPadSz):(d1+kPadSz),(1+kPadSz):(d2+kPadSz));
    w=w((1+kPadSz):(d1+kPadSz),(1+kPadSz):(d2+kPadSz),:);
end
return
