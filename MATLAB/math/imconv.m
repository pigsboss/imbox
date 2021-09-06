function g=imconv(varargin)
%IMCONV image convolution
%
%Syntax:
%   g = imconv(f,h)
%   g = imconv(f,h,boundary)
%   g = imconv(f,h,boundary,padsize)
%
%Input arguments:
%   f          is input 2D image
%   h          is input 2D discrete filter (PSF)
%   boundary   (optional) can be either 'circular' or 'symmetric', or
%   an arbitrary scalar.
%   padsize    (optional) can be either scaler or vector
f = varargin{1};
typef = str2func(class(f));
h = varargin{2};
[H,W,L] = size(f);
[~,~,lenPSF] = size(h);
g = zeros(H,W,L,lenPSF);
for m = 1:L
  for n = 1:lenPSF
    g(:,:,m,n) = coreimconv(double(f(:,:,m)),double(h(:,:,n)),...
      varargin{3:nargin});
  end
end
if L==1
  g = reshape(g,[H,W,lenPSF]);
end
if lenPSF==1
  g = reshape(g,[H,W,L]);
end
% if lenPSF > 1
%   g = h;
%   for k = 1:lenPSF
%     g(:,:,k) = coreimconv(f,h(:,:,k),varargin{3:nargin});
%   end
% else
%   g = coreimconv(varargin{:});
% end
g = typef(g);
return

function g=coreimconv(varargin)
narginchk(2,4)
f = varargin{1};
h = varargin{2};
sizeF = size(f);
sizeH = size(h);
if sizeH(1) > sizeF(1) || sizeH(2) > sizeF(2)
    error('filter must not be larger than input data.')
end
if sizeH(1) < sizeF(1) || sizeH(2) < sizeF(2)
    h = ifftshift(h);
    h = [h(1:round(sizeH(1)/2),:);...
        zeros(sizeF(1)-sizeH(1),sizeH(2));...
        h((1+round(sizeH(1)/2)):sizeH(1),:)];
    sizeH = size(h);
    h = [h(:,1:round(sizeH(2)/2)),...
        zeros(sizeH(1),sizeF(2)-sizeH(2)),...
        h(:,(1+round(sizeH(2)/2)):sizeH(2))];
%     sizeH = size(h);
    h = fftshift(h);
end
boundary = [];boundary_d = 'circular';
padsize = [];
switch nargin
    case 3
        boundary = varargin{3};
    case 4
        boundary = varargin{3};
        padsize = varargin{4};
end
if isempty(boundary),   boundary = boundary_d;end
if strcmp(boundary,'circular')
    H = psf2otf(h);
    g = real(ifft2(H.*fft2((f))));
    return
end
if isempty(padsize),    padsize = kernelSize(h);end
if isscalar(padsize),   padsize = padsize*ones(1,2);end
f = padarray(f,padsize,boundary,'both');
h = padarray(h,padsize,0,'both');
H = psf2otf(h);
g = real(ifft2(H.*fft2((f))));
g = g((1+padsize(1)):(sizeF(1)+padsize(1)), (1+padsize(2)):(sizeF(2)+padsize(2)));
return

