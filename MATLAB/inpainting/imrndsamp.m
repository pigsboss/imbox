function [samps,m] = imrndsamp(img,pct,method)
%IMRNDSAMP Random samplings of image.
%INPUTS
% img    is the input image.
% pct    is the sampled percentage of all pixels.
% method is the sampling method.
%
%SUPPORTED METHOD
% 'uniform'
% 'tv'
%
%RETURNS
% samps is sampled portion of the input image.
% m     is the mask of the sampling.
method_d = 'uniform';
beta = 0.5;
if nargin==2
  method = method_d;
end
typeI = str2func(class(img));
sizeI = size(img);
ndims = length(sizeI);
H = sizeI(1);
W = sizeI(2);
pct = 1 - min(max(pct,0),100)/100;
if ndims == 3
  D = sizeI(3);
else
  D = 1;
end
img = double(img);
switch method
  case {'uniform'}
    x = rand(H,W);
    y = zeros(H,W,D);
    for k = 1:D
      y(:,:,k) = x;
    end
    m = (y>=pct);
    samps = typeI(img.*double(m));
  case {'tv'}
    x = rand(H,W);
    y = zeros(H,W,D);
    v = zeros(H,W);
    for k = 1:D
      v = v+dlt2(sdwt2(img(:,:,k),1)).^2;
    end
    v = (v/max(v(:))).^(1/2);
    x = (1+beta*v).*x/(1+beta);
    for k = 1:D
      y(:,:,k) = x;
    end
    c = sort(x(:));
    idx = min(max(1,round(pct*numel(x))),numel(x));
    m = (y>=c(idx));
    samps = typeI(img.*double(m));
  otherwise
    error([method,' sampling is not supported.'])
end
return
