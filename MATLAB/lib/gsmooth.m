function S = gsmooth(I,scale,boundary)
%GSMOOTH Smooth input image with gaussian filter.
%Syntax
%  S = gsmooth(I, scale)
%  S = gsmooth(I, scale, boundary)
%
%Description
%  I                   is the original data. I can be either 1D, 2D, or 3D.
%  scale               is the scale of Gaussian kernel, which is measured
%  with the sigma parameter of Gaussian kernel function.
%  boundary (optional) can be 'circular', 'symmetric', or an arbitrary
%  scalar. The default is 'symmetric'.
%
%See Also
%imconv, fspecial

if nargin == 2
  boundary = 'symmetric';
end
szI = size(I);
if scale > 0
  S = imconv(I, fspecial('gauss',szI(1:2),scale), boundary);
else
  S = I;
end
return