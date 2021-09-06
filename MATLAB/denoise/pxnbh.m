function p = pxnbh(img,vc,radius,angle)
%PXNBH Neighborhood of given pixel and radius.
%
%INPUTS
% img    is the image. A neighborhood is a subset of the image.
% vc     is position vector of central pixel, in [row, col].
% radius is the radius of the square neighborhood.
%
%RETURN
% p      is the neighborhood as a (2*radius + 1)^2 patch.
xgv = (vc(2)-radius):(vc(2)+radius);
ygv = (vc(1)-radius):(vc(1)+radius);
if nargin == 3
  p = subim(img,ygv,xgv,'s');
else
  [x,y] = meshgrid(xgv-vc(2),ygv-vc(1));
  Rmat = [cos(angle),-sin(angle);sin(angle),cos(angle)];
  r = Rmat * [x(:)';y(:)'];
  X = reshape(r(1,:)+vc(2),size(x));
  Y = reshape(r(2,:)+vc(1),size(y));
  p = subim(img,Y,X,'p');
end
return
