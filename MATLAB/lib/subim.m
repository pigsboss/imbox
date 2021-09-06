function [p,rgv,cgv] = subim(img,rgv,cgv,ext)
%SUBIM Return part of the given image
%
%INPUTS
% img is the image.
% rgv is the row grid vector.
% cgv is the column grid vector.
% ext is the method of extension.
%   'symmetric'
%   'periodic'
%   'clip'
%
%RETURN
% p is selected part of the image.

sizeI = size(img);
if nargin < 4
  ext = 's';
end
ext = lower(ext);
rgv = round(rgv);
cgv = round(cgv);
switch ext
  case {'s','sym','symmetric'}
    rgv = mod(rgv-1,sizeI(1)*2)+1;
    cgv = mod(cgv-1,sizeI(2)*2)+1;
    rrev = logical(rgv>sizeI(1));
    crev = logical(cgv>sizeI(2));
    rgv(rrev) = 2*sizeI(1)+1-rgv(rrev);
    cgv(crev) = 2*sizeI(2)+1-cgv(crev);
  case {'p','per','periodic'}
    rgv = mod(rgv-1,sizeI(1))+1;
    cgv = mod(cgv-1,sizeI(2))+1;
  case {'c','clip','clipping'}
    rgv = rgv((rgv<=sizeI(1)) & (rgv>=1));
    cgv = cgv((cgv<=sizeI(2)) & (cgv>=1));
  otherwise
    error('Unsupported extension.')
end
% rgv = round(rgv);
% cgv = round(cgv);
if isvector(rgv)
  p = img(uint32(round(rgv)),uint32(round(cgv)));
else
  szP = size(rgv);
  szI = size(img);
  idx = rgv+(cgv-1)*szI(1);
  p = reshape(img(uint64(round(idx(:)))),szP);
end
return
