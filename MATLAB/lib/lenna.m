function img = lenna(varargin)
%LENNA Read Lenna portrait.
%SYNTAX:
% lenna         reads the original color image.
% lenna(N)      reads the color image and resizes it to NxN.
% lenna(N,type) reads the color image, resizes it to NxN, and converts it to the
%               specified type. Supported type: 'r', 'g', 'b', 'rgb', 'gray',
%               'bw'.
% lenna(type)   reads the color image and converts it to the specified type.

narginchk(0,2)
switch nargin
  case 0
    N = [];
    type = '';
  case 1
    if isnumeric(varargin{1})
      N = varargin{1};
      type = '';
    else
      N = [];
      type = varargin{1};
    end
  case 2
    N = varargin{1};
    type = varargin{2};
end
img = imread('Lenna.png');
if ~isempty(N)
  if isscalar(N)
    N = [N,N];
  end
  img = imresize(img,N);
end
type = lower(type);
[H,W,~] = size(img);
switch type
  case {''}
  case {'gray','grey','bw'}
    img = mean(double(img),3)/255;
  otherwise
    if isempty(strfind(type,'r'))
      img(:,:,1) = uint8(zeros(H,W));
    end
    if isempty(strfind(type,'g'))
      img(:,:,2) = uint8(zeros(H,W));
    end
    if isempty(strfind(type,'b'))
      img(:,:,3) = uint8(zeros(H,W));
    end
end
return