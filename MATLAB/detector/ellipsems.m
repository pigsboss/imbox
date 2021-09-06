function varargout = ellipsems(varargin)
%ELLIPSEMS detect multi-scale blobs in given image and use elliptical areas
% to fit them.
%
%Syntax:
%   ellipsems
%   ellipsems(I)
%   ellipsems(I,scale)
%   ellipsems(I,scale,threshold)
%   ellipsems(..., 'polyopts', line_properties, elliopts, line_properties)
%   [bmes, bimg] = ellipsems(...)
%
%Description:
%   I is an image contains bright blobs. Default is the image in the current
%   figure object.
%
%   scale is the size of blobs to be extracted.
%
%   threshold is used to rule out artificials.
%
%   polyopts is passed to MATLAB built-in function line() to draw border polygons
%   of extracted blobs.
%
%   elliopts is passed to MATLAB built-in function line() to draw border
%   ellipses of extracted blobs.
%   
%Returns:
%   bmes is structure of blob measures. Flux, centre location, elliptical
%   parameters of each blob are included.
%
%   bimg is structure of blob images. Foreground image, background image,
%   mask and border polygon vertices are included.
%
%Overall Design:
% +-----------+                +---------------------------+
% | I (image) |                | B (blobs map)             |
% | scale     +-----BLOBMS---->| L (multiscale laplacians) |
% | threshold | Detect blobs.  | other variables           |
% +-----------+                +-------------+-------------+
%                                            |
%                                            v
% +--BORDER---------------+    +--BLOBFILTER---------------+
% |Detect border vertices |    |Remove invalid blobs:      |
% |of blobs.              |    |1. Close to the boundaries |
% |                       |<---+2. ...                     |
% |                       |    |                           |
% +-----------------------+    +---------------------------+
%
%See Also:
%border, blobms, ellifitbf, fajita
%


%Detect multiscale blobs:
[B, L, I, xgv, ygv, scale, varargin] = blobms(varargin{:});

%Parse optional inputs:
[EDGE, PRECISION, polyopts, elliopts] = parse_inputs(varargin{:});

%Construct xy meshgrid:
[X, Y] = meshgrid(xgv, ygv);

%Filter blobs:
[B, idx, NB] = blobfilter(X,Y,B);

%Initialize blob structures:
[bmes,bimg] = initblobstruct(NB);

if nargout == 0
  figure;imagesc(I);axis image;colormap('gray')
end

for k=1:NB
  %Find scale index of the k-th blob:
  [~, scidx] = min(abs(scale - B(idx(k))));

  %Find border vertices of the k-th blob:
  [rhob,phib,wgtb,msk,rect] =...
    border(X(idx(k)),Y(idx(k)),L(:,:,scidx),B(idx(k)),EDGE,PRECISION);

  img = imcrop(I,rect);

  %Elliptical fitting of the border vertices:
  [ell, geoparam] = ellifitbf(rhob,phib,wgtb);

  xb = rhob.*cos(phib)+X(idx(k));
  yb = rhob.*sin(phib)+Y(idx(k));
  ell(1) = ell(1)+X(idx(k));
  ell(2) = ell(2)+Y(idx(k));

  %Set blob parameters and images.

  %x and y coordinates of border polygon, i.e., a sequence of vertices.
  bimg(k).xbdply = xb; 
  bimg(k).ybdply = yb;

  %Blob mask is a true-and-false (1-and-0) matrix.
  bimg(k).mask = msk;
  %Blob image is a subimage of the input image.
  bimg(k).image = img;

  %Background image:
  bimg(k).bgrdimg = min(img,fajita(img.*double(~msk),(~msk),'linear'));
  %Foreground image:
  bimg(k).fgrdimg = max(0,img - bimg(k).bgrdimg);

  %Maximum pixel value of the blob:
  bmes.peak(k) = max(bimg(k).fgrdimg(:));
  blobX = imcrop(X,rect);
  blobY = imcrop(Y,rect);
  bmes.flux(k) = sum(bimg(k).fgrdimg(:));

  %Subpixel coordinate of blob center.
  if bmes.flux(k)>0
    bmes.xctrsubpix(k) = sum(blobX(:).*bimg(k).fgrdimg(:)) /...
      bmes.flux(k);
    bmes.yctrsubpix(k) = sum(blobY(:).*bimg(k).fgrdimg(:)) /...
      bmes.flux(k);
  end
  %Coordinate of blob center.
  bmes.xctr(k) = X(idx(k));
  bmes.yctr(k) = Y(idx(k));
  bmes.scale(k) = B(idx(k));

  %The rectangle is a 4-element vector defined as [x y width height].
  bmes.rect(k,:) = rect;
  bmes.elliparam(k,:) = ell;
  bmes.posangle(k) = ell(5);
  bmes.semiminor(k) = ell(4);
  bmes.semimajor(k) = ell(3);
  bmes.eccentricity(k) = geoparam(2);
end
%draw polygon and ellipse borders of blobs
if nargout == 0
  for k = 1:NB
    xb = bimg(k).xbdply; 
    yb = bimg(k).ybdply;
    ell = bmes.elliparam(k,:);
    line([xb(:);xb(1)],[yb(:);yb(1)],polyopts{:})
    ellidraw(ell,elliopts{:})
  end
  drawnow
end
varargout{1} = bmes;
varargout{2} = bimg;
return

function [B,idx,NB] = blobfilter(X,Y,B)
%remove blobs close to boundaries.
idx = find(B>0);
NB = length(idx);
sizeI = size(B);
for k = 1:NB
  radius = pi*B(idx(k));
  if X(idx(k))<radius || X(idx(k))>(sizeI(2)-radius) || ...
      Y(idx(k))<radius || Y(idx(k))>(sizeI(1)-radius)
    B(idx(k)) = 0;
  end
end
%update blob number and index
idx = find(B>0);
NB = length(idx);
return

function [EDGE, PRECISION, polyopts, elliopts] = parse_inputs(varargin)
EDGE_d = 0;
PRECISION_d = 1e-3;
polyopts_d = {'color','r','linestyle','-'};
elliopts_d = {'color','g','linestyle','-'};
[EDGE, PRECISION, polyopts, elliopts] =...
  getoptval(varargin, {'edge','precision','polyopts','elliopts'});
if isempty(EDGE), EDGE = EDGE_d; end
if isempty(PRECISION), PRECISION = PRECISION_d; end
if isempty(polyopts), polyopts = polyopts_d; end
if isempty(elliopts), elliopts = elliopts_d; end
return
