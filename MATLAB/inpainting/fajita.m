function varargout = fajita(f,m,ORDER,PRECISION,MAXLOOP)
%FAJITA Image inpainting by Fast AdJacent InpainTing Agorithm.
%
%SYNTAX
% [f,g] = fajita(f,m,ORDER,PRECISION,MAXLOOP)
%
%INPUT
% f                    is the input image with sockets.
%
% m         (optional) is the image mask. Default is logical(f > 0).
%
% ORDER     (optional) is order of discrete Laplacians. Default is 0.
%                      FAJITA uses direct inpainting from nearest pixels for
%                      0-order.
%                      For higher orders FAJITA uses total variation
%                      minimization iterations.
%
% PRECISION (optional) is criteria of convergence. Default is 1e-3.
%
% MAXLOOP   (optional) is max number of loops. Default is 100.
%
%RETURN
% f is inpainted image.
%
% g is gradient of the steepest descent algorithm.

% global VERBOSE
narginchk(1,5)
imgname = '';
if ischar(f)
  imgname = f;
  f = imread(imgname);
  imgname = [imgname,'.ipt.png'];
end
fi = f;
typef = str2func(class(f));
f = double(f);
ORDER_d = 0;
PRECISION_d = 1e-3;
MAXLOOP_d = 100;
m_d = (f>0);
switch nargin
  case 1
    m = m_d;
    ORDER = ORDER_d;
    PRECISION = PRECISION_d;
    MAXLOOP = MAXLOOP_d;
  case 2
    if ischar(m) || isscalar(m)
      ORDER = m;
      m = m_d;
    else
      ORDER = ORDER_d;
    end
    PRECISION = PRECISION_d;
    MAXLOOP = MAXLOOP_d;
  case 3
    PRECISION = PRECISION_d;
    MAXLOOP = MAXLOOP_d;
  case 4
    MAXLOOP = MAXLOOP_d;
end
clear ORDER_d PRECISION_d MAXLOOP_d m_d

if ischar(ORDER)
  switch ORDER
    case {'nearest'}
      ORDER = 0;
    case {'linear','bilinear'}
      ORDER = 1;
    case {'cubic','bicubic'}
      ORDER = 3;
  end
end
ff = fajitanearest(f,m);
g = dlt2(dlt2(f,ORDER*2),ORDER*2);
if ORDER > 0
  g2 = dlt2(dlt2(ff,ORDER*2),ORDER*2);
  if max(abs(g(~m))) > max(abs(g2(~m)))
      f = fajitanearest(f,m);
  end
  meanf = mean(f(:));
  delta = PRECISION * meanf;
  for k = 1:MAXLOOP
    l = dlt2(f,ORDER*2);
    g = dlt2(l,ORDER*2);
    h = dlt2(g,ORDER*2);
    gamma = sum(l(:).*h(:)) / sum(h(:).^2);
    d = rand(1)*2*gamma*g.*double(~m);
    f = f - d;
    if max(d(:)) <= delta
%       if VERBOSE >= 1
%         disp(['Convergence criteria has been met by precision at step ',...
%           int2str(k)])
%       end
      break
    end
  end
else
  f = ff;
end
f = typef(f);
switch nargout
  case 0
    switch ndims(f)
      case 2
        figure('Name','Input image');
        imagesc(double(fi));axis image
        figure('Name',['Inpainted image, order=',int2str(ORDER)]);
        imagesc(double(f));axis image
      case 3
        figure('Name','Input image');
        imshow(fi);
        figure('Name',['Inpainted image, order=',int2str(ORDER)]);
        imshow(f);
        if ~isempty(imgname)
          imwrite(f,imgname,'png')
        end
    end
  otherwise
    varargout{1} = f;
    varargout{2} = g;
end
return

function Y = fajitanearest(X,M)
%FAJITANEAREST FAst adJacent InpainTing Algorithm.
%INPUTS:
%  X           is image with sockets.
%  M (boolean) is mask of the image X. false in M labels sockets.
%RETURN:
%  Y           is inpainted image.
typex = str2func(class(X));
Y = double(X);
while any(~M(:))
  [Y,M] = adjaverage(Y,M);
end
Y = typex(Y);
return

function [Y,M] = adjaverage(X,M)
%ADJAVERAGE Adjacent pixel averaging.
sft = [-1,-1; -1,0; -1,1; 0,-1; 0,1; 1,-1; 1,0; 1,1];
wgt = sqrt(sum(sft.^2,2));
sizeX = size(X);
ndX = length(sizeX);
sizeM = size(M);
ndM = length(sizeM);
if ndM==2 && ndX==3
  W = zeros(sizeX);
  for k = 1:sizeX(3)
    W(:,:,k) = M;
  end
  M = logical(W);
end
W = zeros(size(M));
Y = zeros(size(X));
for k = 1:8
  W = W+zeroshift(M,sft(k,:))*wgt(k);
  Y = Y+zeroshift(X,sft(k,:))*wgt(k);
end
Y((W>0.7)&(~M)) = Y((W>0.7)&(~M))./W((W>0.7)&(~M));
Y(M) = X(M);
M = (W>0.7);
return
