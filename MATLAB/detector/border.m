function [rhoBMask,phiBMask,wgtBMask,pBMask,rectP] =...
  border(xc,yc,L,scale,varargin)
%BORDER border of dark blob
%
%Syntax
%   xc and yc is center of the current blob.
%   L         is the laplacian.
%   scale     is scale of the current blob.
%   EDGE      is width of edge beyond the detected vertex.
%   PRECISION is passed to findroot().
%
%Algorithms
% +--------------------------------------------------------+
% |1.    Select a patch centred at the given blob.         |
% |                                                        |
% |2.    Transform Laplacian from cartesian framework to   |
% |      polar coordinate system.                          |
% |                                                        |
% |3.    For each radius L(rho|phi=phi_k), find the border |
% |      vertex:                                           |
% |                                                        |
% |3.1   Find the first root of L(rho) = 0. This root is   |
% |      an unstable solution. Refer to this root as rho_a.|
% |                                                        |
% |3.2   Find the first concave point of L(rho) on the     |
% |      right of rho_a (rho>=rho_a). The concave point is |
% |      defined as:                                       |
% |                                                        |
% |3.2.a       L(rho)<=0                                   |
% |                                                        |
% |3.2.b       L'(rho)==0                                  |
% |                                                        |
% |3.2.c       L"(rho)>=0                                  |
% |                                                        |
% |      Refer to the concave point as rho_b. rho_b is the |
% |      vertex requested.                                 |
% |                                                        |
% +--------------------------------------------------------+
%
%See also
%blobms, ellipsems
%
%References:
% 1. imbox/Notes/ellipsems.pptx

%Set parameters and initial values
narginchk(4,6)
PRECISION_d = 1e-1;
EDGE_d = 0;
EDGE = [];
PRECISION = [];
switch nargin
  case 5
    if varargin{1} >=1
      EDGE = varargin{1};
    else
      PRECISION = varargin{1};
    end
  case 6
    EDGE = varargin{1};
    PRECISION = varargin{2};
end
if isempty(EDGE); EDGE = EDGE_d;end
if isempty(PRECISION); PRECISION = PRECISION_d;end

%Select the patch.
%rB is maximum radius of a blob with the given scale.
%Use symmetric boundary to keep the blob centre on the centre of the
%patch.
szL = size(L);
rB = round(4*pi*scale);
xgv = (xc-rB):(xc+rB);
ygv = (yc-rB):(yc+rB);
[pL,rgv,cgv] = subim(L,ygv,xgv,'sym');
szP = size(pL);
rectP = [cgv(1) rgv(1) max(cgv)-min(cgv) max(rgv)-min(rgv)];

%Patch of laplacian in polar coordinate system.
[ppL,rho,phi] = xy2polar(rect2square(pL));

%Patch of blob mask in polar coordinate system.
%A mask is a logical matrix. Its elements tell whether the corresponding
%pixels are inside the current blob. For example, if I and M are the image
%and the blob mask, (I & M) is the image of the current blob.
ppBMask = zeros(size(ppL));
[dphi,drho] = size(ppL);
rhoBMask = zeros(dphi,1);
phiBMask = phi(:,1);
wgtBMask = zeros(dphi,1);

for k=1:dphi
  %Smooth the slope of the blob along the k-th direction.
  ppL(k,:) = gsmooth(ppL(k,:),scale,'symmetric');

  %Normalized derivative of ppL:
  dppL = normalize(imconv(ppL(k,:),[0.5,0,-0.5],'symmetric'));
  %Normalized 2nd-order derivative of ppL:
  d2ppL = normalize(imconv(ppL(k,:),[1,-2,1],'symmetric'));

  %Special points:
  %1. rho_a: first zero point (unstable solution of L(rho) == 0).
  %2. rho_b: first concave point (solution of L'(rho) == 0).
  %3. rho_c: second zero point (stable solution of L(rho) == 0).
  %4. rho_d: first convex point on the right of rho_b.
  %5. rho_e: second concave point.
  %6. rho_f: second zero point (unstable solution of L(rho) == 0).
  
  %Find all roots of L(rho) == 0:
  rho_0 = max(min(findroot(1:drho,ppL(k,:),PRECISION),drho),1);
  if rho_0 == drho
    rhoBMask(k) = rho(drho);
    wgtBMask(k) = 0;
    ppBMask(k,:) = double((1:drho) <= drho);
    continue
  end
  if isempty(rho_0)
    warning('imbox:border:noRootFound',...
      'No root of L(x)=0 found with current precision.')
    rhoBMask(k) = rho(drho);
    wgtBMask(k) = 0;
    ppBMask(k,:) = double((1:drho) <= drho);
    continue
  end
  rho_a = rho_0(1);
  
  %Find all roots of L'(rho) == 0:
  if std(dppL(ceil(rho_a):drho)) >= eps
    rho_1 = max(min(findroot(ceil(rho_a):drho,dppL(ceil(rho_a):drho),...
      PRECISION),drho),1);
  else
    warning('imbox:border:noRootFound',...
      'No root of L''(x)=0 found with default precision.')
    rhoBMask(k) = rho(drho);
    wgtBMask(k) = 0;
    ppBMask(k,:) = double((1:drho) <= drho);
    continue
  end
  if isempty(rho_1)
    warning('imbox:border:noRootFound',...
      ['No root of L''(x)=0 found with current precision.',...
      'Try default precision.'])
    rho_1 = max(min(findroot(ceil(rho_a):drho,dppL(ceil(rho_a):drho),...
      PRECISION_d),drho),1);
    if isempty(rho_1)
      warning('imbox:border:noRootFound',...
        'No root of L''(x)=0 found with default precision.')
      rhoBMask(k) = rho(drho);
      wgtBMask(k) = 0;
      ppBMask(k,:) = double((1:drho) <= drho);
      continue
    end
  end
  %Merit function of concave points:
  %s(x) = f"(x) - f(x)
  score_ccv = d2ppL(round(rho_1))-ppL(k,round(rho_1));
  rho_ccv = rho_1(score_ccv>0);
  if isempty(rho_ccv)
    warning('imbox:border:noRootFound',...
      'No concave point of L(x) found with current precision.')
    rhoBMask(k) = rho(drho);
    wgtBMask(k) = 0;
    ppBMask(k,:) = double((1:drho) <= drho);
    continue
  end
  rho_b = rho_ccv(1);
  rhoBMask(k) = rho(1,min(round(rho_b+EDGE),drho));
  wgtBMask(k) = 1;
  ppBMask(k,:) = double((1:drho) <= (rho_b+EDGE));
end
rhoBMask = gsmooth(rhoBMask,pi*scale,'circular');

%Rule out extension of the patch.
%The patch is extended using symmetric boundary to keep the blob centre
%on the geometric centre of the patch, i.e., the origin of the polar
%coordinate system. This is mandatory to run xy2polar.
pBMask = square2rect(logical(polar2xy(ppBMask) > 0),szP);
xidx_head = find(xgv>=1,1,'first');
xidx_tail = find(xgv<=szL(2),1,'last');
yidx_head = find(ygv>=1,1,'first');
yidx_tail = find(ygv<=szL(1),1,'last');
pBMask = pBMask(yidx_head:yidx_tail,xidx_head:xidx_tail);
return