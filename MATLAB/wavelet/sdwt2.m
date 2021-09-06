function [c,w] = sdwt2(varargin)
%SDWT2 sparse discrete wavelet transform for 2D array.

% parsing input arguments:
J = [];
h = [1/16 1/4 3/8 1/4 1/16];
narginchk(1,3)
switch nargin
    case 1
        x=varargin{1};
    case 2
        x=varargin{1};
        J=varargin{2};
    case 3
        x=varargin{1};
        J=varargin{2};
        h=varargin{3};
end
if isvector(h)
    h = h(:)*(h(:)');
end

% setup initial parameters:
[d1,d2]=size(x);
if isempty(J) % the max_scale is determined automatically
    J=lastpow2(min(d1,d2));
end
if J<=0
    c=x;
    w=zeros(size(c));
    return
end
sizeh = size(h);
lenh = prod(sizeh);
[IX,IY] = meshgrid(1:sizeh(2), 1:sizeh(1));
IX = IX - ceil((sizeh(2)+1)/2);
IY = IY - ceil((sizeh(1)+1)/2);
c = x;
w = zeros(d1,d2,J);
for k = 0:(J-1)
    cnext = zeros(d1,d2);
    for l = 1:lenh
%         disp([IY(l),IX(l)])
        cnext = cnext + symmshift(c,[IY(l),IX(l)])*h(l);
    end
    w(:,:,k+1) = c - cnext;
    c = cnext;
    IX = IX*2;
    IY = IY*2;
end
return

function y=symmshift(x,vec)
%SYMMSHIFT 2D array shift with symmetric boundary.
sizex = size(x);
vec = mod(vec, 2*sizex);
xud = flipud(x);
xlr = fliplr(x);
xrot = rot90(x,2);
if vec(1) < sizex(1)
    if vec(2) < sizex(2)
        y = [xrot((sizex(1)-vec(1)+1):sizex(1), (sizex(2)-vec(2)+1):sizex(2))...
            xud((sizex(1)-vec(1)+1):sizex(1), 1:(sizex(2)-vec(2)));...
            xlr(1:(sizex(1)-vec(1)), (sizex(2)-vec(2)+1):sizex(2)),...
            x(1:(sizex(1)-vec(1)), 1:(sizex(2)-vec(2)))];
        return
    else
        vec(2) = vec(2) - sizex(2);
        y = [xud((sizex(1)-vec(1)+1):sizex(1), (sizex(2)-vec(2)+1):sizex(2))...
            xrot((sizex(1)-vec(1)+1):sizex(1), 1:(sizex(2)-vec(2)));...
            x(1:(sizex(1)-vec(1)), (sizex(2)-vec(2)+1):sizex(2)),...
            xlr(1:(sizex(1)-vec(1)), 1:(sizex(2)-vec(2)))];
        return
    end
else
    vec(1) = vec(1) - sizex(1);
    if vec(2) < sizex(2)
        y = [xlr((sizex(1)-vec(1)+1):sizex(1), (sizex(2)-vec(2)+1):sizex(2))...
            x((sizex(1)-vec(1)+1):sizex(1), 1:(sizex(2)-vec(2)));...
            xrot(1:(sizex(1)-vec(1)), (sizex(2)-vec(2)+1):sizex(2)),...
            xud(1:(sizex(1)-vec(1)), 1:(sizex(2)-vec(2)))];
        return
    else
        vec(2) = vec(2) - sizex(2);
        y = [x((sizex(1)-vec(1)+1):sizex(1), (sizex(2)-vec(2)+1):sizex(2))...
            xlr((sizex(1)-vec(1)+1):sizex(1), 1:(sizex(2)-vec(2)));...
            xud(1:(sizex(1)-vec(1)), (sizex(2)-vec(2)+1):sizex(2)),...
            xrot(1:(sizex(1)-vec(1)), 1:(sizex(2)-vec(2)))];
        return
    end
end
