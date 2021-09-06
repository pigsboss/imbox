function varargout = ellidraw(varargin)
%ELLIDRAW draw ellipse define by center (xc,yc), major radius xr, minor
%radius yr, and phase angle of major axis phi_major.
%
%Syntax:
%   ellidraw(xc, yc, xr, yr, phi_major)
%   ellidraw(ElliParms)
%   ellidraw(..., NumPts)
%   ellidraw(..., LineOpts)
%   [X, Y] = ellidraw(...)
%Input arguments:
%   ElliParms is a 5 x NumEll matrix. Each column is a group of elliptical parameters.
%   NumPts    is the number of points used to draw an ellipse.
%   LineOpts  are options to pass to matlab builtin line function.
%Returns:
%   X, Y      are vertices of polygons used to approximate ellipses.

[xc,yc,xr,yr,phi_major,N,opts,NumEll] = parse_inputs(varargin{:});
phi = (2*pi*(1:N)/N)';
rho = (ones(N,1)*(xr.*yr))./sqrt((cos(phi)*yr).^2 + (sin(phi)*xr).^2);
x = rho.*(cos(phi)*ones(1,NumEll));
y = rho.*(sin(phi)*ones(1,NumEll));
X = x.*(ones(N,1)*cos(phi_major))-y.*(ones(N,1)*sin(phi_major));
Y = y.*(ones(N,1)*cos(phi_major))+x.*(ones(N,1)*sin(phi_major));
clear x y
X = X+(ones(N,1)*xc);
Y = Y+(ones(N,1)*yc);
nargoutchk(0,2)
switch nargout
    case 0
        for k = 1:NumEll
            if isempty(opts)
                linvar = {[X(:,k);X(1,k)],[Y(:,k);Y(1,k)]};
            else
                linvar = [{[X(:,k);X(1,k)]},{[Y(:,k);Y(1,k)]},opts{:}];
            end
            line(linvar{:})
        end
    case 2
        varargout{1} = X;
        varargout{2} = Y;
end
return

function [xc,yc,xr,yr,phi_major,N,opts,NumEll] = parse_inputs(varargin)
xc = []; xc_d = 0;
yc = []; yc_d = 0;
xr = []; xr_d = 2;
yr = []; yr_d = 1;
phi_major = []; phi_major_d = 0;
N = []; N_d = 40;
opts = [];
if ~isvector(varargin{1})
    xc = varargin{1}(1,:);
    yc = varargin{1}(2,:);
    xr = varargin{1}(3,:);
    yr = varargin{1}(4,:);
    phi_major = varargin{1}(5,:);
    if nargin >= 2
        if isnumeric(varargin{2})
            N = varargin{2};
            opts = {varargin(3:nargin)};
        else
            opts = {varargin(2:nargin)};
        end
    end
else
    if nargin == 1
        xc = varargin{1}(1);
        yc = varargin{1}(2);
        xr = varargin{1}(3);
        yr = varargin{1}(4);
        phi_major = varargin{1}(5);
    elseif isvector(varargin{2}) && isnumeric(varargin{2}) && nargin >= 5
        xc = varargin{1};
        yc = varargin{2};
        xr = varargin{3};
        yr = varargin{4};
        phi_major = varargin{5};
        if nargin >= 6
            if isnumeric(varargin{6})
                N = varargin{6};
                opts = {varargin(7:nargin)};
            else
                opts = {varargin(6:nargin)};
            end
        end
    else
        xc = varargin{1}(1);
        yc = varargin{1}(2);
        xr = varargin{1}(3);
        yr = varargin{1}(4);
        phi_major = varargin{1}(5);
        if isnumeric(varargin{2})
            N = varargin{2};
            opts = {varargin(3:nargin)};
        else
            opts = {varargin(2:nargin)};
        end
    end
end
if isempty(xc), xc = xc_d;end
if isempty(yc), yc = yc_d;end
if isempty(xr), xr = xr_d;end
if isempty(yr), yr = yr_d;end
if isempty(phi_major), phi_major = phi_major_d;end
if isempty(N),  N = N_d;end
NumEll = length(xr); % number of ellipses
xc = xc(:)';
yc = yc(:)';
xr = xr(:)';
yr = yr(:)';
phi_major = phi_major(:)';
N = round(max(N,3));
return