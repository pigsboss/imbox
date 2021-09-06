function varargout = equ2ga(varargin)
%GA2EQU Galactic coordinates to equatorial coordinates (J2000.0)
narginchk(1,2)
switch nargin
  case 1
    ra = varargin{1}{1};
    dec = varargin{1}{2};
  case 2
    ra = varargin{1};
    dec = varargin{2};
end
pole_ra = 192.859508;
pole_dec = 27.128336;
posangle = 122.932-90.0;
b = asind(cosd(ra-pole_ra).*cosd(dec)*cosd(pole_dec)+sind(dec)*sind(pole_dec));
l = atan2((sind(dec)-sind(b)*sind(pole_dec)),...
  (cosd(dec).*sind(ra-pole_ra)*cosd(pole_dec)))*180/pi+posangle;
l = mod(l,360);
if any(ra(:)<0)
  l(l>180)=l(l>180)-360;
end
nargoutchk(0,2)
if nargout>1
    varargout{1} = l;
    varargout{2} = b;
else
    varargout{1}={l,b};
end
return