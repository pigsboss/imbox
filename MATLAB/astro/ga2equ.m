function varargout = ga2equ(varargin)
%GA2EQU Galactic coordinates to equatorial coordinates (J2000.0)
narginchk(1,2)
switch nargin
  case 1
    l = varargin{1}{1};
    b = varargin{1}{2};
  case 2
    l = varargin{1};
    b = varargin{2};
end
pole_ra = 192.859508;
pole_dec = 27.128336;
posangle = 122.932-90.0;
ra = atan2((cosd(b).*cosd(l-posangle)),...
  (sind(b)*cosd(pole_dec)-sind(pole_dec)*cosd(b).*sind(l-posangle)))*180/pi...
  +pole_ra;
dec = asind(cosd(pole_dec)*cosd(b).*sind(l-posangle) + sind(b)*sind(pole_dec));
ra = mod(ra,360);
if any(l(:)<0)
  ra(ra>180)=ra(ra>180)-360;
end
nargoutchk(0,2)
if nargout>1
    varargout{1} = ra;
    varargout{2} = dec;
else
    varargout{1}={ra,dec};
end
return