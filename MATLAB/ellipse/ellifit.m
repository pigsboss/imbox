function [xc,yc,xr,yr,phi_major,EoF]=ellifit(x,y)
% center of the ellipse
xc = mean(x);
yc = mean(y);
% shift to the center
x = x-xc;
y = y-yc;
% distance to the center
r = sqrt(x.^2 + y.^2);
% phase angle, in [-pi/2, pi/2)
phi = atan2(y,x);
phi(phi >= (pi/2)) = phi(phi >= (pi/2)) - pi;
phi(phi < (-pi/2)) = phi(phi < (-pi/2)) + pi;
% phase angle of major axis
phi_major = sum((r.^10).*phi)/sum(r.^10);
% rotate major axis to x-axis
X = x.*cos(phi_major)+y.*sin(phi_major);
Y = y.*cos(phi_major)-x.*sin(phi_major);
x = X;y = Y;clear X Y
x2 = sum(x.^2);y2 = sum(y.^2);x4 = sum(x.^4); y4 = sum(y.^4);x2y2 = sum((x.^2).*(y.^2));
xr = (x2*y4 - y2*x2y2)/(x4*y4 - x2y2^2);
yr = (y2*x4 - x2*x2y2)/(y4*x4 - x2y2^2);
if xr>0 && yr>0
    xr = sqrt(1/xr);
    yr = sqrt(1/yr);
else
    disp('bad fit.')
%     display(x)
%     display(y)
    xr = (max(x) - min(x))/2;
    yr = (max(y) - min(y))/2;
end
clear x2 y2 x2y2
EoF = std((x/xr).^2+(y/yr).^2);
return
