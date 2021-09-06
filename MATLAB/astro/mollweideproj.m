function varargout=mollweideproj(phi,theta,z)
%MOLLWEIDEPROJ Mollweide projection.
% x: from -2*sqrt(2) to 2*sqrt(2)
% y: from -sqrt(2) to sqrt(2)
% theta: from -pi/2 to pi/2, 
% phi: from -pi to pi
%
% [x,y]=mollweideproj(phi,theta)
% (x,y) is mollweide projection of (phi,theta) on unit sphere.
%
% [x,y,z,c]=mollweideproj(phi,theta,z)
% (x,y) is mollweide projection of (phi,theta) on unit sphere. z(x,y) =
% z(phi,theta). c is coverage (redundancy) of the projection.

dim=size(theta);
x=zeros(size(theta));
y=zeros(size(theta));
for n=1:numel(x)
    lambda=newtonraphson(theta(n));
    x(n)=2*sqrt(2)*phi(n)*cos(lambda)/pi;
    y(n)=sqrt(2)*sin(lambda);
end
x=reshape(x,dim);
y=reshape(y,dim);
y=y*pi*0.5/sqrt(2);
if nargin==2
    varargout(1)={x};
    varargout(2)={y};
    return
else
    X=x;
    Y=y;
    DIM=dim;
    DIM(1)=round(0.4*DIM(1));
    DIM(2)=round(0.8*DIM(2));
    y=round((y-min(min(y)))*(DIM(1)-1)/(max(max(y))-min(min(y)))+1);
    x=round((x-min(min(x)))*(DIM(2)-1)/(max(max(x))-min(min(x)))+1);
    x=x(:);
    y=y(:);
    z=z(:);
    c=zeros(DIM);
    Z=zeros(DIM);
    for n=1:numel(x)
        c(y(n),x(n))=c(y(n),x(n))+1;
        Z(y(n),x(n))=Z(y(n),x(n))+z(n);
    end
    Z=Z./(c+double(c<=0))-double(c<=0);
    Z=flipud(Z);
    c=flipud(c);
    varargout(1)={X};
    varargout(2)={Y};
    varargout(3)={Z};
    varargout(4)={c};
end
return

function lambda=newtonraphson(theta)
%NEWTONRAPHSON Newton-Raphson iteration for solving equation in the main
%function.
if abs(abs(theta)-pi/2)<=eps
    lambda=sign(theta)*pi/2;
    return
end
lambda=theta;
delta=pi;
while abs(delta)>1e-12
    delta=(2*lambda+sin(2*lambda)-pi*sin(theta))/(2+2*cos(2*lambda));
    lambda=lambda-delta;
end
return