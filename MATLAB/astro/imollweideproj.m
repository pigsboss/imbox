function [phi,theta]=imollweideproj(x,y)
y=y*sqrt(2)/pi*2;
lambda=asin(y/sqrt(2));
theta=asin((2*lambda+sin(2*lambda))/pi);
if abs(abs(lambda)-pi/2)<=eps
    phi=0;
else
    phi=pi*x/2/sqrt(2)/cos(lambda);
end
% phi=phi+pi;
% theta=pi/2-theta;
return