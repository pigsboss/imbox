close all
clear
ecc = 0.8;
phi = pi/4;
b = 1;
N = 100;
theta = ((1:N)/N*2*pi)';
nlvl = 0.01;
rho = normrnd(b./sqrt(1 - ecc^2*cos(theta(:)-phi).^2),nlvl);
rho(1:40)=normrnd(rho(1:40),0.3);
x = rho(:).*cos(theta(:));
y = rho(:).*sin(theta(:));
scatter(x,y,'*')
wgt = ones(size(rho));
[~,~,~,~,eof] = ellifitbf(rho,theta,wgt,10);
% ellifitbf(rho,theta,wgt,10);
ellifitbf(rho,theta,1./eof.^2,10);
