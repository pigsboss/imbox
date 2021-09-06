function varargout = ellifitbf(rho,theta,wgt,MAXLOOPS)
%ELLIFITBF Elliptical fitting by brute force.
narginchk(2,4)
switch nargin
    case 2
        wgt = ones(size(rho));
        MAXLOOPS = 10;
    case 3
        if isscalar(wgt)
            MAXLOOPS = wgt;
            wgt = ones(size(rho));
        else
            MAXLOOPS = 10;
        end
end
N = 20;
t = 2;
emax = 1;
emin = 0;
pmax = pi;
pmin = 0;
J = zeros(MAXLOOPS,1);
for l=1:MAXLOOPS
    egv = ((1:N) - 1) / N * (emax-emin) + emin;
    pgv = ((1:N) - 1) / N * (pmax-pmin) + pmin;
    [S,b,ecc,phi] = ellmsemat(rho,theta,wgt,egv,pgv);
    S = sdwt2(S,1);
    [J(l),idx] = min(S(:));
    ectr = ecc(idx);
    pctr = phi(idx);
    bctr = b(idx);
    emax = min(ectr + (egv(N)-egv(1))*0.5/t,1);
    emin = max(ectr - (egv(N)-egv(1))*0.5/t,0);
    pmax = min(pctr + (pgv(N)-pgv(1))*0.5/t,pi);
    pmin = max(pctr - (pgv(N)-pgv(1))*0.5/t,0);
    if l>1
        if abs(J(l)-J(l-1)) <= 1e-8
            break
        end
    end
end
actr = bctr/sqrt(1-ectr^2);
ell = [0,0,actr,bctr,pctr];
nargoutchk(0,6)
switch nargout
    case 0
        ellidraw(ell,100,'color','r')
        disp(['eccentricity:',num2str(ectr)])
        disp(['semi-minor axis:',num2str(bctr)])
        disp(['polar angle of major axis:',num2str(pctr)])
    case 1
        varargout{1} = ell;
    case 2
        varargout{1} = ell;
        varargout{2} = [bctr;ectr;pctr];
    case 3
        varargout{1} = bctr;
        varargout{2} = ectr;
        varargout{3} = pctr;
    case 4
        varargout{1} = bctr;
        varargout{2} = ectr;
        varargout{3} = pctr;
        varargout{4} = J;
    case 5
        varargout{1} = bctr;
        varargout{2} = ectr;
        varargout{3} = pctr;
        varargout{4} = J;
        varargout{5} = ellmse(rho,theta,bctr,ectr,pctr);
    case 6
        varargout{1} = 0;
        varargout{2} = 0;
        varargout{3} = actr;
        varargout{4} = bctr;
        varargout{5} = pctr;
        varargout{6} = min(J(J(:)>0));
end
return

function [S,b,ecc,phi] = ellmsemat(rho,theta,wgt,egv,pgv)
[ecc,phi] = meshgrid(egv,pgv);
S = zeros(size(ecc));
b = zeros(size(ecc));
wgts = sum(wgt(:));
for l = 1:numel(ecc)
    d = 1./sqrt(1-ecc(l)^2*cos(theta(:)-phi(l)).^2);
    b(l) = sum(rho(:).*wgt(:)) ./ sum(d(:).*wgt(:));
    S(l) = sum((b(l).*d - rho(:)).^2.*...
        wgt(:))/wgts;
end
return

function S = ellmse(rho,theta,b,ecc,phi)
S = abs(b./sqrt(1-ecc^2*cos(theta-phi).^2) - rho);
return

% 
% function G = ellmsegrad(rho,theta,wgt,b,ell,phi)
% wgts = sum(wgt(:));
% d = (1 - ell^2*cos(theta(:)+phi).^2);
% gb = 2*sum(((b./sqrt(d)-rho(:))./sqrt(d)).*wgt(:))./wgts;
% gp = 2*b*ell^2*sum((b./sqrt(d)-rho(:)).*(d.^(-3/2)).*cos(theta(:)+phi).*...
%     sin(theta(:)+phi).*wgt(:))./wgts;
% ge = -2*b*ell*sum((b./sqrt(d)-rho(:)).*(d.^(-3/2)).*cos(theta(:)+phi).*...
%     cos(theta(:)+phi).*wgt(:))./wgts;
% G = [gb;ge;gp];
% return
