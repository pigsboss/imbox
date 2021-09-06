function [sigma,sigma_sigma]=mssigma(N,p)
%MSSIGMA standard deviations for multi-scale (b-spline wavelet) data
if N <= 10
    sigma=[...
      0.890544;...
      0.200866;...
      0.085688;...
      0.041390;...
      0.020598;...
      0.010361;...
      0.005256;...
      0.002706;...
      0.001412;...
      0.000784];
    sigma=sigma(1:N);
    sigma_sigma=[];
    return
end
if nargin<2
    sigma=sigmadec(N);
    sigma_sigma=[];
else
    if p<1
        % p is relative precision requirement of estimate
        sigma=sigmadec(N);
        sigma_sigma=ones(size(sigma));
        k=1;
        sigma2=sigma.^2;
        while max(sigma_sigma./sigma)>p
            sigma_next=sigmadec(N);
            sigma=(sigma*k+sigma_next)/(k+1);
            sigma2=sigma2+sigma_next.^2;
            sigma_sigma=sqrt(sigma2/(k+1) - sigma.^2)/sqrt(k+1);
            disp(['round ',int2str(k+1),', the relative deviations:'])
            disp(sigma_sigma./sigma)
            k=k+1;
        end
    else
        % p is sample size to do the estimate
        sigma=zeros(N,p);
        parfor k=1:p
            sigma(:,k)=sigmadec(N);
        end
        sigma_sigma=std(sigma,0,2);
    end
end
return

function sigma=sigmadec(N)
%SIGMADEC sigma decompose
rng('shuffle')
x=normrnd(0,1,2^N);
[~,w]=sdwt2(x);
[~,~,J]=size(w);
sigma=zeros(J,1);
for k=1:J
    x=w(:,:,k);
    sigma(k)=std(x(:));
end
return
