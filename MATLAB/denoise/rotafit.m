function [Fxy,Fp,Res] = rotafit(fxy,order)
%ROTAFIT Rotation-symmetric fitting of 2-D data.
%INPUTS
% fxy is the sampled 2-D data.
%RETURN
% Fxy is the rotation-symmetric model that fit the data.
% Fp  is the model in polar coordinate grid.

if nargin==1
  order = -1;
end
% find extrapval
[dy,dx] = size(fxy);
twinx = tukeywin(dx,0.1);
twiny = tukeywin(dy,0.1);
twin = twiny*(twinx)';
twin = min(max(1-twin,0),1);
clear twinx twiny
extrapval = sum(fxy(:).*twin(:))/sum(twin(:));
% disp(['extrapolation value: ',num2str(extrapval)])
clear twin

% elliptical fitting
[fp,rho,theta,~,~] = xy2polar(fxy,[],[],'nearest',extrapval);
[dphi,drho] = size(fp);
nrhoidx = zeros(dphi,1);
for k = 1:dphi
    nrhoidx(k) = find(fp(k,:)<=fp(k,1)*0.5,1);
%     nrhoidx(k) = round(sum(abs(fp(k,:)).*(1:drho)) / sum(abs(fp(k,:))));
end
[b,ecc,phi,~,eof] = ellifitbf(rho(1,nrhoidx)',theta(:,1));
sl = zeros(1,drho);
if order>=0
  nrhoidx = round(sdwt2(nrhoidx,order));
  for k = 1:dphi
  %     sl = sl+interp1(rho(k,:),fp(k,:),...
  %         rho(k,:)/sqrt(1-ecc^2*cos(theta(k,1)-phi)^2),...
  %         'nearest','extrap')/eof(k);
      sl = sl+interp1(rho(k,:),fp(k,:),...
          rho(k,:)/sqrt(1-ecc^2*cos(theta(k,1)-phi)^2),...
          'nearest','extrap');
  end
  sl = sl/dphi;
  % sl = sl/sum(1./eof);
  Fp = zeros(dphi,drho);
  for k = 1:dphi
      Fp(k,:) = interp1(rho(k,:),sl,...
          rho(k,:)/rho(1,nrhoidx(k))*b,'nearest','extrap');
  end
else
  for k = 1:dphi
%     sl = sl+interp1(rho(k,:),fp(k,:),...
%         rho(k,:)/sqrt(1-ecc^2*cos(theta(k,1)-phi)^2),...
%         'nearest','extrap')/eof(k);
    sl = sl+interp1(rho(k,:),fp(k,:),...
        rho(k,:)/sqrt(1-ecc^2*cos(theta(k,1)-phi)^2),...
        'nearest','extrap');
  end
  sl = sl/dphi;
  % sl = sl/sum(1./eof);
  Fp = zeros(dphi,drho);
  for k = 1:dphi
      Fp(k,:) = interp1(rho(k,:),sl,...
          rho(k,:)*sqrt(1-ecc^2*cos(theta(k,1)-phi)^2),'nearest','extrap');
  end
end
Fxy = polar2xy(Fp);
Res = fxy - Fxy;
return
