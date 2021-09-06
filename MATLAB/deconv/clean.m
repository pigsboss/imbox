function [C,R] = clean(I,PSF,beta,MAXNUMIT)
%CLEAN A simple implementation of the CLEAN algorithm.
if nargin<4
  MAXNUMIT=1e4;
end
R = I;
% sigma = sigmaclipping(R);
% k = 2;
C = zeros(size(I));
PSF = ifftshift(PSF);
PSF = PSF/max(PSF(:));
sumPSF = sum(PSF(:));
maxR = max(R(:));
meanR = mean(R(:));
l = 1;
sk = skewness(R(:));
while (sk>=0) && (maxR > 0) && (maxR > meanR) && (l<=MAXNUMIT)
    sk = skewness(R(:));
    disp(skewness(R(:)))
    [~,x] = max(max(R,[],1),[],2);
    [maxR,y] = max(max(R,[],2),[],1);
%     sigma = std(R(:));
    meanR = mean(R(:));
    R = R-beta*(maxR-meanR)*circshift(PSF,[y-1,x-1]);
    C(y,x) = C(y,x)+beta*(maxR-meanR)*sumPSF;
    l = l+1;
end
return
