function [Amtv,xmtv,ymtv,A,TV] = findposmtv(I,P,xrng,yrng)
%FINDPOSMTV Find position of an object with given profile through total-variance
%minimization.
%INPUT
%  I      is the image contains the object as well as a background galaxy.
%  P      is the profile of the object.
%  xrng   is the range of x coordinate of the searching area.
%  yrng   is the range of y coordinate of the searching area.
%RETURN
%  Amtv   is the best estimate of the amplitude of the object.
%  xmtv   is the best estimate of the x coordinate of the object.
%  ymtv   is the best estimate of the y coordinate of the object.
%  A      is the amplitude map on the searching area.
%  TV     is the total variance map on the searching area.
[rngX,rngY] = meshgrid(xrng,yrng);
A = zeros(size(rngX));
TV = A;
disp(['Searching from (',num2str(min(xrng)),',',...
  num2str(min(yrng)),') to (',...
  num2str(max(xrng)),',',num2str(max(yrng)),')...'])
for k = 1:numel(A);
  [Apos,TVpos] = findampmtv(I,P,rngX(k),rngY(k));
  A(k) = Apos(length(Apos));
  TV(k) = TVpos(length(Apos));
end
[~,idx] = min(TV(:));
xmtv = rngX(idx);
ymtv = rngY(idx);
Amtv = A(idx);
disp(['Best estimate of amplitude: ',num2str(Amtv)])
disp(['Best estimate of x coordinate: ',num2str(xmtv)])
disp(['Best estimate of y coordinate: ',num2str(ymtv)])
return
