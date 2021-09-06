function [A,TV] = findampmtv(I,P,x,y)
%FINDAMPMTV Find amplitude of an object with given profile at given position through
%total-variance-minimization.
%INPUTS
%  I        is the image contains a point source and a smooth background object.
%  P        is the profile of the object.
%  x and y  are coordinates of the point source.
%RETURN
%  A is the best estimate of the amplitude of the point source.
%  TV is the entropy (Skilling and Gull) of the remenant background object.
epsilon = 1e-3;
ML = 20;
P = subshift(P,[y,x]);
A = zeros(ML,1);
A(1) = 0;
TV = A;
TV(1) = totalvar(I-A(1)*P);
G = I - A(1)*P;
for l = 2:ML
  f0 = -sum(sum(dlt2(G).*dlt2(P)));
  f1 =  sum(sum(dlt2(P).^2));
  A(l) = A(l-1) - f0/f1;
  G = I - A(l)*P;
  TV(l) = totalvar(G);
  if abs(A(l) - A(l-1)) <= epsilon
    break
  end
end
A = A(1:l);
TV = TV(1:l);
return