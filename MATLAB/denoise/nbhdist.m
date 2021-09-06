function d = nbhdist(p1,p2,wt)
%NBHDIST Observation-independent distance between two neighborhoods.
%INPUTS:
% p1 and p2 are the two patches (neighborhoods).
% wt        is the Gaussian weighing filter.
%
%RETURN:
% d         is the observation-independent weighted Euclidean distance between the
%           two patches.
d = zeros(4,1);
d(1) = sum(((p1(:)-p2(:)).^2).*wt(:));
for n=1:3
  p2 = rot90(p2);
  d(n+1) = sum(((p1(:)-p2(:)).^2).*wt(:));
end
d = min(d);
return
