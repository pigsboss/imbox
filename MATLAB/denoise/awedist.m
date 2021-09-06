function d = awedist(p1,img,wt)
%AWEDIST Accelerated weighted Euclidean distance between two patches
p1 = rot90(p1,2);
wt = rot90(wt,2);
% img = rot90(img,2);
d = sqrt(abs(sum((p1(:).^2).*wt(:))+...
  imconv(img.^2,wt,'circular')-...
  2*imconv(img,wt.*p1,'circular')));
return
