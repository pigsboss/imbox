function d = wedist(p1,img,wt)
%WEDIST Weighted Euclidean distance between two patches
sizeI = size(img);
sizeP = size(p1);
d = zeros(sizeI);
rgv = (1:sizeP(1)) - (sizeP(1)+1)/2;
cgv = (1:sizeP(2)) - (sizeP(2)+1)/2;
for m = 1:sizeI(1)
  for n = 1:sizeI(2)
    p2 = subim(img,m+rgv,n+cgv,'symmetric');
    d(m,n) = sqrt(nbhdist(p1,p2,wt));
  end
end
return
