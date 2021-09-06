function rect = square2rect(sq,sizeR)
sizeS = size(sq);
D = sizeS - sizeR;
C = ceil(D/2) + 1;
rect = imcrop(sq,[C(2) C(1) sizeR(2)-1 sizeR(1)-1]);
return
