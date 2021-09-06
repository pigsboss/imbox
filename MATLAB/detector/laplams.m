function L = laplams(f,sc)
%LAPLAMS Multi-scale Laplacians of given image.
[H,W] = size(f);
NS = length(sc);
L = zeros(H,W,NS);
for k = 1:NS
    L(:,:,k) = sc(k)*dlt2(gsmooth(f,sc(k)));
end
return
