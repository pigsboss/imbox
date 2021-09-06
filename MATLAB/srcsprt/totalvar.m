function V = totalvar(I)
%TOTALVAR Total variance of input image I.
V = (dlt2(I)).^2;
V = sum(V(:));
return
