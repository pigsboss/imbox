function t = fftthreshold(sz)
nElem = prod(sz);
nOps = 0;
for k = 1:length(sz)
    nffts = nElem/sz(k);
    nOps  = nOps + sz(k)*log2(sz(k))*nffts;
end
t = eps*nOps;
return