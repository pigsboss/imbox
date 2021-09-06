function sz = kernelSize(h,threshold)
if nargin == 1
  threshold = 1e-8*max(h(:));
end
[y,x] = find(h>=threshold*max(h(:)));
sz = [max(y)-min(y)+1, max(x)-min(x)+1];
return
