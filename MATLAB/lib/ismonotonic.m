function t=ismonotonic(A,threshold)
%ISMONOTONIC monotonic test function for matrix or vector A.
if nargin==1
  threshold=1e-14;
end
sz=size(A);
t=true;
for d=1:length(sz)
  if sz(d)>1
    x=diff(A,1,d);
    if std(x(:))>threshold
      t=false;
    end
  end
end
return
