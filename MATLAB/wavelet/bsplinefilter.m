function h=bsplinefilter(sz)
%BSPLINEFILTER generate multi-scale b-spline filters with NxN sizes
f=zeros(sz);
f(1,1)=1;
f=fftshift(f);
[c,w]=dwt2(f);
[~,~,J]=size(w);
h=zeros(sz(1),sz(2),J+1);
h(:,:,J+1)=c;
for k = 2:J
    h(:,:,k)=c+sum(w(:,:,k:J),3);
end
h(1,1,1) = 1;
h(:,:,1) = fftshift(h(:,:,1));
return
