function g=imcorr(f,h)
%IMCORR Cross correlation for image processing.
% Dimension of f must not exceed dimension of h.
% f=f-mean(f(:));
% h=h-mean(h(:));

[df1,df2]=size(f);
[dh1,dh2]=size(h);
if df1>dh1 || df2>dh2
    error('Dimension of f must not exceed dimension of h.');
end
F=zeros(size(h));
F(1:df1,1:df2)=f;
g=real(ifft2(fft2((h)).*conj(fft2(F))));
dg1=dh1-df1+1;
dg2=dh2-df2+1;
g=g(1:dg1,1:dg2);
return