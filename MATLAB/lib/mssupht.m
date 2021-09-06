function M=mssupht(w,sigma,h)
%MSSUPHT multi-resolution support function, with 3-sigma hard threshold
% w - multiresolution wavelet data
% sigma - multiresolution noise levels of the wavelet data
%
[H,W,J]=size(w);
if length(sigma)~=J
    error('Maximum scales of the wavelet data and its noise level do not match.')
end
M=zeros(size(w));
K=sqrt(2*log10(H*W));
% g=fspecial('gauss',[H W],2);
for k=1:J
    M(:,:,k)=double(w(:,:,k) >= (K*sigma(k)));
    M(:,:,k)=imconv(M(:,:,k),h(:,:,k));
end
return
