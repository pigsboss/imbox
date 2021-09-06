function [J,lnapd,mse]=deconvmap(I,PSF,BG,NUMIT,analyse_mode)
%DECONVMAP Deconvolution by maximum a posteriori estimate.
if nargin<5
    analyse_mode=0;
end
lnapd=zeros(NUMIT,1);
mse=zeros(NUMIT,1);
J=I;
IPSF=rot90(PSF,2);
H=psf2otf(PSF);
IH=psf2otf(IPSF);
%sigma=zeros(NUMIT,1);
%tic
for k=1:NUMIT
    Jp=J;
    J=J.*exp(real(ifft2(fft2(I./real(ifft2(H.*fft2(J)))-1).*IH)));
    J=J.*double(J>=BG)+BG.*double(J<BG);
    if analyse_mode~=0
        lnapd(k)=logpostpdfpois(I,PSF,J);
        mse(k)=sum((Jp-J).^2)/numel(J);
    end
%    sigma(k)=std(I(:)-reshape(imconv(J,PSF),numel(I),[]));
end
%toc
return