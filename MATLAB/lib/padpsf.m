function PSF = padpsf(PSF,sizeI)
[~,~,lenPSF] = size(PSF);
sizePSF = size(PSF);
if any(sizeI ~= sizePSF)
  PSFext = zeros(sizeI);
  for k = 1:lenPSF
    sizePSF = size(PSF);
    PSF(:,:,k) = ifftshift(PSF(:,:,k));
    PSFtmp = [PSF(1:round(sizePSF(1)/2),:,k);...
        zeros(sizeI(1)-sizePSF(1),sizePSF(2));...
        PSF((1+round(sizePSF(1)/2)):sizePSF(1),:,k)];
    sizePSF = size(PSFtmp);
    PSFext(:,:,k) = [PSFtmp(:,1:round(sizePSF(2)/2)),...
        zeros(sizePSF(1),sizeI(2)-sizePSF(2)),...
        PSFtmp(:,(1+round(sizePSF(2)/2)):sizePSF(2))];
    PSFext(:,:,k) = fftshift(PSFext(:,:,k));
  end
  PSF = PSFext;
end
return
