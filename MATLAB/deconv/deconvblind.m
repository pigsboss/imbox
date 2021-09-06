function [J,PSF]=deconvblind(I,NUMIBD,NUMPSF)
%DECONVBLIND Blind deconvolution

MAX_PSF_SZ = [64 64]; % maximum psf size

path('./lib',path)
path('./lib/wavelet',path)

% check image size
sizeI = size(I);
if length(sizeI) ~= 2
    error('Right now this program supports only 2D image.');
end

% creat a window
winfun = tukeywin(sizeI(1))*(tukeywin(sizeI(2)))';

% wavelet transform
% I = I.*winfun;
[~,wI]=dwt2(I.*winfun, 1);

% get initial PSF
PSF = max(imcorr(wI(:,:,1), padarray(wI(:,:,1),MAX_PSF_SZ/2)),0);
sizePSF = size(PSF);
PSF = padarray(PSF, sizeI - sizePSF, 0, 'post');
PSF = circshift(PSF, floor(0.5*(sizeI - sizePSF)));
PSF = PSF/sum(PSF(:));
clear wI

cellI = {I/sum(I(:)), I/sum(I(:))};
cellPSF = {I/sum(I(:)), PSF};
for n = 1:NUMIBD % number of blind deconvolution iterations
    cellPSF = deconvml(cellPSF,cellI{2},NUMPSF,[],[],[],[],[],[],3);
%     cellPSF{2} = cellPSF{2}.*winfun;
    cellI = deconvml(cellI,cellPSF{2},1,[],[],[],[],[],[],3);
end
J = cellI{2};
PSF = cellPSF{2};
return