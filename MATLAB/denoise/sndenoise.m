close all
clear
clc
rect = [400 400 511 511];
%%
xoffset1to2 = -5;
yoffset1to2 = -52;
%% load fits file
sn1 = fitsread('12ij/sn2012ij_B_001_130109_c.fit');
sn2 = fitsread('12ij/sn2012ij_B_002_130218_c.fit');
%% estimate psf from fits file 1
[btab,bimg] = blobex(sn1,20,170);
psf1 = estpsf(sqrt(max(min(sn1,500),0)),btab,bimg);
sizeP = size(psf1);
if numel(psf1)<64^2
  psf1 = padarray(psf1,[32-round(sizeP(1)/2),32-round(sizeP(2)/2)],0,'pre');
  sizeP = size(psf1);
  psf1 = padarray(psf1,[64-sizeP(1),64-sizeP(2)],0,'post');
end
if numel(psf1)>64^2
  sizeP = size(psf1);
  psf1 = imcrop(psf1,[sizeP(2)/2 sizeP(1)/2 63 63]);
end
%% estimate psf from fits file 2
[btab,bimg] = blobex(sn2,20,170);
psf2 = estpsf(sqrt(max(min(sn2,500),0)),btab,bimg);
sizeP = size(psf2);
if numel(psf2)<64^2
  psf2 = padarray(psf2,[32-round(sizeP(1)/2),32-round(sizeP(2)/2)],0,'pre');
  sizeP = size(psf1);
  psf2 = padarray(psf2,[64-sizeP(1),64-sizeP(2)],0,'post');
end
if numel(psf2)>64^2
  sizeP = size(psf2);
  psf2 = imcrop(psf2,[sizeP(2)/2 sizeP(1)/2 63 63]);
end
%% load denoised images
dsn1 = inlmerge('sn1');
dsn1 = imcrop(dsn1,[250+xoffset1to2, 300+yoffset1to2, 63, 63]);
dsn2 = inlmerge('sn2');
dsn2 = imcrop(dsn2,[250 300 63 63]);
%% find position and amplitude of SN in image 1
[A1,x1,y1,~,tvmap1] = findposmtv(sdwt2(dsn1,1),sdwt2(psf1,1),(-1:0.25:1)-9,(-1:0.25:1)-12);
%% find position and amplitude of SN in image 2
[A2,x2,y2,~,tvmap2] = findposmtv(sdwt2(dsn2,1),sdwt2(psf1,1),(-1:0.25:1)-9,(-1:0.25:1)-12);
