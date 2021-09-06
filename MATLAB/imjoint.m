clear
close all

%% load image
im1 = 'IMG.bmp';
im2 = 'IMG_0001.bmp';
NA = 20;
im1data = double(imread(im1))/255;
im2data = double(imread(im2))/255;
[H1,W1,C1] = size(im1data);
[H2,W2,C2] = size(im2data);
im1mask = ones(H1,W1);
im2mask = ones(H2,W2);
if C1~=C2
  error('The two image must be of the same type.')
end

%% set anchors in image 1
aimg = cell(NA,1);
apos = cell(NA,1);
xa1 = zeros(NA,1);
ya1 = zeros(NA,1);
for n = 1:NA
  disp(['Select anchor ',int2str(n),' of ',int2str(NA),':'])
  figure('name','Select an anchor');[aimg{n},apos{n}] = imcrop(im1data);
  figure('name',['Selected anchor ',int2str(n)]);imshow(aimg{n})
  xa1(n) = round(apos{n}(1));
  ya1(n) = round(apos{n}(2));
end

%% detect shifts of anchors in image 2
x = zeros(NA,C1);
y = zeros(NA,C1);
cmax = zeros(NA,C1);
xa2 = zeros(NA,1);
ya2 = zeros(NA,1);
wta = zeros(NA,1);
for n=1:NA
  for c=1:C1
    [x(n,c),y(n,c),cmax(n,c)]=adetobj(double(1-aimg{n}(:,:,c)),double(1-im2data(:,:,c)),2,4);
  end
  [~,idx] = max(cmax(n,c));
  xa2(n) = x(n,idx);
  ya2(n) = y(n,idx);
  wta(n) = 1/(std(x(n,:)).^2 + std(y(n,:)).^2 + 1e-6);
  disp(['Goodness of detection: ',num2str(wta(n))]);
end

%% calculate shift and rotation angle
x01 = sum(xa1.*wta)/sum(wta);
y01 = sum(ya1.*wta)/sum(wta);
rho1 = sqrt((xa1-x01).^2 + (ya1-y01).^2);
phi1 = angle((xa1-x01)+1i*(ya1-y01));
x02 = sum(xa2.*wta)/sum(wta);
y02 = sum(ya2.*wta)/sum(wta);
rho2 = sqrt((xa2-x02).^2 + (ya2-y02).^2);
phi2 = angle((xa2-x02)+1i*(ya2-y02));
phi = sum((phi1-phi2).*(rho1+rho2).*wta)/sum((rho1+rho2).*wta);

%% rotate, shift and joint
im3data = zeros(H1+H2,W1+W2,C1);
im1data = double(padarray(im1data,[H2,W2],0,'post'));
im2data = double(padarray(im2data,[H1,W1],0,'post'));
im1mask = double(padarray(im1mask,[H2,W2],0,'post'));
im2mask = double(padarray(im2mask,[H1,W1],0,'post'));
x0 = ceil(0.5*(W1+W2+1));
y0 = ceil(0.5*(H1+H2+1));
im1mask = subshift(im1mask,[y0-y01,x0-x01]);
im2mask = imrotate(subshift(im2mask,[y0-y02,x0-x02]),-phi*180/pi,...
  'nearest','crop');
xwin = 1:(W1+W2);
xwin = gsmooth(double(xwin <= 0.5*(W1+W2)),1);
win = min(1,max(0,ones(H1+H2,1)*xwin));
im1mask = im1mask.*win;
im2mask = im2mask.*(1-win);
im3mask = im1mask + im2mask;
for c=1:3
  im1data(:,:,c) = subshift(im1data(:,:,c),...
    [y0-y01,x0-x01]);
  im2data(:,:,c) = imrotate(subshift(im2data(:,:,c),...
    [y0-y02,x0-x02]),-phi*180/pi,'bilinear','crop');
  d = im1data(:,:,c).*im1mask+im2data(:,:,c).*im2mask;
  d(im3mask>=1) = d(im3mask>=1)./im3mask(im3mask>=1);
  im3data(:,:,c) = d;
end

%% draw joint image
figure;imshow(im3data)

%% finish
im3 = input('output filename:','s');figure;imwrite(imcrop(im3data),im3,'png')
