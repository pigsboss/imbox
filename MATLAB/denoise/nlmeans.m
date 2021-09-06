function varargout = nlmeans(img,pradius,alpha,h)
%NLMEANS Non-local Means Method
%
%References:
% [1] Buades, et al., On image denoising methods, Int J Comput Vis, 76: 123 -
%     139, 2008.
% [2] Buades and Morel, A non-local algorithm for image denoising, in IEEE
%     international conference on computer vision and pattern recognition, 2005.
sizeI = size(img);
psize = [2*pradius+1 2*pradius+1];
disp(['size of image: ',int2str(sizeI(1)),' x ',int2str(sizeI(2))])
disp(['size of patch: ',int2str(psize(1)),' x ',int2str(psize(2))])
wt = fspecial('gauss',psize,alpha);
disp(['standard deviation of Gaussian weighted Euclidean distance: ',num2str(alpha)])
disp(['standard deviation of Gaussian weighted similarity: ',num2str(h)])
denimg = img;
NumPix = numel(img);
h2 = h^2;
N1 = sizeI(1);
parfor  n1 = 1:NumPix
  r1 = mod(n1,N1);
  c1 = (n1-r1)/N1 + 1;
  p1 = pxnbh(img,[r1,c1],pradius);
  s = exp(-1*awedist(p1,img,wt).^2 / h2);
  denimg(n1) = sum(img(:).*s(:))/sum(s(:));
end
if nargout == 0
  figure('Name','Input image');
  imagesc(img);colormap('jet');axis image;drawnow
  figure('Name','Denoised image');
  imagesc(denimg);colormap('jet');axis image;drawnow
else
  varargout{1} = denimg;
end
return
