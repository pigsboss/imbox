function psf = estpsf(I,btab,bimg)
%ESTPSF Estimate PSF from field stars in given image.
%
%INPUTS:
% I             is the input image, which can be scaled to display properly.
% btab and bimg are blob table and blob images returned by blobex.m. To estimate
%               PSF precisely, the input image of blobex.m cannot be scaled.
if ~isempty(I)
  fh = figure('Name','Estimate PSF');
  imagesc(I);axis image
else
  fh = gcf;
end
figure(fh);
ellidraw(btab(2:7,:),'color','r');
waituserinput = true;
[~,numblobs] = size(btab);
bselected = false(numblobs,1);
while waituserinput == true
  disp('Click inside the ellipse one by one to toggle field stars (red = unselected, green = selected).')
  disp('Press ESC to finish.')
  figure(fh);
  [x,y,button] = ginput(1);
  disp(['x = ',num2str(x)]);
  disp(['y = ',num2str(y)]);
  % find closest blob:
  d = sqrt((btab(2,:) - x).^2 + (btab(3,:) - y).^2);
  [dmin,idx] = min(d);
  disp(['closest blob: ',int2str(idx)])
  if dmin <= min(btab(4:5,idx))
    disp(['blob ',int2str(idx),' is toggled.'])
    bselected(idx) = ~(bselected(idx));
  else
    disp(['click outside blob ',int2str(idx),'.'])
    waituserinput = false;
  end
  figure(fh);
  if any(~(bselected)) == true
    ellidraw(btab(2:7,~(bselected)),'color','r')
  end
  if any(bselected) == true
    ellidraw(btab(2:7,bselected),'color','g')
  end
  if button == 27
    break
  end
end
psf = bloboverlay(btab(:,bselected),bimg(:,:,bselected));
psf = psf/sum(psf(:));
return

function bprofile = bloboverlay(btab,bimg)
[H,W,L] = size(bimg);
[X,Y] = meshgrid(1:W,1:H);
for l = 1:L
  tf = sum(sum(bimg(:,:,l)));
  xc = sum(sum(X.*bimg(:,:,l)))/tf;
  yc = sum(sum(Y.*bimg(:,:,l)))/tf;
  bimg(:,:,l) = subshift(bimg(:,:,l),[(H+1)/2-yc,(W+1)/2-xc]);
  bimg(:,:,l) = sqrt(btab(1,l))*bimg(:,:,l)/sum(sum(bimg(:,:,l)));
end
bprofile = mean(bimg,3);
return