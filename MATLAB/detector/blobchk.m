function [bmes,bimg] = blobchk(bmes,bimg)
%BLOBCHK

DUP_THRES_PIX = 1.5; %threshold for duplicate test, in number of pixels.

% check duplicate blobs
NB = numel(bmes.flux);
for k = 1:NB
  % distance between sub-pixel centres:
  dctr = sqrt((bmes.xctrsubpix - bmes.xctrsubpix(k)).^2 + ...
    (bmes.yctrsubpix - bmes.yctrsubpix(k)).^2);
  % duplicate indices:
  didx = find(dctr<DUP_THRES_PIX);
  if numel(didx)>1
    % found duplication
    disp('Found duplication.')
    [~, m] = max(bmes.flux(didx));
    ridx = true(size(didx));
    ridx(m) = false;
    ridx = didx(ridx); % indices of blobs to remove
    bmes.flux(ridx) = -1; % set flux of discarded duplicates to -1.
  end
end
% remove duplicates:
idx = logical(bmes.flux>0);
bmes.xctr = bmes.xctr(idx);
bmes.yctr = bmes.yctr(idx);
bmes.scale = bmes.scale(idx);
bmes.rect = bmes.rect(idx,:);
bmes.xctrsubpix = bmes.xctrsubpix(idx);
bmes.yctrsubpix = bmes.yctrsubpix(idx);
bmes.flux = bmes.flux(idx);
bmes.peak = bmes.peak(idx);
bmes.eccentricity = bmes.eccentricity(idx);
bmes.semiminor = bmes.semiminor(idx);
bmes.semimajor = bmes.semimajor(idx);
bmes.posangle = bmes.posangle(idx);
bmes.elliparam = bmes.elliparam(idx,:);
bimg = bimg(idx);
% NB = numel(bmes.flux);

% check blended blobs

return