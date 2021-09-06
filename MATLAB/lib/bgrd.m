function LL = bgrd(I,PSF,METHOD)
%BGRD Use specified method to obtain background of specified image.
%
%References:
%   [1] http://en.wikipedia.org/wiki/Full_width_at_half_maximum
%

switch upper(METHOD)
  case {'MANUAL','MAN','M'}
    figure;imagesc(I);drawnow;
    disp('Please select a region to estimate the background.')
    [~,LL] = imflux;
    disp(['The background level estimate is ',num2str(LL)])
  case {'INTRINSIC','INT','I'}
    IPSF = rot90(PSF,2);
    LL = imconv(I,IPSF,'symmetric');
  case {'SPLINE','SPL','S'}
    LL = sdwt2(I,ceil(log2(fwhm(PSF))-1));
  case {'GAUSSIAN','GAUSS','GAU','G'}
    LL = gsmooth(I,fwhm(PSF)/2.35482);
  case {'CLEAN','CLN','C'}
    IPSF = rot90(PSF,2);
    [~,LL] = clean(imconv(I,IPSF,'symmetric'),imconv(PSF,IPSF,'symmetric'),0.2);
  otherwise
    error('Unsupported parameter BACKGROUND. Please run "help dd" to refer to the document.')
end
return


function W = fwhm(X)
[x,y] = find(X>0.5*max(X(:)));
W = max(max(x)-min(x),max(y)-min(y));
return