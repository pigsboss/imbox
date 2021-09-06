function [varargout]=imflux(varargin)
%IMFLUX calculate image flux
%Syntax:
% imflux
% imflux=(x)
% imflux=(x,rect)
% [TF,MI,IS,NUMPIX,MIS]=imflux
% [TF,MI,IS,NUMPIX,MIS]=imflux(x)
% [TF,MI,IS,NUMPIX,MIS]=imflux(x,rect)
%
%Input arguments:
% x (optional)    - image array or image filename
% rect (optional) - four-element position parameter defines a region in x.
%   format: [x_topleft y_topleft width height]
%
%Output:
% TF     - total flux
% MI     - mean intensity
% IS     - standard deviation of intensity
% MIS    - standard deviation of mean intensity
% NUMPIX - number of pixels
rect=0;
switch nargin
    case 0
        x=imcrop;
    case 1
        x=varargin{1};
    case 2
        x=varargin{1};
        rect=varargin{2};
    otherwise
        error('Syntax error')
end
if ischar(x)
    x=double(imread(x));
end
x=mean(x,3);
if isempty(rect)
    h=imagesc(x);axis image
    x=imcrop(h);
elseif ~isscalar(rect)
    x=x(rect(2):(rect(2)+rect(4)), rect(1):(rect(1)+rect(3)));
end
TF=sum(x(:));
MI=mean(x(:));
NUMPIX=numel(x);
MIS=std(x(:))/sqrt(NUMPIX);
IS=std(x(:));
switch nargout
    case 1
        varargout{1}=TF;
    case 2
        varargout{1}=TF;
        varargout{2}=MI;
    case 3
        varargout{1}=TF;
        varargout{2}=MI;
        varargout{3}=IS;
    case 4
        varargout{1}=TF;
        varargout{2}=MI;
        varargout{3}=IS;
        varargout{4}=NUMPIX;
    case 5
        varargout{1}=TF;
        varargout{2}=MI;
        varargout{3}=IS;
        varargout{4}=NUMPIX;
        varargout{5}=MIS;
    otherwise
        disp(['Number of pixels: ',int2str(NUMPIX)])
        disp(['Total flux: ',num2str(TF)])
        disp(['MAX intensity: ',num2str(max(x(:)))])
        disp(['MIN intensity: ',num2str(min(x(:)))])
        disp(['Average intensity: ',num2str(MI)])
        disp(['Standard deviation of intensity: ', num2str(IS)])
        disp(['Error of average intensity: ',num2str(MIS)])
end
return
