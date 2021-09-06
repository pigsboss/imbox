function varargout=imsect(varargin)
%IMSECT show selected cross-section of image
%
%Syntax:
%  imcsect
%  imcsect(I)
%  imcsect(NP)
%  imcsect(METHOD)
%  imcsect(I, NP)
%  imcsect(I, METHOD)
%  imcsect(NP, METHOD)
%  imcsect(METHOD, NP)
%  imcsect(I, NP, METHOD)
%  sect = imcsect(...)
%
%Available METHOD:
%  'nearest' Nearest neighbor interpolation
%  'linear' Linear interpolation (default)
%  'spline' Cubic spline interpolation
%  'cubic' Cubic interpolation, as long as data is uniformly-spaced. 
%   interp2 switches to 'spline' and issues a warning if the data is not equally spaced.
%
%Input arguments:
%  I - 2D image
%
%Returns:
%  sect - 1D cross-section data (interpolated)
%  XI - x coordinates of sect
%  YI - y coordinates of sect
I=[];
NP=[];
METHOD=[];
hf=[];
ha=[];
canvas=[];
switch nargin
    case 0
    case 1
        if ischar(varargin{1})
            if ~isnan(str2double(varargin{1}))
                NP=round(str2double(varargin{1}));
            else
                METHOD=varargin{1};
            end
        else
            if isscalar(varargin{1})
                NP=round(varargin{1});
            else
                I=varargin{1};
            end
        end
    case 2
        if ischar(varargin{1})
            if ~isnan(str2double(varargin{1}))
                NP=round(str2double(varargin{1}));
                METHOD=varargin{2};
            else
                METHOD=varargin{1};
                NP=round(str2double(varargin{2}));
            end
        else
            if isscalar(varargin{1})
                NP=round(varargin{1});
                METHOD=varargin{2};
            else
                I=varargin{1};
                if ischar(varargin{2})
                    METHOD=varargin{2};
                else
                    NP=varargin{2};
                end
            end
        end
    case 3
        I=varargin{1};
        NP=round(varargin{2});
        METHOD=varargin{3};
    otherwise
        error('Syntax error.')
end
if isempty(I)
    hf=get(0,'CurrentFigure');
    if ~isempty(hf)
        ha=get(hf,'CurrentAxes');
        if ~isempty(ha)
            canvas=findobj(ha,'Type','image');
            if ~isempty(canvas)
                I=get(canvas,'CData');
            end
        end
    end
end
if isempty(I)
    error('No image data specified.')
end
sz=size(I);
if length(sz)~=2
    error('Specified image data must be two-dimensional.')
end
if isempty(canvas)
    hf=figure;ha=axes;imagesc(I);drawnow
%     canvas=get(ha,'Children');
end
if isempty(NP)
    NP=100;
end
if isempty(METHOD)
    METHOD='linear';
end
hline=findobj(ha,'Tag','imline');
if isempty(hline)
    disp('Please click on the image twice to input the end points of the cross section.')
    figure(hf);axes(ha);[xpt,ypt]=ginput(2);
    imline(ha,xpt,ypt);
    drawnow
else
    xpt=zeros(1,2);
    ypt=xpt;
    xpt(1)=get(findobj(hline,'Tag','end point 1'),'XData');
    xpt(2)=get(findobj(hline,'Tag','end point 2'),'XData');
    ypt(1)=get(findobj(hline,'Tag','end point 1'),'YData');
    ypt(2)=get(findobj(hline,'Tag','end point 2'),'YData');
end
xstep=(xpt(2)-xpt(1))/NP;
ystep=(ypt(2)-ypt(1))/NP;
XI=xpt(1):xstep:xpt(2);
YI=ypt(1):ystep:ypt(2);
[X,Y]=meshgrid(1:sz(2),1:sz(1));
sect = interp2(X,Y,I,XI,YI,METHOD);
figure;plot(XI,sect);axis tight;xlabel('X');ylabel('CData')
switch nargout
    case 1
        varargout{1}=[sect;XI;YI];
    case 2
        varargout{1}=sect;
        varargout{2}=[XI;YI];
    case 3
        varargout{1}=sect;
        varargout{2}=XI;
        varargout{3}=YI;
end
return
