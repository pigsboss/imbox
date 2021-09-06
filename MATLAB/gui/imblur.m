function imblur(varargin)
%IMBLUR image blur tool
%
%Syntax:
%  imblur(I,filter_parameters)
%  imblur(filter_parameters)
%  I = imblur(...)
%
%Input arguments:
%  I - input image
%  filter_parameters - parameters will be passed to fspecial calling
%
%Returns:
%  I - blurred image

I=[];
if ischar(varargin{1})
    filter_offset=0;
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
elseif isnumeric(varargin{1})
    I=varargin{1};
    filter_offset=1;
end
if isempty(I)
    error('Blur target does not exist.')
end
switch length(varargin)-filter_offset
    case 1
        h=fspecial(varargin{filter_offset+1});
    case 2
        if ischar(varargin{filter_offset+2})
            param1=str2double(varargin{filter_offset+2});
        else
            param1=varargin{filter_offset+2};
        end
        h=fspecial(varargin{filter_offset+1},param1);
    case 3
        if ischar(varargin{filter_offset+2})
            param1=str2double(varargin{filter_offset+2});
        else
            param1=varargin{filter_offset+2};
        end
        if ischar(varargin{filter_offset+3})
            param2=str2double(varargin{filter_offset+3});
        else
            param2=varargin{filter_offset+3};
        end
        h=fspecial(varargin{filter_offset+1},param1,param2);
    otherwise
        error('Unsupported calling of fspecial.')
end
figure;imagesc(imfilter(I,h,'circular'));axis image
return