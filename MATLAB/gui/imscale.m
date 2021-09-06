function imscale(varargin)
%IMSCALE
%
%Syntax:
%   imscale(I, min:max, scale_method)
%   imscale min:max scale_method
%

[I,canvas,scale,LL,UL] = parse_inputs(varargin{:});
if ~isempty(LL), I = max(I,LL);end
if ~isempty(UL), I = min(I,UL);end
switch scale
    case 'linear'
    case 'log'
        if min(I(:)) > 0
            I = log10(I);
            cbarstr = 'log to base 10';
        else
            error('Logarithm scaling not allowed for this image.')
        end
    case 'sqrt'
        if min(I(:)) >= 0
            I = sqrt(I);
            cbarstr = 'square root';
        else
            error('Square-root scaling not allowed for this image.')
        end
    case 'square'
            I = (I).^2;
            cbarstr = 'square';
    case 'exp'
            I = 10.^(I);
            cbarstr = 'exponential to base 10';
end
if isempty(canvas)
    imagesc(I);axis image
else
    figure;imagesc(I);axis image;cbax = colorbar;drawnow
    set(get(cbax,'YLabel'),'string',cbarstr,...
    'fontname','Helvetica',...
    'fontsize',16)
end
return

function [I, canvas, scale, LL, UL] = parse_inputs(varargin)
I = [];
canvas=[];
LL = [];
UL = [];
scale = [];scale_d = 'linear';
narginchk(0,3);
for k = 1:nargin
    switch varargin{k}
        case {'linear','log','sqrt','square','exp'}
            scale = varargin{k};
        otherwise
            if ischar(varargin{k})
                [LL, UL] = str2limits(varargin{k});
            elseif isvector(varargin{k})
                LL = min(varargin{k});
                UL = max(varargin{k});
            end
    end
end
if isempty(I),      [I, canvas] = getcanvas;end
if isempty(I),      error('No image data specified.');end
if isempty(scale),  scale = scale_d;end
return

function [LL, UL] = str2limits(str)
LL = [];
UL = [];
dlm = strfind(str,':');
strlen = length(str);
if length(dlm) == 1
    LL = str2double(str(1:(dlm-1)));
    UL = str2double(str((dlm+1):strlen));
    if isnan(LL), LL = [];end
    if isnan(UL), UL = [];end
end
return