function [I, xgv, ygv, canvas] = getcanvas
%GETCANVAS get both CData (the image) and the handle of the current canvas object.
I = [];
xgv = [];
ygv = [];
canvas = [];
hf=get(0,'CurrentFigure');
if ~isempty(hf)
	ha=get(hf,'CurrentAxes');
	if ~isempty(ha)
		canvas=findobj(ha,'Type','image');
		if ~isempty(canvas)
			I=get(canvas,'CData');
      xgv = get(canvas,'XData');
      ygv = get(canvas,'YData');
      if length(xgv)==2
        xgv = xgv(1):xgv(2);
      end
      if length(ygv)==2
        ygv = ygv(1):ygv(2);
      end
      if numel(xgv)*numel(ygv)~=numel(I)
        sizeI = size(I);
        xgv=1:sizeI(2);
        ygv=1:sizeI(1);
      end
		end
	end
end
return
