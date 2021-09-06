function [xoffset,yoffset,cmax] = adetobj(obj,scene,J,M)
%ADETOBJ Accelerated object detection
% M is number multiple-resolution detections.
if nargin <3
    J = 1;
    M = 1;
end
if nargin <4
    M = 1;
end
[do1,do2]=size(obj);
[ds1,ds2]=size(scene);
if do1>ds1 || do2>ds2
    error('object size exceeds the scene.');
end
J = min(J,lastpow2(min(do1,do2)));
M = min(M,lastpow2(min(do1,do2)));
xoffset = 0;
yoffset = 0;
for m=1:M
    resz = 2^(m-M);
    if m == 1
        y0 = 1;
        x0 = 1;
        part = scene;
    else
        y0 = max(1,round(yoffset-2^(1+M-m)));
        x0 = max(1,round(xoffset-2^(1+M-m)));
        part = scene(y0:min(ds1,round(yoffset+do1+2/resz)),...
            x0:min(ds2,round(xoffset+do2+2/resz)));
    end
    if resz < 1
        [x,y,cmax] = detobj(imresize(obj,resz,'nearest'),...
            imresize(part,resz,'nearest'),...
            max(1,J-(M-m)));
    else
        [x,y,cmax] = detobj(obj,part,J);
    end
    xoffset = x0+floor((x-1)/resz);
    yoffset = y0+floor((y-1)/resz);
end
return
