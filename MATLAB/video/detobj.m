function [xoffset,yoffset,coefmax,wobj,wscene,coef]=detobj(obj,scene,J)
%DETECTOBJ detect input object in the given scene

[do1,do2]=size(obj);
[ds1,ds2]=size(scene);
if do1>ds1 || do2>ds2
    error('object size exceeds the scene.');
end
[~,wobj]=sdwt2(obj,J);
[~,wscene]=sdwt2(scene,J);
coef=imcorr(wobj(:,:,J),wscene(:,:,J));
[~,xoffset]=max(max(coef,[],1)); % x coordinate of offset
[~,yoffset]=max(max(coef,[],2)); % y coordinate of offset
coefmax=coef(yoffset,xoffset);  % maximum coef
return
