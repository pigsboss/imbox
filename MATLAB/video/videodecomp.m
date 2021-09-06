function videodecomp(ifname,ofname,withdisplay)
%VIDEOCONVERT convert input video into uncompressed avi or mj2 video.
if nargin==2
    withdisplay=0;
end
extname=ofname((length(ofname)-3):length(ofname));
if ~strcmpi(extname,'.avi')
    ofname=[ofname,'.avi'];
end
ivid=VideoReader(ifname);
ovid=VideoWriter(ofname,'Uncompressed AVI');
open(ovid);
for k=1:ivid.NumberOfFrames
    frame=read(ivid,k);
    writeVideo(ovid,frame);
    if withdisplay
        imshow(frame);
        drawnow
    end
end
close(ovid);
return