function videoconvert(ifname,ofname,type)
%VIDEOCONVERT convert input video into uncompressed avi or mj2 video.
extname=ofname((length(ofname)-3):length(ofname));
defaulttype='avi';
if nargin==2
    if strcmpi(extname,'.avi')
        type='avi';
    elseif strcmpi(extname,'.mj2')
        type='mj2';
    else
        type=defaulttype;
    end
end
ivid=VideoReader(ifname);
if strcmpi(type,'avi')
    ovid=VideoWriter(ofname,'Uncompressed AVI');
elseif strcmpi(type,'mj2')
    ovid=VideoWriter(ofname,'Archival');
else
    error('Unsupported type.')
end
open(ovid);
for k=1:ivid.NumberOfFrames
    frame=read(ivid,k);
    writeVideo(ovid,frame);
end
close(ovid);
return