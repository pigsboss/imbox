function Btab = deconvnull(PSF,hdeconv,NUMIT,SmoothScale,NUMMC)
%DECONVNULL compute histogram of null distribution of blobs deconvolved by
%specified program.
%
%Input arguments:
%   PSF           is the point spread function used to generate null observed
%                 data.
%   hdeconv       is the handle of specified deconvolution program.
%   NUMIT         is the number of iterations required for deconvolution.
%   SmoothScale   is the scale of smooth filter used before blob extraction.
%   NUMMC         is number of Monte-Carlo loops.
%
%Available deconvolutors:
%   deconvml
%   deconvmsme

Bgrd = 100;
sizePSF = size(PSF);
rng('shuffle')
Btabs = cell(1,NUMMC);
ImFiles = cell(1,NUMMC);
filename = [func2str(hdeconv),'_',int2str(NUMIT),'x',num2str(SmoothScale),...
    'x',int2str(NUMMC)];
matname = [filename,'.hist.mat'];
tarname = [filename,'.profile.tar'];
hfmont = figure('Name','Null deconvolution benchmark, all blobs');

hbar = waitbar(0,'Running benchmark...','Name','Null deconvolution benchmark');
t = tic;

for m = 1:NUMMC
    g = poissrnd(imconv(ones(sizePSF(1:2))*Bgrd,PSF));
    f = hdeconv(g,PSF,NUMIT,mean(g(:)));
    [Btabs{m},Bimg] = blobex(f,SmoothScale);
    [H,W,~] = size(Bimg);
    figure(hfmont)
    hmont = montage(reshape(Bimg,H,W,1,[])*255,hot(256));drawnow
    ImFiles{m} = [filename,'_mon_',int2str(m),'.png'];
    imwrite(uint8(get(hmont,'CData')),hot(256),...
        ImFiles{m},'png')
    tElapsed = toc(t);
    tRemain = ((NUMMC-m)/m)*tElapsed;
    waitbar(m/NUMMC,hbar,['Running benchmark, ',...
        num2str(m/NUMMC*100),'% completed, ',...
        int2str(round(tRemain)),' seconds left...']);
end
if exist('hbar','var')==1,  close(hbar);end
Btab = cell2mat(Btabs);
save(matname)
tar(tarname,ImFiles)
gzip(tarname)
delete(tarname)
close(hfmont)
for m = 1:NUMMC
    delete(ImFiles{m})
end
return