function [Stab,Btab] = deconvstd(k,PSF,hdeconv,NUMIT,SmoothScale,NUMMC)
%DECONVNULL compute histogram of null distribution of blobs deconvolved by
%specified program.
%
%Input arguments:
%   k             is the significance measure of the standard point source,
%   namely, the flux of simulated standard point source is k * std(BG).
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

sizePSF = size(PSF);
Bgrd = 100;
Mdl = ones(sizePSF(1:2))*Bgrd;
Mdl(1,1) = sqrt(Bgrd)*k/max(PSF(:))+Bgrd;
Mdl = fftshift(Mdl);
rng('shuffle')
Btabs = cell(1,NUMMC); % all blobs
Stabs = cell(1,NUMMC); % source
ImFiles = cell(1,NUMMC);
filename = [func2str(hdeconv),'_',num2str(k),'sig_',int2str(NUMIT),'x',num2str(SmoothScale),...
    'x',int2str(NUMMC)];
matname = [filename,'.hist.mat'];
tarname = [filename,'.profile.tar'];
hfmont = figure('Name',[num2str(k),' sigma source, all blobs']);
m = 1;
NUMIM = 0;

hbar = waitbar(0,'Running benchmark...','Name',[num2str(k),' sigma source deconvolution benchmark']);
t = tic;

while m <= NUMMC
    g = poissrnd(imconv(Mdl,PSF));
    f = hdeconv(g,PSF,NUMIT,mean(g(:)));
    [Btabs{m},Bimg] = blobex(f,SmoothScale);
    if ~isempty(Btabs{m})
        [~,idx] = max(Btabs{m}(1,:)); % search for brightest blob
        Stabs{m} = Btabs{m}(:,idx); % the brightest blob is taken as the source
        [H,W,~] = size(Bimg);
        figure(hfmont)
        hmont = montage(reshape(Bimg,H,W,1,[])*255,hot(256));drawnow
        ImFiles{m} = [filename,'_mon_',int2str(m),'.png'];
        imwrite(uint8(get(hmont,'CData')),hot(256),...
            ImFiles{m},'png')
        m = m+1;
        NUMIM = NUMIM+1;
    end
    tElapsed = toc(t);
    tRemain = ((NUMMC-m)/m)*tElapsed;
    waitbar(m/NUMMC,hbar,['Running benchmark, ',...
        num2str(m/NUMMC*100),'% completed, ',...
        int2str(round(tRemain)),' seconds left...']);
end
if exist('hbar','var')==1,  close(hbar);end
Btab = cell2mat(Btabs(1:NUMIM));
Stab = cell2mat(Stabs(1:NUMIM));
save(matname)
tar(tarname,ImFiles(1:NUMIM))
gzip(tarname)
delete(tarname)
close(hfmont)
for m = 1:NUMIM
    delete(ImFiles{m})
end
return
