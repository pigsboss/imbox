function ddh0wonlm(NUMDD,prefix)
if nargin == 1
  prefix='./';
end
load HEPSF
BGLVL = 1000; % counts
NUMPX = 512;
NUMIT = 50;
filestr = [prefix,'wo_nlm_x',int2str(NUMDD),'_',datetimefilename];
% NUMDD = 3;
scale = 1:0.1:3;
h = imresize(mean(PSF512.PSF,3),0.5);
h = h/sum(h(:));

rng('shuffle')
mcout = cell(NUMDD,1); % mc output, as {{bmes, bimg}, {bmes, bimg}, ...}
hbar = waitbar(0,'Running MC...');
tMC = tic;
diary([filestr,'.log'])
diary on
for m = 1:NUMDD
  disp(['MC ',int2str(m)])
  g = poissrnd(ones(NUMPX)*BGLVL);
  f = dd(g,h,NUMIT,'LL',BGLVL,'verbose',0);
  [bmes, bimg] = ellipsems(f,scale);
  [bmes, bimg] = blobchk(bmes, bimg);
  mcout{m} = {bmes, bimg};
  tElapsed = toc(tMC);
  tRemain = ((NUMDD-m)/m)*tElapsed;
  waitbar(m/NUMDD,hbar,['Running MC, ',int2str(round(tRemain)),' seconds left...']);
end
diary off
close(hbar)
save([filestr,'.mat'],'mcout','-v7.3')
return