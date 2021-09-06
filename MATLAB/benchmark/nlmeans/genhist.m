%% initialize & configure
fontsize = 20;
fontname = 'helvetica';
siglvl = 0.95;
nbins = 40;
facecolor = [0.8,0.8,0.8];
edgecolor = 'black';

%% nlm examples
g = poissrnd(1000*ones(512));
g10 = inlmeans(g,'reps',10,'rot',false);
g100 = inlmeans(g,'reps',100,'rot',false);
g1k = inlmeans(g,'reps',1000,'rot',false);
g5k = inlmeans(g,'reps',5000,'rot',false);

%% without nlm
bmes = mergemcout('./wo_nlm/*.mat');
flux = sort(bmes.flux);
NB = numel(flux);

%% 10 reps
bmes10 = mergemcout('./10reps/*.mat');
flux10 = sort(bmes10.flux);
NB10 = numel(flux10);

%% 100 reps
bmes100 = mergemcout('./100reps/*.mat');
flux100 = sort(bmes100.flux);
NB100 = numel(flux100);

%% 1k reps
bmes1k = mergemcout('./1000reps/*.mat');
flux1k = sort(bmes1k.flux);
NB1k = numel(flux1k);

%% 5k reps
bmes5k = mergemcout('./5000reps/*.mat');
flux5k = sort(bmes5k.flux);
NB5k = numel(flux5k);

%% make list
flst = {flux,flux10,flux100,flux1k,flux5k}; % flux list
tlst = {'No NLMeans','10 reps NLMeans','100 reps NLMeans',...
  '1000 reps NLMeans','5000 reps NLMeans'}; % title list
nlst = [NB, NB10, NB100, NB1k, NB5k]; % number of blobs list
clst = zeros(size(nlst)); % 0.95 cut list
slst = zeros(size(nlst)); % sigma (noise level) list
plst = {'wo_nlm','10reps','100reps','1kreps',...
  '5kreps'}; % print file prefix list
rlst = [0, 10, 100, 1000, 5000]; % repeats numbers list
glst = {g, g10, g100, g1k, g5k}; % observed data (w/o NLMeans) list

%% draw histograms
close all
for k=1:numel(nlst)
  figure('Name',[tlst{k},' H0 distribution'])
  hist(flst{k},nbins)
  h = findobj(gca,'Type','patch');
  set(h,'FaceColor',facecolor,'EdgeColor',edgecolor)
  set(gca,'fontsize',fontsize,'fontname',fontname)
  axis tight
  xlabel('counts','fontsize',fontsize,'fontname',fontname)
  ylabel('frequency','fontsize',fontsize,'fontname',fontname)
  title(tlst{k},'fontsize',fontsize+4,'fontname',fontname)
  v=axis;
  xmid=0.5*(v(1)+v(2));
  xrng=v(2)-v(1);
  v(1)=xmid-xrng*0.5*1.1;
  v(2)=xmid+xrng*0.5*1.1;
  v(4)=v(4)*1.1;
  axis(v)
  clst(k) = flst{k}(round(0.95*nlst(k)));
  line([clst(k) clst(k)],[0,v(4)],'color','r')
  text(clst(k)+0.02*xrng,0.5*v(4),['95% cut: ',int2str(clst(k)),' cts'],...
    'HorizontalAlignment','left',...
    'fontname',fontname,...
    'fontsize',fontsize-4,'color','k')
  print('-depsc2',[plst{k},'.eps'])
  print('-dpng',[plst{k},'.png'])
end

for k=1:numel(nlst)
  figure('Name',[tlst{k},' noise distribution'])
  hist(glst{k}(:),nbins)
  h = findobj(gca,'Type','patch');
  set(h,'FaceColor',facecolor,'EdgeColor',edgecolor)
  set(gca,'fontsize',fontsize,'fontname',fontname)
  axis tight
  xlabel('counts','fontsize',fontsize,'fontname',fontname)
  ylabel('frequency','fontsize',fontsize,'fontname',fontname)
  title(tlst{k},'fontsize',fontsize+4,'fontname',fontname)
  v=axis;
  xmid=0.5*(v(1)+v(2));
  xrng=v(2)-v(1);
  v(1)=xmid-xrng*0.5*1.1;
  v(2)=xmid+xrng*0.5*1.1;
  v(4)=v(4)*1.1;
  axis(v)
  slst(k) = std(glst{k}(:));
  line(1000+[slst(k) slst(k)],[0,v(4)],'color','r')
  text(1000+slst(k)+0.02*xrng,0.75*v(4),['\sigma: ',num2str(slst(k),2),' cts'],...
    'HorizontalAlignment','left',...
    'fontname',fontname,...
    'fontsize',fontsize-4,'color','k')
  print('-depsc2',[plst{k},'_ndist.eps'])
  print('-dpng',[plst{k},'_ndist.png'])
end

%% draw sensitivities plot
figure('Name','NLMeans denoise')
semilogx(max(1,rlst),slst,'k')
set(gca,'fontsize',fontsize,'fontname',fontname)
xlabel('NLMeans reps','fontsize',fontsize,'fontname',fontname)
ylabel('Noise level, cts','fontsize',fontsize,'fontname',fontname)
title('Noise levels with NLMeans denoise',...
  'fontsize',fontsize+4,'fontname',fontname)
print('-depsc2','nlm_nlvl.eps')
print('-dpng','nlm_nlvl.png')

figure('Name','NLMeans sensitivity')
semilogx(max(1,rlst),clst,'k')
set(gca,'fontsize',fontsize,'fontname',fontname)
xlabel('NLMeans reps','fontsize',fontsize,'fontname',fontname)
ylabel('sensitivity, cts','fontsize',fontsize,'fontname',fontname)
title({'Point source detection sentivities ','with NLMeans denoise'},...
  'fontsize',fontsize+4,'fontname',fontname)
print('-depsc2','nlm_95c.eps')
print('-dpng','nlm_95c.png')
