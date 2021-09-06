function showImStat(TitleStr,I,BG)
disp( '++++++++++++++++++++++++++++++++++')
disp(TitleStr);
disp( '----------------------------------')
disp(['                 max: ',num2str(max(I(:)))]);
disp(['                 min: ',num2str(min(I(:)))]);
disp(['                mean: ',num2str(mean(I(:)))]);
disp(['  standard deviation: ',num2str(std(I(:)))]);
disp(['                 SNR: ',num2str(snr(I,BG)),' dB']);
disp( '++++++++++++++++++++++++++++++++++')
return

function r = snr(I,BG)
s = I-BG;
mu = mean(s(:));
sigma = sigmaclipping(I,3);
if mu<=0
  r = -inf;
elseif sigma<=0
  r = inf;
else
  r = 10*log10(mu/sigma);
end
return

