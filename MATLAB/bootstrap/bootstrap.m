function sigma = bootstrap(N,M)
%BOOTSTRAPDEMO Bootstrap method demonstration
%INPUT
%  N  is the size of the sample
%  M  is the number of bootstrap samples
s = normrnd(0,1,[N 1]);
b = zeros(M,1);
for m = 1:M
  b(m) = mean(roulette(s,[],[N 1]));
end
sigma = std(b);
% disp(['mean of estimator: ',num2str(mean(b))]);
% disp(['variance of estimator: ',num2str(sigma^2)]);
% disp(['expected mean: ',num2str(0)]);
% disp(['expected variance: ',num2str(1/N)]);
return
