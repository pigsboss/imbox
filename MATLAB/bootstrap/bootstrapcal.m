function sigma = bootstrapcal(K,N,L)
%BOOTSTRAPCAL Bootstrap precision calibration.
M = 2.^(1:L);
sigma = zeros(L,K);
tic
for l=1:L
  m = M(l);
  sigma_row = zeros(1,K);
  parfor k=1:K
    sigma_row(k) = bootstrap(N,m);
  end
  sigma(l,:) = sigma_row;
end
t = toc;
disp(['time consumed: ',num2str(t),' seconds.'])
save([int2str(K),'x',int2str(N),'x',int2str(L)])
return
