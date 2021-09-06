function sigma=sigmaclipping(x,k)
%SIGMACLIPPING
% Reference:
% J.-L. Starck and F. Murtagh, Astronomical Image and Data Analysis, 2nd
% Edition, Springer, p40-p41
x = x-mean(x(:));
if nargin<2
    k=3;
end
MAX_LOOP=500;
x=x(:)-mean(x(:));
sigma=std(x(:));
for n = 1:MAX_LOOP
    l=logical(abs(x(:))<=(k*sigma));
    if ~isempty(x(l))
      sigma_new=std(x(l));
    else
      sigma = 0;
      return
    end
    if sigma_new >= sigma
      return
    end
    sigma=sigma_new;
end
disp('Iteration limit reached without convergence.');
return