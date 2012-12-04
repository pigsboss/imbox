function p=lastpow2(A)
%LASTPOW2 Last lower (or equal) power of 2
% 2^p <= abs(A)
p=floor(log(double(abs(A)))/log(2.00));
return