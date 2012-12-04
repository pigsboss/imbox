function p=nextpow2(A)
%NEXTPOW2 Next higher (or equal) power of 2
% 2^p >= abs(A)
p=ceil(log(double(abs(A)))/log(2.00));
return