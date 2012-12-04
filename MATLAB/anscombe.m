function g=anscombe(x)
%ANSCOMBE Anscombe transform
% transform poissonian data so that the variance is approximately 1 whatever the
% mean. It's valid provided that the mean of the poissonian data is larger than
% 4.
%
% x -> 2 * sqrt (x + 3/8)
% Reference:
% http://en.wikipedia.org/wiki/Anscombe_transform
g=2*sqrt(x+3/8);
return