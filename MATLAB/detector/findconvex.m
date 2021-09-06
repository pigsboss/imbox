function [xc, cidx, score] = findconvex(x,y,PRECISION)
%FINDCONVEX Find convex and concave (downward convex) points of 1D
%function.
%
%Syntax
% [xc, cidx, score] = findconvex(x, y, PRECISION)
% [xc, cidx, score] = findconvex(x, y)
% [xc, cidx, score] = findconvex(y)
%
% x is 1D vector. Elements of x are equally spaced.
% y has the same length as x.
% PRECISION is a positive scalar.
%
%Description
% The point (x_c, f(x_c)) on curve y = f(x) is a convex point of the curve,
% if and only if
%  / f(x_c) >= f(x), for all a <= x <= b,
% <  f'(x_c) == 0, and
%  \ f"(x_c) < 0.
%
%Return
% xc is 1D vector. Elements of xc are x values of convex (concave) points.
% cidx is indices of xc in the original vector x.
% score is convexity score of convex points.
% xc, cidx and score have the same length.
%
%
%See Also
%findroot, normalize

narginchk(1,3)
PRECISION_d = 1e-6; %Default precision
switch nargin
  case 1
    y = x;
    x = 1:length(y);
    PRECISION = PRECISION_d;
  case 2
    PRECISION = PRECISION_d;
end
N = length(x);
f1 = normalize(imconv(y,[0.5,0,-0.5],'symmetric') /...
    mean(diff(x))); %Derivative series
f2 = normalize(imconv(y,[1,-2,1],'symmetric') /...
    mean(diff(x)));%2nd order derivative series
%find root of equation f'(x) == 0:
[x1, idx] = findroot(x,f1,PRECISION);
if isempty(x1)
    disp('No convex or concave point found with current precision.')
    xc = [];
    cidx = [];
    score = [];
    return
end
nidx = min(max(round(idx),1),N);
score = y(nidx)-f2(nidx); %Score of convexity
[score, sidx] = sort(score, 'descend');
xc = x1(sidx);
cidx = idx(sidx);
return