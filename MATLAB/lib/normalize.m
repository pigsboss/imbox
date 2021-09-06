function g = normalize(f)
%NORMALIZE Normalize input signal.
%
%Syntax
% g = normalize(f)
%
%Description
% f is input signal.
%
%Return
%            f
% g = ---------------
%      f_max - f_min

fmax = max(f(:));
fmin = min(f(:));
frng = fmax - fmin;
if frng >= eps
  g = f/frng;
else
  error('Input signal cannot be normalized.');
end
return