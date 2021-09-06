function d = roulette(varargin)
%ROULETTE Monte-Carlo resample
%SYNTAX
% d = roulette(sample,weight,resample_size)
% d = roulette(sample,[],resample_size)
% d = roulette(sample,weight)
% d = roulette(sample)
%INPUT
% sample            is the sample.
% weight (optional) is the weight of the sample.
% Default is uniform.
% resample_size     is the size vector of resampled data.
%RETURN
% d is resampled data.
%
narginchk(1,3);
sample = [];
weight = [];
resample_size = [];
switch nargin
  case 1
    sample = varargin{1};
  case 2
    sample = varargin{1};
    weight = varargin{2};
  case 3
    sample = varargin{1};
    weight = varargin{2};
    resample_size = varargin{3};
end
if isempty(sample), d = []; return; end;
sample = sample(:);
N = length(sample);
if isempty(resample_size), resample_size = 1; end;
if isempty(weight)
  idx = round(rand(resample_size)*N+0.5);
  d = sample(idx);
else
  if any(weight<0), error('Weight must be non-negative.');end
  weight = cumsum(weight(:));
  ul = weight;
  ll = [0;weight(1:(N-1))];
  clear weight
  idx = min(max(rand(resample_size)*max(ul),0),max(ul));
  idx = idx(:);
  L = length(idx);
  d = zeros(L,1);
  for l=1:L
    d(l) = sample((idx(l)>=ll) & (idx(l)<ul));
  end
end
return
