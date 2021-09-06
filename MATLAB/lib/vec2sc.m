function s = vec2sc(v,m,x)
%VEC2SC convert vector array v to scalar array s, according to the mapping
%map.
%
%Syntax
% s = vec2sc(v,map)
% s = vec2sc(v,map,x)
%
% v   is a vector array. Each element of v is a vector with L scalar
% elements.
% s   is a scalar array. s contains the same number of scalar elements as
% the number of vectors v does. All the scalars are in a discrete range.
% x   is a vector defines the range of scalars. If x only contains two
% elements then the two elements are considered as x_min and x_max which
% are the lower and upper limits of the range. The stepsize is calculated
% according to size of m. Otherwise x represents every nodes in the
% discrete range.
% m   is a NxL array. L is the length of vectors in v while N is the length
% of the discrete range.
%
%Description
% m maps a scalar in the discrete range to a vector, e.g., a colormap
% maps a scalar to a RGB color vector. To find the original scalar with the
% vector as well as the mapping known, we find the closest vector in the
% mapping by calculating an Euclidean distance defined in the vector space
% first.
%
%Return
%See Also
szBufMB = 2048;
szV = size(v);
szM = size(m);
if nargin==2
  x = 1:szM(1);
end
if szV(ndims(v)) ~= szM(2)
  error('The number of columns of map must equal to the size of the vector.')
end
lenVec = szV(ndims(v));
szS = szV(1:(ndims(v)-1));
v = reshape(v,[],lenVec);
[numVec,~] = size(v);
[lenMap,~] = size(m);
disp(['Size of vector: ', int2str(lenVec)])
disp(['Length of map: ', int2str(lenMap)])
disp(['I/O buffer size: ', int2str(szBufMB), 'MB'])
szBuf = szBufMB*1024*1024/8;
lenBuf = floor(szBuf/(lenMap*lenVec));
numBuf = ceil(numVec/lenBuf);
s = zeros(numVec,1);
mBuf = reshape(reshape(m.',[],1)*ones(1,min(numVec,lenBuf)),lenVec,lenMap,[]);
for n=1:numBuf
  curNumVec = min((n*lenBuf),numVec);
  vBuf = zeros(lenVec,lenMap,curNumVec-((n-1)*lenBuf));
  for c=1:lenVec
    vBuf(c,:,:) = ...
      reshape(ones(lenMap,1)*...
      reshape(v((1+(n-1)*lenBuf):curNumVec,c),1,[]),1,lenMap,[]);
  end
  szVB = size(vBuf);
  if ismatrix(vBuf)
    d = reshape(sum((vBuf - mBuf(:,:,1)).^2,1),lenMap,[]);
  else
    d = reshape(sum((vBuf - mBuf(:,:,1:szVB(3))).^2,1),lenMap,[]);
  end
  [~,idx] = min(d,[],1);
  s((1+(n-1)*lenBuf):min((n*lenBuf),numVec)) = idx;
end
s = reshape(x(s),szS);
return