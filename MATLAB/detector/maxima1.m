function M = maxima1(f,dim)
%MAXIMA1 Find all local maxima along the specified dimension.
%INPUT
% f    is the n-dim array.
% dim  is the dimension specified.
%RETURN
% M    is a boolean matrix of the same size as f.
if nargin==1
    dim = 1;
end
sz = size(f);
if dim>length(sz)
    error('Specified dimension exceeds upper limit.')
end
sftp = zeros(size(sz));
sftp(dim) = 1;
sftm = zeros(size(sz));
sftm(dim) = -1;
M = (f>=circshift(f,sftm));
M = ((f>=circshift(f,sftp)) & M);
return