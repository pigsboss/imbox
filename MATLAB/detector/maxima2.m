function M = maxima2(f,dim)
%MAXIMA2 Find all local maxima of planes parallel with the specified 
%reference one.
%INPUT
% f    is the n-dim array.
% dim  is the reference plane. For example dim = [1 3] stands for the plane
% parallel with the 1st dimension and the 3rd dimension of the array.
%RETURN
% M    is a boolean matrix of the same size as f.
if nargin==1
    dim = [1 2];
end
sz = size(f);
if max(dim)>length(sz)
    error('Specified dimension exceeds upper limit.')
end
sft = zeros(8,length(sz));
k = 1;
M = true(size(f));
for m = -1:1
    for n = -1:1
        if (m~=0) || (n~=0)
            sft(k,dim(1)) = m;
            sft(k,dim(2)) = n;
            M = M & (f >= circshift(f,sft(k,:)));
            k = k+1;
        end
    end
end
return