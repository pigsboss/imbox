function M = maxima3(f)
%MAXIMA3 Find all local maxima of given 3-d array.
%INPUT
% f    is the 3-d array.
%RETURN
% M    is a boolean matrix of the same size as f.
sft = zeros(8,3);
k = 1;
M = true(size(f));
for l = -1:1
    for m = -1:1
        for n = -1:1
            if (l~=0) || (m~=0) || (n~=0)
                sft(k,1) = l;
                sft(k,2) = m;
                sft(k,3) = n;
                M = M & (f >= circshift(f,sft(k,:)));
                k = k+1;
            end
        end
    end
end
return
