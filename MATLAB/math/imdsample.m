function X = imdsample(X,N)
%IMDSAMPLE
for n = 1:N
    szX = size(X);
    X = padarray(X,ceil(szX/2)*2-szX,'replicate','post');
    X = X(:,1:2:end)+X(:,2:2:end);
    X = X(1:2:end,:)+X(2:2:end,:);
end
return