function Y = fajita(X,M)
%FAJITA FAst adJacent InpainTing Algorithm.
%INPUTS:
%  X           is image with sockets.
%  M (boolean) is mask of the image X. false in M labels sockets.
%RETURN:
%  Y           is inpainted image.
typex = str2func(class(X));
Y = double(X);
while any(~M(:))
  [Y,M] = adjaverage(Y,M);
end
Y = typex(Y);
return

function [Y,M] = adjaverage(X,M)
%ADJAVERAGE Adjacent pixel averaging.
sft = [-1,-1; -1,0; -1,1; 0,-1; 0,1; 1,-1; 1,0; 1,1];
wgt = sqrt(sum(sft.^2,2));
sizeX = size(X);
ndX = length(sizeX);
sizeM = size(M);
ndM = length(sizeM);
if ndM==2 && ndX==3
  W = zeros(sizeX);
  for k = 1:sizeX(3)
    W(:,:,k) = M;
  end
  M = logical(W);
end
W = zeros(size(M));
Y = zeros(size(X));
for k = 1:8
  W = W+zeroshift(M,sft(k,:))*wgt(k);
  Y = Y+zeroshift(X,sft(k,:))*wgt(k);
end
Y((W>0.7)&(~M)) = Y((W>0.7)&(~M))./W((W>0.7)&(~M));
Y(M) = X(M);
M = (W>0.7);
return