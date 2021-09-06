function sq = rect2square(rect)
sizeR = size(rect);
if sizeR(1)>sizeR(2)
  D = sizeR(1) - sizeR(2);
  if mod(D,2) == 0
    sq = padarray(rect,[0 D/2]);
  else
    D1 = (D+1)/2;
    D2 = D - D1;
    sq = padarray(rect,[0 D1],0,'pre');
    sq = padarray(sq,[0 D2],0,'post');
  end
elseif sizeR(1)<sizeR(2)
  D = sizeR(2) - sizeR(1);
  if mod(D,2) == 0
    sq = padarray(rect,[D/2 0]);
  else
    D1 = (D+1)/2;
    D2 = D - D1;
    sq = padarray(rect,[D1 0],0,'pre');
    sq = padarray(sq,[D2 0],0,'post');
  end
else
  sq = rect;
end
return