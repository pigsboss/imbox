function [txt,dict] = ocr(page,x,y,w,h,dx,dy,NumChars,NumLines,dict)
%OCR Optical character recognition
txt = char(ones(NumLines,NumChars));
page = normalize(page);
B = blob(max(dlt2(page,2),1));
[yi,xi] = find(B);
if isempty(x), x = min(xi);end
if isempty(y), y = min(yi);end
if isempty(w)
  BR = max(B,[],1);
  lenBR = length(BR);
  w = find(BR((x+1):lenBR)==0,1);
end
% if isempty(h)
%   BC = max(B,[],2);
%   lenBC = length(BC);
%   h = find(BC((y+1):lenBC)==0,1);
% end
for l = 1:NumLines
%   txt{l} = char(ones(1,NumChars));
  for n = 1:NumChars
    % current letter in scanned image
    letter = page((y+(l-1)*(h+dy)):(y+(l-1)*(h+dy)+h-1),...
      (x+(n-1)*(w+dx)):(x+(n-1)*(w+dx)+w-1));
    [asc,~] = lookupdict(dict,letter);
    if isempty(asc)
      [dict,asc] = updatedict(dict,letter);
    end
    txt(l,n) = asc;
%     imagesc(letter);axis image;drawnow;colormap gray
  end
end
return

function [asc,acc] = lookupdict(dict,letter)
if isempty(dict)
  asc = [];
  acc = [];
  return;
end
NumEnt = length(dict{1}); % number of entries
letter = normalize(letter);
t = 1e-2;
acc = zeros(NumEnt,1);
sizeL = size(letter);
NP = numel(letter);
for l = 1:NumEnt
  if any(sizeL ~= size(dict{2}{l}))
    acc(l) = sum(sum((imresize(dict{2}{l},sizeL)-letter).^2))/NP;
  else
    acc(l) = sum(sum((dict{2}{l}-letter).^2))/NP;
  end
end
[maxacc,idx] = min(acc);
if maxacc < t
  asc = dict{1}(idx);
else
  asc = [];
end
return

function [dict,asc] = updatedict(dict,letter)
letter = normalize(letter);
if isempty(dict), dict = cell(2,1);end
NumEnt = length(dict{1});
dict{1} = [dict{1};'A'];
dict{2} = [dict{2};{letter}];
imagesc(letter);axis image;colormap gray;drawnow
userinput = inputdlg({'What is the character?'},'New character found');
if ~isempty(userinput)
  asc = userinput{1};
%   if any(dict{1} == asc)
%     hf=figure('Name','Compare');
%     subplot(1,2,1);imagesc(dict{2}{dict{1}==asc});axis image;colormap gray;title('Exist letter');drawnow
%     subplot(1,2,2);imagesc(letter);axis image;colormap gray;title('Current letter');drawnow
%   end
  dict{1}(NumEnt+1) = asc;
else
  error('User canceled.')
end
return

function I = normalize(I)
maxI = max(I(:));
minI = min(I(:));
if maxI ~= minI
  I = (I-minI)/(maxI-minI);
else
  I = I - minI;
end
return
