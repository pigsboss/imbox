function [bmes, bimg] = mergemcout(files, fields)
%MERGECOUT

if nargin == 1
  fields = {};
end
flist = ls(files);
if isempty(flist)
  error([files,' not found.'])
end
prefix=files(1:find((files=='/')|(files=='\'),1,'last'));
[nfiles,~] = size(flist);
disp([int2str(nfiles),' files have been found.'])
nb = 0;
disp('Pass 1: scanning for blobs.')
for n = 1:nfiles
  file = [prefix,flist(n,:)];
  disp(['loading ',file])
  s = load(file);
  disp([file, ' loaded.'])
  if isfield(s, 'mcout')
    disp(['MC output has been found in ',file])
    [nmc, ~] = size(s.mcout);
    for m = 1:nmc
      nb = nb + length(s.mcout{m}{1}.flux);
    end
  else
    disp([file, ' is skipped.'])
    continue
  end
end
disp([int2str(nb), ' blobs have been found.'])
[bmes,bimg] = initblobstruct(nb);
fn1 = fieldnames(bmes);
fn2 = fieldnames(bimg);
if ~isempty(fields)
  fn1_idx = false(size(fn1));
  fn2_idx = false(size(fn2));
  for k = 1:numel(fields)
    fn1_idx = fn1_idx | strcmpi(fields{k},fn1);
    fn2_idx = fn2_idx | strcmpi(fields{k},fn2);
  end
else
  fn1_idx = true(size(fn1));
  fn2_idx = true(size(fn2));
end
disp('Selected fields to merge:')
disp(fn1(fn1_idx))
disp(fn2(fn2_idx))
disp('Pass 2: collecting blobs.')
nb = 0;
for n = 1:nfiles
  file = [prefix,flist(n,:)];
  disp(['loading ',file])
  s = load(file);
  disp([file, ' loaded.'])
  if isfield(s, 'mcout')
    disp(['MC output has been found in ',file])
    [nmc, ~] = size(s.mcout);
    for m = 1:nmc
      cnb = length(s.mcout{m}{1}.flux);
      % merge bmes
      for k = 1:numel(fn1)
        if fn1_idx(k)
          bmes.(fn1{k})((nb+1):(nb+cnb),:) = s.mcout{m}{1}.(fn1{k})(:,:);
        end
      end
      % merge bimg
      for k = 1:numel(fn2)
        if fn2_idx(k)
          [bimg((nb+1):(nb+cnb)).(fn2{k})] = s.mcout{m}{2}.(fn2{k});
        end
      end
      nb = nb + cnb;
    end
  else
    disp([file, ' is skipped.'])
    continue
  end
end
return
