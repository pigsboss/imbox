function varargout = inlmerge(matname)
%INLMERGE Merge inlmeans output mat files.
files = dir([matname,'*.mat']);
numfiles = length(files);
mat = load(files(1).name);
sizeI = size(mat.img);
s = zeros(sizeI);
z = zeros(sizeI);
for l = 1:numfiles
  disp(['file: ',files(l).name])
  mat = load(files(l).name);
  if mat.mc == true
    disp('denoised by monte-carlo method.')
    disp(['repeats: ',int2str(mat.reps)])
%     disp(['method noise: ',num2str(std(mat.dimg(:) - mat.img(:)))])
    s = s + mat.simg;
    z = z + mat.z;
  end
end
d = s./z;
if nargout == 0
  figure('Name','Denoised image, merged');imagesc(d);axis image;
  figure('Name','Similarity map, merged');imagesc(log10(z));axis image;
else
  varargout{1} = d;
  varargout{2} = s;
  varargout{3} = z;
end
% disp(['method noise: ',num2str(std(d(:) - mat.img(:)))])
return