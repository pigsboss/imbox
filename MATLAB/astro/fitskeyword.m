function [keyval,keycomment]=fitskeyword(fitsfilename,keyname)
%FITSKEYWORD get keyword value and comment of input fits file.
info=fitsinfo(fitsfilename);
keywords=info.PrimaryData.Keywords;
k = strcmpi(keyname,keywords(:,1));
keyval = keywords{k,2};
keycomment = keywords{k,3};
return