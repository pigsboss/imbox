function printbtab(Btab,TABNAME)
%PRINTBTAB print blob table.
% Btab is returned by blobex.m
[NUMPARM,NUMBLOB] = size(Btab);
RLAB = cell(1,NUMBLOB);
for l = 1:NUMBLOB
    RLAB{l} = [int2str(l),' '];
end
RLAB = cell2mat(RLAB);
switch NUMPARM
    case 3
        printmat(Btab',TABNAME,RLAB,'Flux X Y');
    case 5
        printmat(Btab',TABNAME,RLAB,'Flux X Y RX RY');
    case 7
        printmat(Btab',TABNAME,RLAB,'Flux X Y RX RY PHI RMSE');
    otherwise
        error('Unsupported blob table.')
end
return
