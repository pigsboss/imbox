subdir = {'lib','wavelet','math','detector','ellipse','gui','deconv',...
    'video','benchmark','ocr','astro','srcsprt','denoise','inpainting',...
    'bootstrap'};
os = computer;
switch os
    case {'PCWIN64','PCWIN'}
        for k = 1:length(subdir)
            path([pwd,'\',subdir{k}],path)
        end
    otherwise
        for k = 1:length(subdir)
            path([pwd,'/',subdir{k}],path)
        end
end
clear k subdir os