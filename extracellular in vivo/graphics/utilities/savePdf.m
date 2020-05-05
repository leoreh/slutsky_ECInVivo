function savePdf(filename, basepath, f, tiff)

% saves figure as pdf with full vector data.
% file saved in basepath/graphics. if /graphics does not exit, creates it.
% can also save tiff version of figure
%
% INPUT
%   filename        name of figure
%   basePath        path of session data (without Graphics)
%   fig             handle to current figure
%   tiff            logical. save tiff or not
%
% 24 nov 18 LH

cd(basepath)
fullpath = [basepath, '\graphics'];
strdate = date;
filename = fullfile(fullpath, [filename, '_', strdate]);

if ~exist(fullpath, 'dir')
    sprintf('Creating Graphics folder in %s', basepath)
    mkdir('graphics')
end

orient(f, 'landscape')
print(f, filename, '-dpdf', '-bestfit', '-painters')

end