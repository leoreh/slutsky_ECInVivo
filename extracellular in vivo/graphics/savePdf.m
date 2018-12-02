function savepdf(filename, basepath, f)

% saves figure as pdf with full vector data.
% file saved in basepath/graphics. if /graphics does not exit, creates it.
%
% INPUT
%   filename        name of figure
%   basePath        path of session data (without Graphics)
%   fig             handle to current figure
%
% 24 nov 18 LH

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