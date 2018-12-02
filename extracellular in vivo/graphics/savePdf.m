function savePdf(filename, filepath, f)

% saves figure as pdf with full vector data.
% file saved in basePath/Graphics. if /Graphics does not exit, creates it.
%
% INPUT
%   filename        name of figure
%   basePath        path of session data (without Graphics)
%   fig             handle to current figure
%
% 24 nov 18 LH

fullpath = [filepath, '\Graphics'];
strdate = date;
filename = fullfile(fullpath, [filename, '_', strdate]);

if ~exist(fullpath, 'dir')
    sprintf('Creating Graphics folder in %s', filepath)
    mkdir('graphics')
end

orient(f, 'landscape')
print(f, filename, '-dpdf', '-bestfit', '-painters')

end