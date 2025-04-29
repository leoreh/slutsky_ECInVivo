% Define base paths and file extensions
basepath1 = 'D:\OneDrive - Tel-Aviv University\Data\baclofen\lh142';
basepath2 = 'E:\Data\lh142';
extensions = {'.fr', '.cell_metrics.cellinfo', '.session', '.sleep_states', '.spikes.cellinfo', '.units'};  % specify desired extensions

% List subfolders in basepath1 (excluding . and ..)
d = dir(basepath1);
subFolders = {d([d.isdir]).name};
subFolders = subFolders(~ismember(subFolders, {'.','..'}));

% Loop through each subfolder
for i = 1:length(subFolders)
    folderName = subFolders{i};
    srcDir = fullfile(basepath2, folderName);
    dstDir = fullfile(basepath1, folderName);
    if ~exist(srcDir, 'dir')
        warning('Source folder %s not found.', srcDir);
        continue;
    end
    for j = 1:length(extensions)
        fileName = [folderName, extensions{j} ,'.mat'];
        srcFile = fullfile(srcDir, fileName);
        dstFile = fullfile(dstDir, fileName);
        if exist(srcFile, 'file')
            copyfile(srcFile, dstFile);
        else
            warning('File %s not found.', srcFile);
        end
    end
end
