function newDir = cp_basepathFiles(basepath, destDir, fileExt)

% Function to copy files with a specific extension from one directory to another
% basepath  : Source directory path (string)
% destDir   : Destination directory path (string)
% fileExt   : string array of file extensions to filter by (e.g., '.txt', '.jpg')

% Ensure the source directory exists
if ~isfolder(basepath)
    error('Source directory does not exist.');
end
[~, basename] = fileparts(basepath);
% Create the destination directory if it doesn't exist

newDir = fullfile(destDir, basename);
if ~isfolder(newDir)
    mkdir(newDir);
end

% Get a list of files with the specified extension
for iext = 1 : length(fileExt)
    filePattern = fullfile(basepath, ['*' fileExt{iext} '*']);
    files = dir(filePattern);

    % Loop through each file and copy it to the destination directory
    for ifile = 1 : length(files)
        srcFile = fullfile(basepath, files(ifile).name);
        destFile = fullfile(newDir, files(ifile).name);
        copyfile(srcFile, destFile);
    end

    fprintf('Copied %d %s files from %s to %s.\n', length(files), fileExt{iext}, basepath, newDir);

end

return

% EOF

% example call with file extensions
destDir = 'C:\Users\Leore\Downloads\ketamine';
fileExt = [...
    ".cell_metrics.";...
    ".datInfo.";...
    ".fr.";...
    ".session.";...
    ".spikes.cellinfo.";...
    ".units.";...
    ".xml.";...
    ".sleep_states.";...
    ".sleep_sig.";...
    ];
[basepaths] = ket_sessions('ket');
for ifile = 1 : length(basepaths)
    basepath = basepaths{ifile};
    cd(basepath)
    newDir = cp_basepathFiles(basepath, destDir, fileExt);
    cd(newDir)
end


