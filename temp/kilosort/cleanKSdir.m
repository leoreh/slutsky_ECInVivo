function cleanKSdir(basepath)

% cleans basepath after spike sorting (kilosort) and initializing mat
% format (cellExplorer). deletes temp_wh.dat file and moves ks output to ks
% folder. note that to run cellExplorer (ProcessCellMetrics) again after
% cleaning, must change relative path in session file and copy xml to ks
% folder.
% 
% INPUT:
%   basepath    string. path to recording folder {pwd}. 
%
% TO DO LIST:
%   # apply to klustakwik
% 
% 24 may 20 LH

fprintf('\nremoving temp files and arranging ks directory \n')

kspath = fullfile(basepath, 'ks');
mkdir(kspath)

npyfiles = dir(fullfile(basepath, '*npy'));
tsvfiles = dir(fullfile(basepath, '*tsv'));
for i = 1 : length(npyfiles)
    movefile(fullfile(basepath, npyfiles(i).name), kspath)
end
for i = 1 : length(tsvfiles)
    movefile(fullfile(basepath, tsvfiles(i).name), kspath)
end

movefile(fullfile(basepath, 'chanMap.mat'), kspath)
movefile(fullfile(basepath, '.phy'), kspath)
movefile(fullfile(basepath, 'phy.log'), kspath)
movefile(fullfile(basepath, 'params.py'), kspath)
movefile(fullfile(basepath, 'rez.mat'), kspath)
delete(fullfile(basepath, 'temp_wh.dat'))

end

% EOF