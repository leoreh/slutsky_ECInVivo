function fileList = as_fileLists(basepaths)

% arranges the fileList cell array for training a network. basepaths is a 1
% x n cell of fullpath folder names. fileList is a 3 x n cell of fullpath
% filenames (eeg, emg, and labels)

for idir = 1 : length(basepaths)
    
    [~, basename] = fileparts(basepaths{idir});    
    labelsmanfile = [basename, '.AccuSleep_labelsMan.mat'];
    eegfile = [basename '.AccuSleep_EEG.mat'];
    emgfile = [basename '.AccuSleep_EMG.mat'];    
    
    fileList{idir, 1} = fullfile(basepaths{idir}, eegfile);
    fileList{idir, 2} = fullfile(basepaths{idir}, emgfile);
    fileList{idir, 3} = fullfile(basepaths{idir}, labelsmanfile);
end

return