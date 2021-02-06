% session info (cell explorer foramt)
session = CE_sessionTemplate(pwd, 'viaGUI', false,...
    'force', true, 'saveVar', true);      
basepath = session.general.basePath;
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;
[~, basename] = fileparts(basepath);

% emg
load([basename, '.EMG2.datInfo.mat'])
emgname = [basename '.emg.dat'];
EMG = double(bz_LoadBinary(emgname, 'duration', Inf,...
    'frequency', datInfo.fs, 'nchannels', 1, 'start', 0,...
    'channels', 1, 'downsample', 1));
% lfp
lfpname = [basename, '.lfp'];
LFP = double(bz_LoadBinary(lfpname, 'duration', Inf,...
    'frequency', 1250, 'nchannels', nchans, 'start', 0,...
    'channels', 13, 'downsample', 1));

% resmaple
newFs = 1000;
EMG = accurateResampling(double(EMG), datInfo.fs, newFs);
LFP = accurateResampling(double(LFP), 1250, newFs);
% ALT 2
EEG = [interp1([1 : length(LFP)] / 1000, LFP, [1 : length(EMG)] / newFs,...
    'pchip')]';

% dan
mkdir('AccuSleep')
save('AccuSleep\EEG.mat', 'EEG')
save('AccuSleep\EMG.mat', 'EMG')
labels = ones(1, round(length(EMG) / newFs / 2.5)) * 4;
save('AccuSleep\labels.mat', 'labels')
AccuSleep_GUI
AccuSleep_viewer(EEG, EMG, newFs, 2.5, [], [basepath, '\AccuSleep'])

% documentation
doc AccuSleep_instructions

% buzsaki 
emg.orig = double(bz_LoadBinary(filename, 'nChannels', 4, 'channels', 2,...
    'precision', 'int16', 'start', 0, 'frequency', datInfo.fs));
emg.data = bz_NormToRange(zscore(emg.orig), [0 1]);
emg.tstamps = [1 : length(emg.orig)] / datInfo.fs;
save([basename '.emg.mat'], 'emg')
emg.timestamps = emg.tstamps;
EMGFromLFP = emg;
save([basename '.EMGFromLFP.LFP.mat'], 'EMGFromLFP')
badch = setdiff([session.extracellular.electrodeGroups.channels{:}],...
    [session.extracellular.spikeGroups.channels{:}]);
SleepScoreMaster(basepath, 'noPrompts', true, 'rejectChannels', badch)
TheStateEditor(fullfile(basepath, basename))