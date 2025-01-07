
% files
basepath = 'D:\Data\RA_WT2';
cd(basepath)
[~, basename] = fileparts(basepath);
eegfile = [basename '.AccuSleep_EEG.mat'];
emgfile = [basename '.AccuSleep_EMG.mat'];

% params
fs = 1000;
SR = 512;
recDur = length(scoring) / 1000;
emg_orig = EEG2_1K;
eeg_orig = EMG2_1K;
boutLen = 2.5;

% cut to length
idx = 1 : 12 * 60 * 60 * fs;
EEG = EEG2_1K(idx);
EMG = EMG2_1K(idx);
labels = scoring(idx);

% remap to accusleep
labels(labels == 6) = 1;
labels(labels == 5) = 3;
labels(labels == 2) = 5;
labels(labels == -1) = 2;
labels(labels == -3) = 5;

% resample lables
labels = [interp1([1 : length(labels)] / 1000, labels,...
    [1 : idx(end) / fs / boutLen] * boutLen, 'linear', 'extrap')]';
labels = round(labels);

% inspect
AccuSleep_viewer(EEG, EMG, 1000, boutLen, labels, [])

% save files
save(eegfile, 'EEG')
save(emgfile, 'EMG')
save([basename, '.AccuSleep_labelsMan.mat'], 'labels')


% sleep states
ss = accusleep_wrapper('basepath', basepath, 'cleanRec', [],...
    'SR', 1000, 'boutLen', 2.5, 'recSystem', 'tdt', 'calfile', [],...
    'lfpCh', [13 : 16], 'emgCh', [], 'viaGui', false,...
    'forceCalibrate', true, 'inspectLabels', true, 'saveVar', true,...
    'forceAnalyze', true, 'forceLoad', true);
%         AccuSleep_viewer(EEG, EMG, ss.fs, 2.5, ss.labels, []);





