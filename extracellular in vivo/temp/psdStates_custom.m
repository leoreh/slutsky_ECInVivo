% gen params --------------------------------------------------------------
session = CE_sessionTemplate(pwd, 'viaGUI', false,...
    'force', true, 'saveVar', true);      
basepath = session.general.basePath;
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;
[~, basename] = fileparts(basepath);

[cfg_colors, cfg_names, ~] = as_loadConfig([]);
sstates = [1, 4, 5];

% arrange data ------------------------------------------------------------
% datInfo and labels
infofile = dir('*datInfo*');
load(infofile.name)
load([basename, '.AccuSleep_labels.mat']) % load labels. assumes 1 s epochLen

% part 1 ------------------------------------------------------------------
start = 0;
dur = floor(datInfo.nsamps(1) / fs);
labelsTemp = labels(1 : dur);

% load
EEG = double(bz_LoadBinary([basename, '.lfp'], 'duration', dur,...
    'frequency', 1250, 'nchannels', nchans, 'start', 0,...
    'channels', [4 : 7], 'downsample', 1));
EEG = mean(EEG, 2);

% psd
[psd1, faxis, ~] = psd_states('eeg', EEG, 'emg', [],...
    'labels', labelsTemp, 'fs', 1250, 'graphics', true, 'sstates', sstates);

% part 2 ------------------------------------------------------------------
start = ceil(datInfo.nsamps(1) / fs);
dur = Inf;
labelsTemp = labels(start : end);

% load
EEG = double(bz_LoadBinary([basename, '.lfp'], 'duration', dur,...
    'frequency', 1250, 'nchannels', nchans, 'start', 0,...
    'channels', [4 : 7], 'downsample', 1));
EEG = mean(EEG, 2);

% psd
[psd2, ~, ~] = psd_states('eeg', EEG, 'emg', [],...
    'labels', labelsTemp, 'fs', 1250, 'graphics', true, 'sstates', sstates);

% arrange comparison and plot ---------------------------------------------

% smooth and normalize
smf = 21;
psdSmooth1 = smoothdata(psd1 ./ sum(psd1, 2), 2, 'movmean', smf);
psdSmooth2 = smoothdata(psd2 ./ sum(psd2, 2), 2, 'movmean', smf);

fh = figure;
for istate = 1 : length(sstates)
    subplot(1, 3, istate)
    ph = plot(faxis, psdSmooth1(istate, :), 'k', 'LineWidth', 1);
    hold on
    ph = plot(faxis, psdSmooth2(istate, :), 'b', 'LineWidth', 1);
    title(cfg_names{sstates(istate)})
    xlim([0 120])
    xlabel('Frequency [Hz]')
    ylabel('norm PSD')
%     set(gca, 'YScale', 'log')
    legend('Before', 'After')
end


