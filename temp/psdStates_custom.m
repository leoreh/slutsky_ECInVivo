
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare two times from the same session
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gen params --------------------------------------------------------------

fs = 1250;
sstates = [1, 4, 5];
[cfg_colors, cfg_names, ~] = as_loadConfig([]);

% load data ---------------------------------------------------------------

% rec
basepath = 'D:\Data\lh93\lh93_210811_102035';
session = CE_sessionTemplate(basepath, 'viaGUI', false,...
    'force', true, 'saveVar', true);      
nchans = session.extracellular.nChannels;
[~, basename] = fileparts(basepath);
cd(basepath)

% labels
load([basename, '.AccuSleep_labels.mat'])

% dat info
infoname = dir('*datInfo*');
load(infoname(end).name)
if isfield(datInfo, 'fs')
    fs = datInfo.fs;
else
    fs = 20000;
end

% data
eeg = double(bz_LoadBinary([basename, '.lfp'], 'duration', Inf,...
    'frequency', 1250, 'nchannels', nchans, 'start', 0,...
    'channels', [1 : 4], 'downsample', 1));
eeg = mean(eeg, 2); 

% separate ----------------------------------------------------------------

idx1 = [1 : floor(datInfo.nsamps(1) / fs - 10)];
idx2 = ceil(datInfo.nsamps(1) / fs + 10) : floor(length(eeg) / 1250);
labels1 = labels(idx1);
labels2 = labels(idx2);
eeg1 = eeg(idx1(1) : idx1(end) * 1250);
eeg2 = eeg(idx2(1) * 1250 : idx2(end) * 1250);

% psd ---------------------------------------------------------------------

[psd1, faxis, ~] = psd_states('eeg', eeg1, 'emg', [],...
    'labels', labels1, 'fs', 1250, 'graphics', true, 'sstates', sstates);
[psd2, ~, ~] = psd_states('eeg', eeg2, 'emg', [],...
    'labels', labels2, 'fs', 1250, 'graphics', true, 'sstates', sstates);

% graphics ----------------------------------------------------------------
xLimit = [0 100];
yLimit_norm = [10^-4 10^-1];
yLimit = [10^1 10^4];

fh = figure;
for istate = 1 : length(sstates)
    subplot(2, 3, istate)
    ph = plot(faxis, psd1(istate, :), 'k', 'LineWidth', 1);
    hold on
    ph = plot(faxis, psd2(istate, :), 'b', 'LineWidth', 1);
    title(cfg_names{sstates(istate)})
    xlim(xLimit)
    ylim(yLimit)
    xlabel('Frequency [Hz]')
    ylabel('PSD [mV^2/Hz]')
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    legend('Before', 'After')
    
    subplot(2, 3, istate + length(sstates))
    ph = plot(faxis, psd1(istate, :) ./ sum(psd1(istate, :), 2), 'k', 'LineWidth', 1);
    hold on
    ph = plot(faxis, psd2(istate, :) ./ sum(psd2(istate, :), 2), 'b', 'LineWidth', 1);
    xlim(xLimit)
    ylim(yLimit_norm)
    xlabel('Frequency [Hz]')
    ylabel('norm PSD')
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')      
end

% arrange to prism --------------------------------------------------------
clear prismData
count = 1;
for istate = 1 : 2 : length(sstates) * 2
    prismData(istate, :) = psd1(count, :);
    prismData(istate + 1, :) = psd2(count, :);
    count = count + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare two different sessions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rec1 --------------------------------------------------------------------
basepath = 'K:\Data\lh95\lh95_210825_080400';
session = CE_sessionTemplate(basepath, 'viaGUI', false,...
    'force', true, 'saveVar', true);      
nchans = session.extracellular.nChannels;
[~, basename] = fileparts(basepath);
cd(basepath)
load([basename, '.AccuSleep_labels.mat'])
infoname = dir('*datInfo*');
load(infoname(end).name)
if isfield(datInfo, 'fs')
    fs = datInfo.fs;
else
    fs = 20000;
end

start = 0;
dur = floor(datInfo.nsamps(1) / fs);
labelsTemp1 = labels(1 : dur);

EEG1 = double(bz_LoadBinary([basename, '.lfp'], 'duration', dur,...
    'frequency', 1250, 'nchannels', nchans, 'start', start,...
    'channels', [1 : 4], 'downsample', 1));
EEG1 = mean(EEG1, 2);    
[psd1, faxis, ~] = psd_states('eeg', EEG1, 'emg', [],...
    'labels', labelsTemp1, 'fs', 1250, 'graphics', true, 'sstates', sstates);

% rec2 --------------------------------------------------------------------
basepath = 'K:\Data\lh95\lh95_210825_080400';
session = CE_sessionTemplate(basepath, 'viaGUI', false,...
    'force', true, 'saveVar', true);      
nchans = session.extracellular.nChannels;
[~, basename] = fileparts(basepath);
cd(basepath)
load([basename, '.AccuSleep_labels.mat'])
infoname = dir('*datInfo*');
load(infoname(end).name)
if isfield(datInfo, 'fs')
    fs = datInfo.fs;
else
    fs = 20000;
end

start = ceil(datInfo.nsamps(1) / fs);
dur = Inf;
labelsTemp2 = labels(start : end);

EEG2 = double(bz_LoadBinary([basename, '.lfp'], 'duration', length(labels),...
    'frequency', 1250, 'nchannels', nchans, 'start', start,...
    'channels', [1 : 4], 'downsample', 1));
EEG2 = mean(EEG2, 2);

[psd2, ~, ~] = psd_states('eeg', EEG2, 'emg', [],...
    'labels', labelsTemp2, 'fs', fs, 'graphics', true, 'sstates', sstates);

% graphics ----------------------------------------------------------------
xLimit = [0 100];
yLimit_norm = [10^-4 10^-1];
yLimit = [10^1 10^5];

fh = figure;
sb1 = subplot(2, 2, 1);
ph = plot(faxis, psd1, 'LineWidth', 3);
xlabel('Frequency [Hz]')
ylabel('PSD [mV^2/Hz]')
set(ph, {'color'}, cfg_colors(sstates))
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlim(xLimit)
ylim(yLimit)
title('TDT (lh95)')

sb2 = subplot(2, 2, 2);
ph = plot(faxis, psd2, 'LineWidth', 3);
xlabel('Frequency [Hz]')
ylabel('PSD [mV^2/Hz]')
set(ph, {'color'}, cfg_colors(sstates))
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlim(xLimit)
ylim(yLimit)
title('OE (lh93)')

sb3 = subplot(2, 2, 3);
ph = plot(faxis, psd1 ./ sum(psd1, 2), 'LineWidth', 3);
xlabel('Frequency [Hz]')
ylabel('norm PSD')
set(ph, {'color'}, cfg_colors(sstates))
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlim(xLimit)
ylim(yLimit_norm)

sb4 = subplot(2, 2, 4);
ph = plot(faxis, psd2 ./ sum(psd2, 2), 'LineWidth', 3);
xlabel('Frequency [Hz]')
ylabel('norm PSD')
set(ph, {'color'}, cfg_colors(sstates))
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
xlim(xLimit)
ylim(yLimit_norm)


