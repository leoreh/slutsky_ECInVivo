basepath = 'I:\Data\Processed\lh58\fepsp\lh58_200916_1024';

cd(basepath)
session = CE_sessionTemplate(pwd, 'viaGUI', false,...
    'force', true, 'saveVar', true);
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;
[~, basename] = fileparts(basepath);

filename = fullfile(basepath, [basename, '.dat']);

data = double(bz_LoadBinary(filename, 'nChannels', nchans, 'channels', [25 : 32],...
    'frequency', fs, 'start', 62, 'duration', 0.5));


close all
fh = figure;

idx = round(0.157 * fs : 0.17 * fs);
b2uv = 0.195;

subplot(1, 2, 1)
hold on
count = 1;
plot(data(idx, 5) * b2uv)
count = count + 1;
plot([0 5 * fs / 1000], [0 0], 'k')
plot([0 0], [0, -b2uv * 500], 'k')
set(gca, 'YTick', [], 'XTick', [], 'XColor', 'w', 'YColor', 'w')
box off
xlabel('5 ms')
ylabel('0.5 mV')

idx = 0.3 * fs : 0.45 * fs;
subplot(1, 2, 2)
hold on
count = 1;
plot(-data(idx, 4) * b2uv)
count = count + 1;
plot([0 50 * fs / 1000], [0 0], 'k')
plot([0 0], [0, -b2uv * 1000], 'k')
box off
set(gca, 'YTick', [], 'XTick', [], 'XColor', 'w', 'YColor', 'w')
xlabel('50 ms')
ylabel('1 mV')




