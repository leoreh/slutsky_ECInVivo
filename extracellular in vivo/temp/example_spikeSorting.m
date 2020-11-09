% demonstrate spike sorting

b2uv = 0.195;
dur = 0.3;
saveFig = false;

basepath = 'D:\VMs\shared\lh70\lh70_201015_0951';
cd(basepath)
[~, basename] = fileparts(basepath);
filename = fullfile(basepath, [basename '.dat']);

session = CE_sessionTemplate(pwd, 'viaGUI', false,...
    'force', true, 'saveVar', true);
basepath = session.general.basePath;
nchans = session.extracellular.nChannels;
fs = session.extracellular.sr;
spkgrp = session.extracellular.spikeGroups.channels;

sig = bz_LoadBinary(filename, 'start', 25212, 'duration', dur,...
    'nChannels', nchans, 'channels', [12 : 15]);
sig = double(sig) * b2uv;
tstamps = 1 / fs : 1 / fs : dur;

hpsig = filterLFP(sig, 'fs', 20000, 'passband', [800 Inf],...
    'order', 50, 'type', 'fir1', 'dataOnly', true, 'graphics', false);

thr = mean(hpsig) + 5 * std(hpsig);

clear stamps
for i = 1 : size(sig, 2)
    s = find(hpsig(:, i) < -thr(i) | hpsig(:, i) > thr(i));
    s(diff(s) < fs * 0.001) = [];
    stamps{i} = s;
end
stamps = sort(vertcat(stamps{:}));
nspks = length(stamps);
snips = zeros(4, 32, nspks);

for i = 1 : nspks
    spkidx(i, :) = [stamps(i) - 15, stamps(i) + 16];
    snips(:, :, i) = sig(spkidx(i, 1) : spkidx(i, 2), :)';
end

% compute PCA
nFeatures = 4 * 3;
fetMat = zeros(nspks, nFeatures);
if ~isempty(snips)
    for ii = 1 : 4
        [~, pcFeat] = pca(permute(snips(ii, :, :), [3, 2, 1]));
        fetMat(:, ii * 3 - 2 : ii * 3) = (pcFeat(:, 1 : 3));
    end
end
fetMat(:, 13) = double(stamps);
dim = [2, 12];

% select specific spikes
clear spk
spk(1) = 12;
spk(2) = 12;
spk(2) = 14;

close all
fh = figure;
clr = 'bgr';

subplot(3, 3, [1 : 3])
yOffset = max(range(sig));
for i = 1 : size(sig, 2)
    plot(tstamps, sig(:, i) + (yOffset / 2) * i, 'k', 'LineWidth', 1)
    hold on
end
yLimit = ylim;
for i = 1 : length(spk)
    fill([spkidx(spk(i), :) fliplr(spkidx(spk(i), :))]' / fs, sort([yLimit yLimit]), clr(i),...
        'FaceAlpha', 0.2, 'EdgeAlpha', 0)
end
ylabel('Amplitude [uV]')
yticks([])
yticklabels('')

subplot(3, 3, [4 : 6])
yOffset = max(range(hpsig));
for i = 1 : size(sig, 2)
    plot(tstamps, hpsig(:, i) + (yOffset / 1.2) * i, 'k', 'LineWidth', 1)
    hold on
    plot(xlim, -[thr(i) thr(i)] + (yOffset / 1.2) * i, '--r', 'LineWidth', 1)
end
yLimit = ylim;
ylim([yLimit(1), yLimit(2)])
yLimit = ylim;
for i = 1 : length(spk)
   fill([spkidx(spk(i), :) fliplr(spkidx(spk(i), :))]' / fs, sort([yLimit yLimit]), clr(i),...
        'FaceAlpha', 0.2, 'EdgeAlpha', 0)
end
xlabel('Time [s]')
ylabel('Amplitude [uV]')
xticks([0 dur]) 
ylabel('Amplitude [uV]')
yticks([]);
yticklabels('')
set(gca, 'TickLength', [0 0])

for ii = 1 : length(spk)
    subplot(3, 3, 6 + ii)
    for i = 1 : 4
        plot([-15 : 16] / fs * 1000, snips(i, :, spk(ii)) + yOffset * i,...
            clr(ii), 'LineWidth', 2)
        hold on
    end
    axis tight
    yLimit = ylim;
    ylim([yLimit(1) - yOffset, yLimit(2) + yOffset])
    xlabel('Time [ms]')
    ylabel('Amplitude [uV]')
    yticks([])
    yticklabels('')
    xlabel('Time [ms]')
end

subplot(3, 3, 9)
scatter(fetMat(:, dim(1)), fetMat(:, dim(2)), 10, 'k', 'filled')
hold on
for i = 1 : length(spk)
    scatter(fetMat(spk(i), dim(1)), fetMat(spk(i), dim(2)), 30, clr(i), 'filled')
end
xlabel(['E' num2str(ceil(dim(1) / 3)) '; PC' num2str(mod(dim(1), 3))])
ylabel(['E' num2str(ceil(dim(2) / 3)) '; PC1'])
yticks([])
yticklabels('')
xticks([])
xticklabels('')

if saveFig
    figpath = fullfile(basepath, 'graphics');
    mkdir(figpath)
    figname = [figpath '\spikeSorting'];
    export_fig(figname, '-tif', '-r300', '-transparent')
end