basepath = 'H:\data\TERI\wal3_060918';
[~, filename] = fileparts(basepath);
filename = [basepath '\' filename '.spk.4'];
clu = [5 10 13];
tet = 4;
c = ['g', 'b', 'r'];

f = figure;

% feature space
subplot(3, 5, [1, 2, 6, 7, 11, 12])
plotFet(fet{tet}, 'saveFig', false, 'clu', clu, 'fet', [4 7], 'newFig', false, 'c', c)
ylabel('PC2')
axis tight
ax = gca;
ax.XTick = [];
ax.YTick = [];
title('')
    
% waveform
for i = 1 : length(clu)
    idx = find(spikes.cluID == clu(i) & spikes.shankID == tet);
    spkidx{i} = find(spikes.times{idx} < 216 * 60);
end
wavPre = getWaveforms(filename, 'clu', clu, 'spkidx', spkidx);
for i = 1 : length(clu)
    idx = find(spikes.cluID == clu(i) & spikes.shankID == tet);
    spkidx{i} = find(spikes.times{idx} > 218 * 60);
end
wavPost = getWaveforms(filename, 'clu', clu, 'spkidx', spkidx);

for i = 1 : length(clu)
    avgPreTERI{i} = squeeze(mean(wavPre{i}));
    stdPreTERI{i} = squeeze(std(double(wavPre{i})));
    avgPostTERI{i} = squeeze(mean(wavPost{i}));
    stdPostTERI{i} = squeeze(std(double(wavPost{i})));
end

for i = 1 : length(clu)
    subplot(3, 5, 2 + i)
    plotWaveform(avgPreTERI{i}, stdPreTERI{i}, c(i))
    subplot(3, 5, 7 + i)
    plotWaveform(avgPostTERI{i}, stdPostTERI{i}, c(i))
end

% ISI hist
for i = 1 : length(clu)    
    idx = spikes.cluID == clu(i) & spikes.shankID == tet;
    binsize = 0.001;
    bins = [0 : binsize : 0.05];
    subplot(3, 5, 12 + i)
    h = histogram([diff(spikes.times{idx})], bins);
    h.EdgeColor = 'none';
    h.FaceColor = c(i);
    box off
    axis tight
    ax = gca;
    ax.XTick = [];
    ax.YTick = []; ax.YColor = 'none';
    line([0.003 0.003], ax.YLim, 'color', 'k', 'LineWidth', 1)
    if i == 1
        xticks([0 0.05])
        xlabel('ISI [ms]')
    end
end

filename = fullfile('D:\Google Drive\PhD\Slutsky\Manuscripts\Styr et al 2018\Figures', 'figS9');
orient(f, 'landscape')
print(f, filename, '-dpdf', '-bestfit', '-painters')