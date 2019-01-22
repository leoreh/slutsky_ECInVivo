clu = [2 9 16];
tet = 1;

f = figure;

% feature space
subplot(4, 3, [1 : 6])
plotFet(fet{tet}, 'saveFig', false, 'clu', clu, 'fet', [1 4], 'newFig', false)

% waveform
c = ['g', 'b', 'r'];
for i = 1 : length(clu)
    subplot(4, 3, 6 + i)
    plotWaveform(spikes.avgWaveform{clu(i)}, spikes.stdWaveform{clu(i)}, c(i))
end