% correct ks for not including a dead time for spike extraction

fs = 20000;
% margin for overlapping spikes. if two spikes are within 250 us they are probably the same spike. 
marg = fs * 0.00025;   
spikes2 = spikes;

for i = 1 : length(spikes.ts)
    
    idx = find(diff(spikes.ts{i}) < marg);
    spikes2.ts{i}(idx + 1) = [];
    spikes2.times{i}(idx + 1) = [];
    spikes2.ids{i}(idx + 1) = [];
    
    overlapping(i) = length(idx);
    
    
    figure
    % ISI histogram - pre
    subplot(1, 2, 1)
    binsize = 0.001;
    bins = [0 : binsize : 0.05];
    h = histogram([diff(spikes.times{i})], bins);
    h.EdgeColor = 'none';
    h.FaceColor = 'k';
    xticks([0 0.05])
    box off
    axis tight
    ax = gca;
    ax.YTick = []; ax.YColor = 'none';
    line([0.003 0.003], ax.YLim, 'color', 'k', 'LineWidth', 1)
    
    % ISI histogram - post
    subplot(1, 2, 2)
    binsize = 0.001;
    bins = [0 : binsize : 0.05];
    h = histogram([diff(spikes2.times{i})], bins);
    h.EdgeColor = 'none';
    h.FaceColor = 'k';
    xticks([0 0.05])
    box off
    axis tight
    ax = gca;
    ax.YTick = []; ax.YColor = 'none';
    line([0.003 0.003], ax.YLim, 'color', 'k', 'LineWidth', 1)
    
end

