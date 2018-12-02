function plotCluster(spikes, varargin)

% INPUT
% required:
%   spikes      struct in bz format including recduration
% optional:
%   graphics    plot figure {1}.
%   win         time window for calculation {all}.
%   binsize     size in s of bins {60}.
%   save        save figure {1}.
%   basePath    recording session path {pwd}
%
% OUTPUT
% spkcount      struct with fields strd, norm, avg, std, spkcount.bins, spkcount.binsize
%
% 24 nov 18 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'savefig', 1, @islogic);
addOptional(p, 'filepath', pwd);

parse(p,varargin{:})
savefig = p.Results.savefig;
filepath = p.Results.filepath;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grpcolor = ['k', 'b', 'r', 'm', 'g', 'y'];

[nchans nsamps] = size(spikes.avgWaveform{1});
nunits = length(spikes.UID);
ngrps = length(unique(spikes.shankID));

plotidx = ceil(sqrt(nunits));
plotidx = [plotidx ceil(nunits * 2 / plotidx)];

wvidx = 1 : 2 : nunits * 2;
histidx = 2 : 2 : nunits * 2;

f = figure;
for i = 1 : nunits
    % waveform
    subplot(plotidx(1), plotidx(2), wvidx(i))
    unitcolor = grpcolor(spikes.shankID(i));
    for j = 1 : nchans
        offset = j * 100;
        errbounds = [spikes.avgWaveform{i}(j, :) + spikes.stdWaveform{i}(j, :);...
            spikes.avgWaveform{i}(j, :) - spikes.stdWaveform{i}(j, :)];
        p = patch([1 : nsamps, nsamps : -1 : 1],...
            [errbounds(1, :), errbounds(2, end : -1 : 1)], unitcolor);
        p.EdgeColor = 'none';
        p.FaceAlpha = 0.5;
        hold on
        l = plot(spikes.avgWaveform{i}(j, :), 'lineWidth', 1, 'Color', unitcolor);
        set(p, 'YData', get(p, 'YData') - offset);
        set(l, 'YData', get(l, 'YData') - offset);
    end
    axis off
    title(int2str(i))
    
    % ISI correlogram
    subplot(plotidx(1), plotidx(2), histidx(i))
    binsize = 0.001;
    bins = [-0.1 : binsize : 0.1] ;
    h = histogram([diff(spikes.times{i}) diff(spikes.times{i}) * -1], bins);
    h.EdgeColor = 'none';
    h.FaceColor = unitcolor;
    axis off
    
end

if savefig
    filename = 'ClusterGraphics';
    saveVectorFig(filename, filepath, f)    
end

end

% EOF


