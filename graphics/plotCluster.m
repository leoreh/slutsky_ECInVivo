function f = plotCluster(varargin)

% INPUT
%   basepath    path to recording
%   spikes      struct (see getSpikes)
%   clu         units to plot. if == nunits than a figure with all units
%               will be plotted in addition to figures of individual units
%   saveFig     save figure {true} or not (false)
%
% CALLS
%   plotWaveform
%
% 24 nov 18 LH. 
% 04 dec 18. added individual units.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'spikes', []);
addOptional(p, 'basepath', pwd);
addOptional(p, 'clu', []);
addOptional(p, 'saveFig', false, @islogical);

parse(p,varargin{:})
spikes = p.Results.spikes;
clu = p.Results.clu;
basepath = p.Results.basepath;
saveFig = p.Results.saveFig;

if isempty(spikes)
    warning('spikes will be loaded from %s', basepath)
    spikes = getSpikes('basepath', basepath);
end
nunits = length(spikes.UID);
if isempty(clu)
    clu = 1 : nunits;
elseif length(clu) > nunits
    error('specified more units to plot than are available in spikes')
end
if ~isfield(spikes, 'lRat'); spikes.lRat = nan(nunits, 1); end
if ~isfield(spikes, 'iDist'); spikes.iDist = nan(nunits, 1); end
if ~isfield(spikes, 'isi'); spikes.isi = nan(nunits, 1); end
if ~isfield(spikes, 'su'); spikes.su = nan(nunits, 1); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grpcolor = ['k', 'b', 'r', 'm', 'g', 'y'];

% plot individual units
for i = 1 : length(clu)
    f = figure;
    
    % waveform
    subplot(2, 2, 1)
    unitcolor = grpcolor(spikes.shankID(clu(i)));
    plotWaveform(spikes.avgWaveform{clu(i)}, spikes.stdWaveform{clu(i)}, unitcolor, 'vert', spikes.samplingRate)
    
    % ISI histogram
    subplot(2, 2, 2)
    binsize = 0.001;
    bins = [0 : binsize : 0.05];
    h = histogram([diff(spikes.times{clu(i)})], bins);
    h.EdgeColor = 'none';
    h.FaceColor = unitcolor;
    xticks([0 0.05])
    box off
    axis tight
    ax = gca;
    ax.YTick = []; ax.YColor = 'none';
    line([0.003 0.003], ax.YLim, 'color', 'k', 'LineWidth', 1)
    
    % raster
    subplot(2, 2, [3 4])
    y = ones(length(spikes.times{clu(i)}) ,1);
    plot(spikes.times{clu(i)} / 60 / 60, y, '.k', 'markerSize', 0.1)
    axis tight
    axis off
    title(sprintf('nSpks = %d', length(spikes.times{clu(i)})), 'FontSize', 14, 'FontWeight', 'norma');
    
    % discriptives
    if spikes.su(clu(i)) == 1
        su = 'SU';
    elseif spikes.su(clu(i)) == 0
        su = 'MU';
    else
        su = 'NaN';
    end
    txt = sprintf('%d - %s; L = %.2f; iDist = %.2f; ISI = %.2f',...
        clu(i), su, spikes.lRat(clu(i)), spikes.iDist(clu(i)), spikes.isi(clu(i)));
    suptitle(txt)
    
    if saveFig
        fullpath = [basepath, '\graphics'];
        strdate = date;
        figname = fullfile(fullpath, ['clu', int2str(i), '_', strdate]);
        
        if ~exist(fullpath, 'dir')
            sprintf('Creating Graphics folder in %s', basepath)
            mkdir('graphics')
        end
        saveas(f, figname, 'png')
    end
    
end
% close all

% plot all units in a grid
if length(clu) == nunits
    
    plotidx = ceil(sqrt(nunits));
    plotidx = [plotidx ceil(nunits * 2 / plotidx)];
    wvidx = 1 : 2 : nunits * 2;
    histidx = 2 : 2 : nunits * 2;
    
    f = figure;
    for i = 1 : nunits
        
        % waveform
        subplot(plotidx(1), plotidx(2), wvidx(i))
        unitcolor = grpcolor(spikes.shankID(i));
        plotWaveform(spikes.avgWaveform{i}, spikes.stdWaveform{i}, unitcolor)
        title(int2str(spikes.cluID(i)))
        
        % ISI histogram
        subplot(plotidx(1), plotidx(2), histidx(i))
        binsize = 0.0005;
        bins = [-0.15 : binsize : 0.15];
        h = histogram([diff(spikes.times{i}) diff(spikes.times{i}) * -1], bins);
        h.EdgeColor = 'none';
        h.FaceColor = unitcolor;
        axis off
        
    end
    
    if saveFig
        savePdf('clusters', basepath, f)
    end
end

end

% EOF


