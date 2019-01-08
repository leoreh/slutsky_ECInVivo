function spkcount = spkCount(spikes, varargin)

% for each unit calculates firing frequency in Hz, defined as spike count
% per binsize divided by binsize {default 1 min}. Also calculates average
% and std of normailzed firing frequency. Normalization is according to
% average or maximum FR within a specified window.
% 
% INPUT
% required:
%   spikes      struct (see getSpikes)
% optional:
%   graphics    plot figure {1}.
%   win         time window for calculation {[1 Inf]}. specified in min.
%   binsize     size in s of bins {60}.
%   saveFig     save figure {1}.
%   basepath    recording session path {pwd}
%   normMethod  normalize to 'max' or 'avg' FR within normWin {'avg'}.
%   normWin     window to calculate avg or max FR when normalizing {[1 Inf]}.
%               specified in min.
% 
% EXAMPLES:     spkCount(spikes, 'normMethod', 'avg', 'normWin', [90 Inf]);
%               will normalize FR according to the average FR between 90
%               min and the end of the recording.
%
% OUTPUT
% spkcount      struct with fields strd, norm, avg, std, bins, binsize,
%               normMethod, normWin
%
% 24 nov 18 LH. updates:
% 05 jan 18 LH  added normMethod and normWin
% 07 jan 18 LH  added disqualify units and debugging

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
validate_win = @(win) assert(isnumeric(win) && length(win) == 2,...
    'time window must be in the format [start end]');

p = inputParser;
addOptional(p, 'binsize', 60, @isscalar);
addOptional(p, 'graphics', 1, @islogical);
addOptional(p, 'win', [1 Inf], validate_win);
addOptional(p, 'saveFig', 1, @islogical);
addOptional(p, 'basepath', pwd);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'normMethod', 'avg', @ischar);
addOptional(p, 'normWin', [1 Inf], validate_win);
addOptional(p, 'disqualify', 0, @islogical);

parse(p,varargin{:})
graphics = p.Results.graphics;
binsize = p.Results.binsize;
win = p.Results.win;
saveFig = p.Results.saveFig;
basepath = p.Results.basepath;
saveVar = p.Results.saveVar;
normMethod = p.Results.normMethod;
normWin = p.Results.normWin;
disqualify = p.Results.disqualify;

% validate windows
if win(1) == 0; win(1) = 1; end
if win(2) == Inf 
    recDur = floor(max(spikes.spindices(:, 1)) / 60);    % [min]
    win(2) = recDur; 
end
if normWin(1) == 0; normWin(1) = 1; end
if normWin(2) == Inf; normWin(2) = recDur; end
win = win(1) : win(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nunits = length(spikes.UID);
nmints = ceil(win(end) - win(1) + 1);

% calculate spike count
spkcount.strd = zeros(nunits, nmints);
for i = 1 : nunits
    for j = 1 : nmints
        % correct for last minute
        if j > win(end)
            binsize = mod(win(2), 60) * 60;
        end
        spkcount.strd(i, j) = sum(ceil(spikes.times{i} / 60) == win(j)) / binsize;
    end
end

% disqualify units that were mute during the baseline
if disqualify
    j = 1;
    for i = 1 : nunits
        mute = sum(~any(spkcount.strd(i, normWin(1) : normWin(2)), 1));
        if mute > diff(normWin) * 0.5 || ...
                mean(spkcount.strd(i, normWin(1) : normWin(2))) < 0.1
            idx(j) = i;
            j = j + 1;
        end
    end
    spkcount.strd(idx, :) = [];
    nunits = size(spkcount.strd, 1);
end

% normalize spike count
spkcount.norm = zeros(nunits, nmints);
for i = 1 : nunits
    switch normMethod
        case 'max'
            bline = max(spkcount.strd(i, normWin(1) : normWin(2)));
        case 'avg'
            bline = mean(spkcount.strd(i, normWin(1) : normWin(2)));
    end
    spkcount.norm(i, :) = spkcount.strd(i, :) / bline;
end

% calculate mean and std of norm spike count
spkcount.avg = mean(spkcount.norm, 1);
spkcount.std = std(spkcount.norm, 0, 1);
errbounds = [abs(spkcount.avg) + abs(spkcount.std);...
    abs(spkcount.avg) - abs(spkcount.std)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if graphics
    
    f = figure;
    x = ([1 : nmints] / 60);
    
    % raster plot of units
    subplot(3, 1, 1)
    hold on
    for i = 1 : nunits
        y = ones(length(spikes.times{i}) ,1) * i;
        plot(spikes.times{i} / 60 / 60, y, '.k', 'markerSize', 0.1)
    end
    axis tight
    ylabel('Unit #')
    title('Raster Plot')
    
    % firing rate across time
    subplot(3, 1, 2)
    hold on
    for i = 1 : nunits
        plot(x, spkcount.strd(i, :))
    end
    axis tight
    ylabel('Frequency [Hz]')
    title('Spike Count')
    
    % normalized firing rate across time
    subplot(3, 1, 3)
    hold on
    for i = 1 : nunits
        plot(x, spkcount.norm(i, :))
    end
    p = patch([x, x(end : -1 : 1)], [errbounds(1 ,:), errbounds(2, end : -1 : 1)], [.5 .5 .5]);
    p.EdgeColor = 'none';
    plot(x, spkcount.avg, 'lineWidth', 3, 'Color', 'k')
    axis tight
    xlabel('Time [h]')
    ylabel('Norm. Frequency')
    title('Norm. Spike Count')
%     ylim([0 2])
    
    if saveFig
        filename = 'spikeCount';
        savePdf(filename, basepath, f)
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saveVar
    spkcount.win = [win(1) win(end)];
    spkcount.binsize = binsize;
    spkcount.normMethod = normMethod;
    spkcount.normWin = normWin;
    
    [~, filename] = fileparts(basepath);
    save([basepath, '\', filename, '.spkcount.mat'], 'spkcount')
end

end

% EOF