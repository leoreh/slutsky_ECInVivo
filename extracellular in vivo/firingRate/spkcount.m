function spkcount = spkcount(spikes, varargin)

% INPUT
% required:
%   spikes      struct in bz format including recduration
% optional:
%   graphics    plot figure {1}.
%   win         time window for calculation {[0 Inf]}. specified in seconds.
%   binsize     size in s of bins {60}.
%   save        save figure {1}.
%   basePath    recording session path {pwd}
%
% OUTPUT
% spkcount      struct with fields strd, norm, avg, std, spkcount.bins, spkcount.binsize
%
% 24 nov 18 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments and initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
validate_win = @(win) assert(isnumeric(win) && length(win) == 2,...
    'time window must be in the format [start end]');

p = inputParser;
addOptional(p, 'binsize', 60, @isscalar);
addOptional(p, 'graphics', 1, @islogic);
addOptional(p, 'win', [0 spikes.sessionDur], validate_win);
addOptional(p, 'savefig', 1, @islogic);
addOptional(p, 'filepath', pwd);

parse(p,varargin{:})
graphics = p.Results.graphics;
spkcount.binsize = p.Results.binsize;
win = p.Results.win;
win = win / 60;
savefig = p.Results.savefig;
filepath = p.Results.filepath;

% nbins = ceil(diff(win) / spkcount.binsize);
nunits = length(spikes.UID);
nmints = ceil(win(2)) - win(1);

spkcount.strd = zeros(nunits, nmints);
spkcount.norm = zeros(nunits, nmints);
spkcount.avg = zeros(1, nmints);
spkcount.std = zeros(1, nmints);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate spike count
for i = 1 : nunits
    binsize = 60;
    for j = 1 : nmints
        % correct for last minute
        if j > spikes.sessionDur / 60
            binsize = mod(spikes.sessionDur, 60);
        end
        %         binsize = diff(spkcount.bins(j, :)) + 1);
        %         spkcount.strd(i, j) = sum(spikes.times{i} > spkcount.bins(j) & spikes.times{i} < spkcount.bins(j + 1));
        spkcount.strd(i, j) = sum(ceil(spikes.times{i} / 60) == j) / binsize;
    end
end

% normalize spike count to maximum count
for i = 1 : nunits
    for j = 1 : nmints
        spkcount.norm(i, j) = spkcount.strd(i, j) / max(spkcount.strd(i, :));
    end
end

% calculate average (norm) spike count, std and error bounds
for j = 1 : nmints
    spkcount.avg(j) = mean(spkcount.norm(:, j));
    spkcount.std(j) = std(spkcount.norm(:, j));
end
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
        y = ones(length(spikes.times{i}) ,1) * spikes.UID(i);
        plot(spikes.times{i} / 60 / 60, y, '.k', 'markerSize', 0.1)
    end
    axis tight
    ylabel('Unit #')
    title('Raster Plot')
    
    subplot(3, 1, 2)
    hold on
    for i = 1 : nunits
        plot(x, spkcount.strd(i, :))
    end
    axis tight
    ylabel('Frequency [Hz]')
    title('Spike Count')
    
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
    
end

if savefig
    filename = 'spikeCount';
    saveVectorFig(filename, filepath, f)
end

end

% EOF