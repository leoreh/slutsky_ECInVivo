function plot_rippleSpks(ripp, varargin)

% plots spiking activity in realtion to ripples. see getRipples and
% getRippleSpks.m
%
% INPUT:
%   ripp          	struct
%   basepath        path to recording {pwd}
%   saveFig         logical
%
% 11 jan 23 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'saveFig', true, @islogical);

parse(p, varargin{:})
basepath    = p.Results.basepath;
saveFig     = p.Results.saveFig;

% files
cd(basepath)
[~, basename] = fileparts(basepath);
rippfile = fullfile(basepath, [basename, '.ripp.mat']);
unitsfile = fullfile(basepath, [basename, '.units.mat']);

% load
if isempty(ripp)
    if exist(rippfile, 'file')
        load(rippfile, 'ripp')
    else
        error('ripp file missing')
    end
end

if exist(unitsfile, 'file')
    load(unitsfile, 'units')
else
    error('units file missing')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMatlabGraphics(true)
fh = figure;

nepochs = size(ripp.epochs, 1);
durPlot = [-50 50] / 1000;
x = durPlot(1) : diff(durPlot) / nepochs : durPlot(2);

if isfield(ripp.spks, 'mu')

    % map of nspks per ripple
    sb1 = subplot(4, 3, [1, 4]);
    repIdx = randperm(nepochs, min([300, nepochs]));
    PlotColorMap(ripp.spks.mu.rippMap(repIdx, :), 1, 'bar','on', 'x', x);
    xlabel('Time [ms]')
    ylabel('Ripple No.')
    title(sprintf('Representative MU\nspike rate'))

    % mean nspks across ripples
    sb2 = subplot(4, 3, 2);
    ydata = ripp.spks.mu.rippMap;
    xdata = linspace(ripp.maps.durWin(1), ripp.maps.durWin(2),...
        size(ydata, 2));
    plot(xdata, mean(ydata), 'k')
    hold on
    patch([xdata, flip(xdata)], [mean(ydata) + std(ydata),...
        flip(mean(ydata) - std(ydata))],...
        'k', 'EdgeColor', 'none', 'FaceAlpha', .2, 'HitTest', 'off')
    axis tight
    yLimit = ylim;
    xlabel('Time [ms]')
    ylabel('MU spike rate')
    title(sprintf('MU mean spike rate\nduring ripples'))
    
    % mean nspks across random non-ripple epochs
    sb3 = subplot(4, 3, 5);
    ydata = ripp.spks.mu.randMap;
    xdata = linspace(ripp.maps.durWin(1), ripp.maps.durWin(2),...
        size(ydata, 2));
    plot(xdata, mean(ydata), 'g')
    hold on
    patch([xdata, flip(xdata)], [mean(ydata) + std(ydata),...
        flip(mean(ydata) - std(ydata))],...
        'g', 'EdgeColor', 'none', 'FaceAlpha', .2, 'HitTest', 'off')
    axis tight
    ylim(yLimit)
    xlabel('Time [ms]')
    ylabel('MU spike rate')
    title(sprintf('MU mean spike rate\nduring non-ripple epochs'))

    % hist of nspks in ripples vs. random epochs
    sb4 = subplot(4, 3, [3, 6]);
    ydata = sum(ripp.spks.mu.rippMap, 2);
    hh = histogram(ydata, 50,...
        'Normalization', 'probability');
    hh.EdgeColor = 'none';
    hh.FaceColor = 'k';
    hh.FaceAlpha = 0.3;
    hold on
    ydata = sum(ripp.spks.mu.randMap, 2);
    hh = histogram(ydata, 50,...
        'Normalization', 'probability');
    hh.EdgeColor = 'none';
    hh.FaceColor = 'g';
    hh.FaceAlpha = 0.3;
    legend({'Ripple', 'Non-Ripple'}, 'Location', 'best')
    ylabel('Probability')
    xlabel('MU spikes')
    title(sprintf('MU spikes distribution\nper ripple / non-ripple epochs'))

end

if isfield(ripp.spks, 'su')
    
    % map of mean rate per unit across ripples
    sb5 = subplot(4, 3, [7, 10]);
    ydata = squeeze(mean(ripp.spks.su.rippMap, 2));
    ydata = [ydata' ./ max(ydata')]';
    PlotColorMap(ydata, 1, 'bar','on', 'x', x);
    xlabel('Time [ms]')
    ylabel('Unit no.')
    title('Norm. SU MFR');

    % mean nspks across RS units and ripples
    sb6 = subplot(4, 3, 8);
    unitIdx = units.clean(1, :);
    ydata = squeeze(mean(ripp.spks.su.rippMap(unitIdx, :, :), 2));
    xdata = linspace(ripp.maps.durWin(1), ripp.maps.durWin(2),...
        size(ydata, 2));
    ydata = [ydata' ./ max(ydata')]';
    plot(xdata, mean(ydata), 'b')
    hold on
    patch([xdata, flip(xdata)], [mean(ydata) + std(ydata),...
        flip(mean(ydata) - std(ydata))],...
        'b', 'EdgeColor', 'none', 'FaceAlpha', .2, 'HitTest', 'off')
    xlabel('Time [ms]')
    ylabel('Norm. SU MFR')
    legend(sprintf('n = %d', sum(unitIdx)))
    ylim([0 1])
    title(sprintf('Norm. MFR during ripples\naveraged across RS units'));
    axis tight

    % mean nspks across FS units and ripples
    sb7 = subplot(4, 3, 11);
    unitIdx = units.clean(2, :);
    ydata = squeeze(mean(ripp.spks.su.rippMap(unitIdx, :, :), 2));
    xdata = linspace(ripp.maps.durWin(1), ripp.maps.durWin(2),...
        size(ydata, 2));
    ydata = [ydata' ./ max(ydata')]';
    plot(xdata, mean(ydata), 'r')
    hold on
    patch([xdata, flip(xdata)], [mean(ydata) + std(ydata),...
        flip(mean(ydata) - std(ydata))],...
        'b', 'EdgeColor', 'none', 'FaceAlpha', .2, 'HitTest', 'off')
    xlabel('Time [ms]')
    ylabel('Norm. SU MFR')
    legend(sprintf('n = %d', sum(unitIdx)))
    ylim([0 1])
    title(sprintf('Norm. MFR during ripples\naveraged across FS units'));
    axis tight

    % plot of nspks in ripples vs. random epochs per unit
    sb8 = subplot(2, 3, 6);
    unitIdx = units.clean(1, :);
    ydata = ripp.spks.su.rippMap(unitIdx, :, :);
    ydata = squeeze(mean(mean(ydata, 2), 3));
    xdata = ripp.spks.su.randMap(unitIdx, :, :);
    xdata = squeeze(mean(mean(xdata, 2), 3));
    plot(xdata, ydata, '.b', 'MarkerSize', 10)
    hold on
    unitIdx = units.clean(2, :);
    ydata = ripp.spks.su.rippMap(unitIdx, :, :);
    ydata = squeeze(mean(mean(ydata, 2), 3));
    xdata = ripp.spks.su.randMap(unitIdx, :, :);
    xdata = squeeze(mean(mean(xdata, 2), 3));
    plot(xdata, ydata, '.r', 'MarkerSize', 10)
    set(gca, 'yscale', 'log', 'xscale', 'log')
    eqLim = [min([ylim, xlim]), max([ylim, xlim])];
    plot(eqLim, eqLim, '--k', 'LineWidth', 1)
    xlim(eqLim)
    ylim(eqLim)
    ylabel('MFR during ripples')
    xlabel('MFR during non-ripple epochs')
    legend({'RS', 'FS'})
    title(sprintf('MFR per unit during\nripples / non-ripple epochs'))

end

sgtitle(basename)

% save figure
if saveFig
    figpath = fullfile(basepath, 'graphics');
    mkdir(figpath)
    figname = fullfile(figpath, sprintf('%s_rippleSpks', basename));
    export_fig(figname, '-png', '-transparent', '-r300')
end

end

% EOF