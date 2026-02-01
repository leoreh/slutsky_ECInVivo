function recRep_plot(varargin)

% creates a figure where each panel (row) is a different plot across time.
% the relation between panels can be determined such that, for example, the
% buttom panel is a zoom-in view of the top panel.
%
% INPUT:
%   basepath    char. path to session folder {pwd}
%   panels2plot string array. determines which plots to include and in what
%               order. can be "spec", "emg", "hypnogram", "raster", "raw"
%   tlayout     vector of 2 elements descrbigin the tiled layout
%   panelTiles  cell describing the spread of panels across tiles
%   th          handle of tile layout
%   xDur        numeric vector where each element describes the x-range for
%               the corrsponding panel
%   xStart      numeric vector where each element describes time point
%               where the panel starts, ie xlim(1)
%   xRelations  cell of zoom-in relations (see recRep_wrapper)
%   tLine       numeric. time point in seconds to add a dashed line
%   sstates     numeric vector representing states to plot [1, 2, 4, 5]
%   rawCh       numeric vector describing the channels to load for raw data
%   emg         numeric vector of emg data
%   spec        struct (see calc_spec.m)
%   boutTimes   cell array (see as_bouts.m)
%   saveFig     logical
%   adjustT     logical. adjust tstamps accordig to tLine
%   xOffset     numeric. value to add to adjustT
%
% OUTPUT
%
% CALLS
%
% TO DO LIST
%
% 05 mar 24 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'panels2plot', []);
addParameter(p, 'tlayout', []);
addParameter(p, 'panelTiles', []);
addParameter(p, 'th', []);
addParameter(p, 'xDur', []);
addParameter(p, 'xStart', []);
addParameter(p, 'xRelations', []);
addParameter(p, 'tLine', []);
addParameter(p, 'rawCh', []);
addParameter(p, 'emg', []);
addParameter(p, 'spec', []);
addParameter(p, 'boutTimes', []);
addParameter(p, 'saveFig', false, @islogical);
addParameter(p, 'adjustT', false, @islogical);
addParameter(p, 'xOffset', 0);

parse(p, varargin{:})
basepath        = p.Results.basepath;
panels2plot     = p.Results.panels2plot;
tlayout         = p.Results.tlayout;
panelTiles      = p.Results.panelTiles;
th              = p.Results.th;
xDur            = p.Results.xDur;
xStart          = p.Results.xStart;
xRelations      = p.Results.xRelations;
tLine           = p.Results.tLine;
rawCh           = p.Results.rawCh;
emg             = p.Results.emg;
spec            = p.Results.spec;
boutTimes       = p.Results.boutTimes;
saveFig         = p.Results.saveFig;
adjustT         = p.Results.adjustT;
xOffset         = p.Results.xOffset;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load necassary data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(basepath)
[~, basename] = fileparts(basepath);

% load data
vars1 = ["session"; "fr"; "datInfo"];
v = basepaths2vars('basepaths', {basepath}, 'vars', vars1);
v2 = basepaths2vars('basepaths', {basepath}, 'vars', ["sr"]);
if isfield(v2, 'fr'), v.sr = v2.fr; elseif isfield(v2, 'sr'), v.sr = v2.sr; end
if isfield(v, 'SleepState'), [v.ss] = v.SleepState; v = rmfield(v, 'SleepState'); end
if isfield(v, 'sleep_states'), [v.ss] = v.sleep_states; v = rmfield(v, 'sleep_states'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

npanels = length(panels2plot);

% state params
cfg = as_loadConfig();
sstates = [2, 1, 4, 5];

% emg params
yLimit_emg = [prctile(emg, 0.05), prctile(emg, 99.95)];
emg_tstamps = [1 : length(emg)];

% raw params
nchans = v.session.extracellular.nChannels;
fs_raw = v.session.extracellular.sr;

% recording params
recDur = floor(v.session.general.duration);
if fs_raw == 20000
    bit2uv = 0.195;      % oe
else
    bit2uv = 1;          % tdt
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fig params
if isempty(th)
    fh = figure;
    th = tiledlayout(tlayout(1), tlayout(2));
    th.TileSpacing = 'compact';
    th.Padding = 'compact';
end
yLimit = nan(npanels, 2);

% plot
for ipanel = 1 : npanels

    % create panel
    axh(ipanel) = nexttile(th, panelTiles{ipanel}(1), [1, length(panelTiles{ipanel})]);
    hold on

    % set plot duration
    plotLim = [xStart(ipanel) - xDur(ipanel); xStart(ipanel) + xDur(ipanel)]';
    if plotLim(1) < 1
        plotLim(1) = 1;
    end
    if plotLim(2) > recDur
        plotLim(2) = recDur
    end

    % plot
    switch panels2plot(ipanel)

        case 'hypnogram'
            for istate = 1 : length(sstates)
                dummy(istate) = plot(nan, nan, 'color', cfg.colors{sstates(istate)}, 'LineWidth', 25);
                winIdx = InIntervals(boutTimes{sstates(istate)}, plotLim);
                sbouts = boutTimes{sstates(istate)}(winIdx, :);
                if ~isempty(sbouts)
                    plot(sbouts', 1 * ones(size(sbouts))',...
                        'color', cfg.colors{sstates(istate)}, 'LineWidth', 25)
                end
            end
            set(gca, 'ytick', [])
            set(gca, 'YColor', 'none')
            lgdh = legend(dummy, cfg.names(sstates), 'Location', 'north',...
                'AutoUpdate', 'off');
            set(axh(ipanel), 'XTick', [], 'XColor', 'none');

        case 'spec'
            winIdx = InIntervals(spec.tstamps, plotLim);
            specTmp = spec;
            specTmp.s = spec.s(winIdx, :);
            specTmp.tstamps = plotLim(1) : plotLim(2) - 1;
            plot_spec(specTmp, 'ch', 1, 'logfreq', true, 'saveFig', false,...
                'axh', axh(ipanel), 'xtime', 1)
            ylim([0 100])
            yticks([1, 10, 100])
            xlabel([])
            ylabel('Freq. (Hz)')

        case 'emg'
            winIdx = InIntervals(emg_tstamps, plotLim);
            emg_xval = [plotLim(1) : plotLim(2)];
            plot(emg_xval, emg(winIdx) ,'k', 'LineWidth', 0.5)
            ylim(yLimit_emg)
            ylabel('EMG RMS')
            hold on

        case 'raster'
            % arrange units
            [~, unitsIdx] = sortrows(v.units.clean', 'descend');
            nunits = sum(v.units.clean, 2);

            % plot all units in black
            spktimes = cellfun(@(x) [x(InIntervals(x, plotLim))]',...
                v.spikes.times, 'uni', false)';
            spktimes = spktimes(unitsIdx(1 : sum(nunits)));
            LineFormat.Color = [0 0 0];
            lineFormat.LineWidth = 2;
            plotSpikeRaster(spktimes,...
                'PlotType', 'vertline', 'LineFormat', LineFormat);

            % plot fs units in brown
            LineFormat.Color = [0.65 0.16 0.16];
            for iunit = 1 : nunits(1)
                spktimes{iunit} = nan;
            end
            plotSpikeRaster(spktimes,...
                'PlotType', 'vertline', 'LineFormat', LineFormat);

            % axis params
            ylabel('Unit #')
            yLimit(ipanel, :) = ylim;
            yticks([0 : 10 : yLimit(ipanel, end)])
            set(axh(ipanel), 'YDir', 'normal')
            hold on

        case 'raw'
            raw = binary_load([basename, '.dat'], 'duration', xDur(ipanel) * 3,...
                'fs', fs_raw, 'nCh', nchans, 'start', xStart(ipanel) - xDur(ipanel),...
                'ch', rawCh, 'downsample', 1);
            raw_tstamps = [1 : size(raw, 1)] / fs_raw + xStart(ipanel) - xDur(ipanel);

            % ypush = median(prctile(raw, 0.01)) * [0 : size(raw, 2) - 1];
            % hard yLimit
            yLimitRaw = [-2000, 400];
            ypush = (-500) * [0 : size(raw, 2) - 1];
            plot(raw_tstamps, (raw + ypush), 'k')
            ylabel('Voltage (uV)')
            ylim(yLimitRaw)
    end

    % add dashed line at time point of interest
    yLimit(ipanel, :) = ylim;
    if ~isempty(tLine)
        if InIntervals(tLine, plotLim)
            if strcmp(get(axh(ipanel), 'YScale'), 'log')
                yLimit(ipanel, 1) = 1;
            end
            plot([tLine, tLine], yLimit(ipanel, :), '--r', 'LineWidth', 2)
        end
    end

    % determine time units for xticks according to xDur
    if xDur(ipanel) > 60 * 60
        xUnits = 60 * 60;
    elseif xDur(ipanel) > 300
        xUnits = 60;
    elseif xDur(ipanel) > 1
        xUnits = 2;
    else
        xUnits = 0.1;
    end

    % adjust xticks and labels according to real time
    if adjustT
        tickStart = xStart(ipanel) - xDur(ipanel);
        tickEnd = xStart(ipanel) + xDur(ipanel);
        ticksBefore = fliplr(tLine : -xUnits : tickStart);
        ticksAfter = tLine : xUnits : tickEnd;
        newTicks = [ticksBefore, ticksAfter(2 : end)];
        oldTicks = xticks;
        xLabelsVal = (newTicks - tLine) + xOffset;
        xLabels = seconds(xLabelsVal);
        if xUnits >= 30 * 60
            xLabels.Format = 'hh:mm';
        elseif xUnits == 0.1
            xLabels.Format = 'hh:mm:ss.SS';
        else
            xLabels.Format = 'hh:mm:ss';
        end
        xLabels = char(xLabels);
    else
        newTicks = xStart(ipanel) - xDur(ipanel) : xUnits : (xStart(ipanel) + xDur(ipanel));
        xLabels = string(newTicks);
    end
    xticks(newTicks);
    xticklabels(xLabels);

    % set x-limit
    xLimit = [xStart(ipanel), xStart(ipanel) + xDur(ipanel)];
    xlim(xLimit)

end

% x-link panels that have the same duration and start time
linkedPanels = {};
cnt = 1;
identifiers = strings(npanels, 1);
for i = 1:npanels
    identifiers(i) = sprintf('%d_%d', xDur(i), xStart(i));
end
[uIdentifiers, ~, ic] = unique(identifiers, 'stable');
for i = 1:length(uIdentifiers)
    idx = find(ic == i);
    if length(idx) > 1
        linkedPanels{cnt} = idx;
        cnt = cnt + 1;
    end
end
for ix = 1:length(linkedPanels)
    linkaxes(axh(linkedPanels{ix}), 'x');
end

% Enable panning to horizontal only
panObj = pan;
panObj.Motion = 'horizontal';
panObj.Enable = 'on';

% set titles and labels
xlabel(th, 'Time (hh:mm:ss.S)')
title(th, basename, 'Interpreter', 'none')

% add marking of next smallest panel
for ipanel = 1 : npanels
    if ~isempty(xRelations(ipanel))
        for iline = 1 : length(xRelations{ipanel})
            panelIdx = xRelations{ipanel}(iline);
            ph(ipanel) = plot(axh(ipanel), [xStart(panelIdx), xStart(panelIdx) + xDur(panelIdx)],...
                [min(yLimit(ipanel, 1)), min(yLimit(ipanel, 1))], 'b', 'LineWidth', 10);
            ylim(axh(ipanel), yLimit(ipanel, :))
            addlistener(axh(panelIdx), 'XLim', 'PostSet', @(src, evnt) updateLine(axh(panelIdx), evnt, ph(ipanel)));

        end
    end
end

% make sure all panels with raw recordings have the same y-limit
% cellfun(@range, (yLimitRaw), 'uni', false)
% yrange = range(yLimit');
% rawIdx = find(strcmp(panels2plot, 'raw'));
% [~, midx] = max(yrange(rawIdx));
% for iraw = rawIdx
%    yLimit(iraw, :) = yLimit(rawIdx(midx), :);
%    ylim(axh(iraw), yLimit(iraw, :))
% end

% save
if saveFig
    set(gcf,'renderer','Painters')
    figname = fullfile(basepath, 'graphics', [basename, '_rep']);
    export_fig(figname, '-pdf', 'renderer', 'painters')
end

end

% EOF

function updateLine(ax_handle, ~, line_handle)

xLimit = ax_handle.XLim;

% update line position to span the current x-axis limits
set(line_handle, 'XData', xLimit);

end