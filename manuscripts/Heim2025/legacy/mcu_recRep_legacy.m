
% CREATE A WRAPPER THAT ADDS THE FIRING RATE THROUGHOUT THE EXPERIMENT AND

% creates a figure where each panel (row) is a different plot across time.
% the relation between panels can be determined such that, for example, the
% buttom panel is a zoom-in view of the top panel.

% determine the order and type of panels
panels2plot =["spec", "emg", "hypnogram", "raster", "raw"];
npanels = length(panels2plot);

saveFig = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select file and windows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% file
basepath = 'F:\Data\lh96\lh96_220121_090213';   % wt bac on

% channels for raw trace
rawCh = [1 : 4];

% xDur determines the window size for each plot and xStart determined the
% timepoint where the plot starts. the actual duration that is plotted is
% 3-times xDur to allow for manual scrolling through the data. All window
% params should be given in seconds. Panels that share the same duration
% and start point will be x-linked
xDur(1) = 6 * 60 * 60;
xDur(2) = xDur(1);
xDur(3) = xDur(1);
xDur(4) = 8;
xDur(5) = 0.4;

xStart(1) = 7 * 60 * 60;
xStart(2) = xStart(1);
xStart(3) = xStart(1);
xStart(4) = 7 * 60 * 60 + 60;
xStart(5) = 7 * 60 * 60 + 60;

% set relationships between panels for lines representing zoom in views
xRelations = nan(size(xDur));

% Step 1: Identify unique elements and their last occurrences
[uVals, ~, uIdx] = unique(xDur, 'stable');
lastIdx = accumarray(uIdx(:), 1 : numel(xDur), [], @max);
for idx = 1 : numel(uVals)
    currVal = uVals(idx);
    lastIndex = lastIdx(idx);
    smallerElems = xDur < currVal;
    smallerIndices = find(smallerElems);
    if ~isempty(smallerIndices)
        nextSmallestIndex = smallerIndices(find(smallerIndices > lastIndex, 1));
        if ~isempty(nextSmallestIndex)
            xRelations(lastIndex) = nextSmallestIndex;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(basepath)
[~, basename] = fileparts(basepath);

% load data
varsFile = ["fr"; "spikes"; "datInfo"; "session";...
    "units"; "sleep_states"];
varsName = ["fr"; "spikes"; "datInfo"; "session";...
    "units"; "ss"];
v = getSessionVars('basepaths', {basepath}, 'varsFile', varsFile,...
    'varsName', varsName);

% load raw data
rawIdx = find(contains(panels2plot, 'raw'));
if ~isempty(rawIdx)

    % session params
    nchans = v.session.extracellular.nChannels;
    fs_raw = v.session.extracellular.sr;

    raw = double(bz_LoadBinary([basename, '.dat'], 'duration', xDur(rawIdx) * 3,...
        'frequency', fs_raw, 'nchannels', nchans, 'start', xStart(rawIdx) - xDur(rawIdx),...
        'channels', rawCh, 'downsample', 1));
    raw_tstamps = [1 : size(raw, 1)] / fs_raw + xStart(rawIdx) - xDur(rawIdx);

end

% state params
cfg = as_loadConfig();
sstates = [5, 4, 2, 1];

% sleep signals
filename = fullfile(basepath, [basename, '.sleep_sig.mat']);
load(filename, 'emg');
load(filename, 'eeg');

% re-calc state epochs for better visualization.
% LSLEEP => NREM; N/REM => REM
minDur = [10, 5, 5, 10, 5, 5];
minDur = 4;
interDur = 10;
labels = v.ss.labels;
labels(labels == 6) = 5;
labels(labels == 3) = 4;
[stateEpochs, ~] = as_epochs('labels', labels,...
    'minDur', minDur, 'interDur', interDur);

% subsample emg
fs_emg = 250;
emg = emg(1 : 1250 / fs_emg : end);
emg = (emg - min(emg)) / (max(emg) - min(emg));
yLimit_emg = [prctile(emg, 0.05), prctile(emg, 99.95)];
emg_tstamps = [1 : length(emg)] / fs_emg;

% subsample eeg and re-calc spec
eeg = eeg(1 : 1250 / fs_emg : end);
spec = calc_spec('sig', eeg, 'fs', fs_emg, 'graphics', false, 'saveVar', false,...
    'padfft', -1, 'winstep', 1, 'logfreq', true, 'ftarget', [],...
    'ch', [{1}], 'force', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fig params
fh = figure;
th = tiledlayout(npanels, 1);
th.TileSpacing = 'compact';
% th.Padding = 'none';

% plot
for ipanel = 1 : npanels

    % create panel
    axh(ipanel) = nexttile;
    hold on

    % plot duration
    plotLim = [xStart(ipanel) - xDur(ipanel); xStart(ipanel) + xDur(ipanel)]';

    % plot
    switch panels2plot(ipanel)

        case 'hypnogram'
            fs = 1;
            for istate = 1 : length(sstates)
                dummy(istate) = plot(NaN, NaN, 'color', cfg.colors{sstates(istate)}, 'LineWidth', 25);
                winIdx = InIntervals(stateEpochs{sstates(istate)}, plotLim);
                sepochs = stateEpochs{sstates(istate)}(winIdx, :);
                if ~isempty(sepochs)
                    plot(sepochs', 1 * ones(size(sepochs))',...
                        'color', cfg.colors{sstates(istate)}, 'LineWidth', 25)
                end
            end
            set(gca, 'ytick', [])
            set(gca, 'YColor', 'none')
            lgdh = legend(dummy, cfg.names(sstates), 'Location', 'best',...
                'AutoUpdate', 'off');

        case 'spec'
            fs = 1;
            winIdx = InIntervals(spec.tstamps, plotLim);
            specTmp = spec;
            specTmp.s = spec.s(winIdx, :);
            specTmp.tstamps = plotLim(1) : plotLim(2) - 1;
            plot_spec(specTmp, 'ch', 1, 'logfreq', true, 'saveFig', false,...
                'axh', axh(ipanel), 'xtime', 1)
            ylim([0 20])
            yticks([1, 4, 16])
            xlabel([])
            ylabel('Freq. (Hz)')

        case 'emg'
            fs = fs_emg;
            winIdx = InIntervals(emg_tstamps, plotLim);
            emg_xval = [plotLim(1) : 1 / fs : plotLim(2)];
            plot(emg_xval, emg(winIdx))
            ylim(yLimit_emg)
            ylabel('EMG (a.u.)')

        case 'raster'
            fs = 1;

            % arrange units
            [~, unitsIdx] = sortrows(v.units.clean', 'descend');
            nunits = sum(v.units.clean, 2);

            % plot all units in blue
            spktimes = cellfun(@(x) [x(InIntervals(x, plotLim))]',...
                v.spikes.times, 'uni', false)';
            spktimes = spktimes(unitsIdx(1 : sum(nunits)));
            LineFormat.Color = [0.1 0.1 0.8];
            lineFormat.LineWidth = 2;
            plotSpikeRaster(spktimes,...
                'PlotType', 'vertline', 'LineFormat', LineFormat);

            % plot fs units in red
            LineFormat.Color = [0.8 0.1 0.1];
            for iunit = 1 : nunits(1)
                spktimes{iunit} = nan;
            end
            plotSpikeRaster(spktimes,...
                'PlotType', 'vertline', 'LineFormat', LineFormat);

            % axis params
            ylabel('Unit #')
            yLimit = ylim;
            yticks([0 : 10 : yLimit(end)])
            set(axh(ipanel), 'YDir', 'normal')
            hold on

        case 'raw'
            ypush = median(prctile(raw, 0.1)) * [0 : size(raw, 2) - 1];
            plot(raw_tstamps, raw + ypush, 'k')
            ylabel('Voltage (mV)')

    end

    % adjust xticklabels according to real time
    % determine time units according to xDur
    if xDur(ipanel) > 60 * 60
        xUnits = 60 * 60;
    elseif xDur(ipanel) > 300
        xUnits = 60;
    elseif xDur(ipanel) > 1
        xUnits = 1;
    else
        xUnits = 0.1;
    end

    newTicks = xStart(ipanel) - xDur(ipanel) : xUnits : (xStart(ipanel) + xDur(ipanel));
    xticks(newTicks)
    xLabels = string([xStart(ipanel) - xDur(ipanel) : xUnits : xStart(ipanel) + xDur(ipanel)]);
    xticklabels(xLabels)
    xLimit = [xStart(ipanel), xStart(ipanel) + xDur(ipanel)];
    xlim(xLimit)

end

% x-link panels that have the same duration
uDur = unique(xDur, 'stable');
cnt = 1;
for idur = 1 : length(uDur)
    idx = find(xDur == uDur(idur));
    if length(idx) > 1
        linkedPanels{cnt} = idx;
        cnt = cnt + 1;
    end
end
for ix = 1 : length(linkedPanels)
    linkaxes(axh(linkedPanels{ix}), 'x')
end
% Enable panning to horizontal only
panObj = pan;
panObj.Motion = 'horizontal';
panObj.Enable = 'on';

% set x label
xlabel('Time (s)')

% add marking of next smallest panel
for ipanel = 1 : npanels
    if ~isnan(xRelations(ipanel))
        yLimit = ylim(axh(ipanel));
        panelIdx = xRelations(ipanel);
        ph(ipanel) = plot(axh(ipanel), [xStart(panelIdx), xStart(panelIdx) + xDur(panelIdx)],...
            [0, 0], 'b', 'LineWidth', 5);
        ylim(axh(ipanel), yLimit)
        addlistener(axh(panelIdx), 'XLim', 'PostSet', @(src, evnt) updateLine(axh(panelIdx), evnt, ph(ipanel)));

    end
end

% save
if saveFig
    set(gcf,'renderer','Painters')
    figname = fullfile(basepath, 'graphics', [basename, '_rep']);
    export_fig(figname, '-pdf', 'renderer', 'painters')
end

% EOF

function updateLine(ax_handle, ~, line_handle)
    
    xLimit = ax_handle.XLim;  
    
    % update line position to span the current x-axis limits
    set(line_handle, 'XData', xLimit);

end