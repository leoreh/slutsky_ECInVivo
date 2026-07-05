function hFig = evt_plotSpks(evtSpks, spkPeth, varargin)
% EVT_PLOTSPKS Visualizes event-modulated spiking activity.
%
%   hFig = EVT_PLOTSPKS(evtSpks, spkPeth, varargin)
%
%   SUMMARY:
%       Generates a comprehensive summary figure:
%       Row 1 (MUA):
%           - PETH Heatmap (Events x Time).
%           - Mean Rates (Event vs Control) with SEM shading.
%           - Spike Count Probability Distribution.
%       Row 2 (SUA):
%           - Normalized PETH Heatmap (Units x Time) sorted by peak latency.
%           - Population Mean Response (Split by Cell Type if available).
%           - Modulation Scatter Plot (Baseline FR vs Event FR).
%
%   INPUTS:
%       evtSpks    - (Struct) Stats from ripp_spks.m. Flat struct with
%                     per-unit fields (.frEvt, .frCtrl, ...).
%       spkPeth     - (Struct) Maps from ripp_spkPeth.m, grouped as:
%                     .mu / .su, each with .ripp/.ctrl/.tstamps.
%       varargin    - Parameter/Value pairs:
%           'basepath'   - (Char) Save location.
%           'flgSaveFig' - (Log)  Save generated figure? (Default: true).
%
%   OUTPUTS:
%       hFig        - (handle) Figure handle.
%
%   DEPENDENCIES:
%       plot_colorMap, plot_stdShade.
%
%   HISTORY:
%       Jan 2026 - Created.
%       Jul 2026 - Dropped PlotColorMap dependency; fixed evtSpks.mu/.su
%                  field access; hardened maps handling and saving.
%

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'evtSpks', @isstruct);
addRequired(p, 'spkPeth', @isstruct);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgSaveFig', true, @islogical);
addParameter(p, 'name', 'evt', @ischar);
addParameter(p, 'lbl', 'Event', @ischar);
parse(p, evtSpks, spkPeth, varargin{:});

basepath = p.Results.basepath;
flgSaveFig = p.Results.flgSaveFig;
name = p.Results.name;
lbl = p.Results.lbl;

%% ========================================================================
%  VALIDATE
%  ========================================================================
% spkPeth must carry .mu and .su, each with .ripp/.ctrl/.tstamps.
subs = {'mu', 'su'};
flds = {'ripp', 'ctrl', 'tstamps'};
for iSub = 1:numel(subs)
    if ~isfield(spkPeth, subs{iSub})
        error('evt_plotSpks:missingField', 'spkPeth.%s is missing.', subs{iSub});
    end
    for iFld = 1:numel(flds)
        if ~isfield(spkPeth.(subs{iSub}), flds{iFld})
            error('evt_plotSpks:missingField', ...
                'spkPeth.%s.%s is missing.', subs{iSub}, flds{iFld});
        end
    end
end

%% ========================================================================
%  SETUP
%  ========================================================================
[~, basename] = fileparts(basepath);

% Number of single units (rows of the SUA map)
nUnits = size(spkPeth.su.ripp, 1);

% Load unit classification (if available and consistent with nUnits)
unitsFile = fullfile(basepath, [basename, '.units.mat']);
flgUnits = false;
idxRS = false(nUnits, 1);
idxFS = false(nUnits, 1);

if isfile(unitsFile)
    try
        s = load(unitsFile, 'units');
        if isfield(s, 'units') && isfield(s.units, 'type')
            uType = s.units.type(:);
            if numel(uType) == nUnits
                flgUnits = true;
                idxRS = uType == 'RS';
                idxFS = uType == 'FS';
                fprintf('[RIPP]: Loaded unit types: %d RS, %d FS.\n', ...
                    sum(idxRS), sum(idxFS));
            else
                warning('evt_plotSpks:unitMismatch', ...
                    'units.type length (%d) ~= nUnits (%d); skipping type split.', ...
                    numel(uType), nUnits);
            end
        end
    catch ME
        warning('evt_plotSpks:unitsLoad', ...
            'Failed to load units file: %s', ME.message);
    end
end

%% ========================================================================
%  PREP FIG
%  ========================================================================

% Colors
clrEvt = [0 0 0];          % Black
clrCtrl = [0.5 0.5 0.5];    % Gray
clrRS   = [0 0 1];          % Blue
clrFS   = [1 0 0];          % Red

hFig = figure('Name', [basename '_evtSpks'], 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.1 0.1 0.6 0.8]);

tl = tiledlayout(hFig, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
title(tl, [basename ' - ' lbl ' Modulation'], 'Interpreter', 'none');

%% ========================================================================
%  MUA (Row 1)
%  ========================================================================

muMaps = spkPeth.mu;
tVec = muMaps.tstamps(:)';

% Collapse the singleton unit dimension: [1 x nEvents x nBins] -> [nEvents x nBins]
rData = collapseUnit(muMaps.ripp);
cData = collapseUnit(muMaps.ctrl);
nEvents = size(rData, 1);

% PETH Map (Events)
% ---------------------
hAx = nexttile(tl);
if nEvents == 0
    text(0.5, 0.5, 'No events', 'Parent', hAx, 'HorizontalAlignment', 'center');
else
    % Subsample if too many events to visualize clearly
    maxPlot = 200;
    if nEvents > maxPlot
        showIdx = sort(randperm(nEvents, maxPlot));
    else
        showIdx = 1:nEvents;
    end

    plot_colorMap(rData(showIdx, :), 'hAx', hAx, 'x', tVec, 'flgBar', true);
end
title(hAx, sprintf('MUA PETH (n=%d)', nEvents));
ylabel(hAx, 'Event #');
xlabel(hAx, 'Time (s)');

% Mean Rates (Event vs Control)
% ---------------------------------
hAx = nexttile(tl); hold(hAx, 'on');

plot_stdShade('hAx', hAx, 'dataMat', rData, 'xVal', tVec, ...
    'clr', clrEvt, 'alpha', 0.2);
plot_stdShade('hAx', hAx, 'dataMat', cData, 'xVal', tVec, ...
    'clr', clrCtrl, 'alpha', 0.2);

axis(hAx, 'tight');
xlabel(hAx, 'Time (s)');
ylabel(hAx, 'Spikes / Bin');
title(hAx, 'MUA Mean Response');
legend(hAx, {lbl, 'Control'}, 'Location', 'best', 'Box', 'off');

% Spiking Probability Distribution
% -----------------------------------
hAx = nexttile(tl); hold(hAx, 'on');

% Sum spikes per event (integrating over time window)
rCounts = sum(rData, 2);
cCounts = sum(cData, 2);

histogram(hAx, rCounts, 'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', 'EdgeColor', clrEvt, 'LineWidth', 2);
histogram(hAx, cCounts, 'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', 'EdgeColor', clrCtrl, 'LineWidth', 2);

xlabel(hAx, 'Spike Count');
ylabel(hAx, 'Probability');
title(hAx, 'MUA Spike Counts');
axis(hAx, 'tight');

%% ========================================================================
%  SUA (Row 2)
%  ========================================================================

suMaps = spkPeth.su;
tVec = suMaps.tstamps(:)';
nBins = size(suMaps.ripp, 3);

% Mean PETH across events: [nUnits x nEvents x nBins] -> [nUnits x nBins]
% (reshape rather than squeeze to stay valid when nUnits == 1)
meanPeth = reshape(mean(suMaps.ripp, 2, 'omitnan'), nUnits, nBins);

% Normalize by per-unit peak (guard zero / all-NaN units)
peakRates = max(meanPeth, [], 2);
peakRates(~isfinite(peakRates) | peakRates == 0) = 1;
normPeth = meanPeth ./ peakRates;

% Sort by peak latency
[~, peakIdx] = max(normPeth, [], 2);
[~, sortOrd] = sort(peakIdx);

% Normalized PETH Map (Sorted Units)
% -------------------------------------
hAx = nexttile(tl);
if nUnits == 0
    text(0.5, 0.5, 'No units', 'Parent', hAx, 'HorizontalAlignment', 'center');
else
    plot_colorMap(normPeth(sortOrd, :), 'hAx', hAx, 'x', tVec, 'flgBar', true);
end
title(hAx, sprintf('SU Norm PETH (n=%d)', nUnits));
ylabel(hAx, 'Unit # (Sorted)');
xlabel(hAx, 'Time (s)');

% Population Mean by Type (RS vs FS)
% -------------------------------------
hAx = nexttile(tl); hold(hAx, 'on');

if flgUnits
    lgd = {};
    if any(idxRS)
        plot_stdShade('hAx', hAx, 'dataMat', normPeth(idxRS, :), 'xVal', tVec, ...
            'clr', clrRS, 'alpha', 0.2);
        lgd{end+1} = sprintf('RS (n=%d)', sum(idxRS));
    end
    if any(idxFS)
        plot_stdShade('hAx', hAx, 'dataMat', normPeth(idxFS, :), 'xVal', tVec, ...
            'clr', clrFS, 'alpha', 0.2);
        lgd{end+1} = sprintf('FS (n=%d)', sum(idxFS));
    end
    if ~isempty(lgd)
        legend(hAx, lgd, 'Location', 'best', 'Box', 'off');
    end
elseif nUnits > 0
    % Plot all if no types
    plot_stdShade('hAx', hAx, 'dataMat', normPeth, 'xVal', tVec, ...
        'clr', clrEvt, 'alpha', 0.2);
    legend(hAx, {sprintf('All Units (n=%d)', nUnits)}, ...
        'Location', 'best', 'Box', 'off');
end

axis(hAx, 'tight');
xlabel(hAx, 'Time (s)');
ylabel(hAx, 'Norm. Firing Rate');
title(hAx, 'SUA Mean Response');

% Modulation Scatter (Event vs Baseline FR)
% ---------------------------------------------
hAx = nexttile(tl); hold(hAx, 'on');

if isfield(evtSpks, 'frEvt') && isfield(evtSpks, 'frCtrl')
    frEvt = evtSpks.frEvt(:);
    frCtrl = evtSpks.frCtrl(:);

    % Log-axis limits from positive, finite rates only
    posFR = [frEvt; frCtrl];
    posFR = posFR(isfinite(posFR) & posFR > 0);
    if isempty(posFR)
        minVal = 0.01; maxVal = 1;
    else
        minVal = min(posFR);
        maxVal = max(posFR);
    end
    if maxVal <= minVal, maxVal = minVal * 10; end

    % Unity line
    plot(hAx, [minVal maxVal], [minVal maxVal], 'k--');

    % Only split by type when the vectors align with the unit indices
    flgSplit = flgUnits && numel(frEvt) == nUnits;
    if flgSplit
        if any(idxRS)
            scatter(hAx, frCtrl(idxRS), frEvt(idxRS), 20, clrRS, 'filled', ...
                'MarkerFaceAlpha', 0.6);
        end
        if any(idxFS)
            scatter(hAx, frCtrl(idxFS), frEvt(idxFS), 20, clrFS, 'filled', ...
                'MarkerFaceAlpha', 0.6);
        end
    else
        scatter(hAx, frCtrl, frEvt, 20, clrEvt, 'filled', 'MarkerFaceAlpha', 0.6);
    end

    set(hAx, 'XScale', 'log', 'YScale', 'log');
    xlim(hAx, [minVal maxVal]);
    ylim(hAx, [minVal maxVal]);
else
    text(0.5, 0.5, 'No FR stats', 'Parent', hAx, 'HorizontalAlignment', 'center');
end

xlabel(hAx, 'Baseline FR (Hz)');
ylabel(hAx, [lbl ' FR (Hz)']);
title(hAx, 'Rate Modulation');

%% ========================================================================
%  SAVING
%  ========================================================================

if flgSaveFig
    figDir = fullfile(basepath, 'graphics');
    if ~exist(figDir, 'dir')
        mkdir(figDir);
    end
    % Hide interactive axes toolbars so they are not baked into the export
    axAll = findall(hFig, 'Type', 'axes');
    for iAx = 1:numel(axAll)
        if ~isempty(axAll(iAx).Toolbar)
            axAll(iAx).Toolbar.Visible = 'off';
        end
    end

    figFile = fullfile(figDir, [basename, '_', name, '_spks.png']);
    try
        exportgraphics(hFig, figFile, 'Resolution', 300);
    catch
        saveas(hFig, figFile);
    end
end

end     % EOF


%% ========================================================================
%  HELPER: COLLAPSEUNIT
%  ========================================================================

function m = collapseUnit(x)
% COLLAPSEUNIT Reshapes a single-unit map [1 x nEvents x nBins] to
% [nEvents x nBins]. Passes 2-D input through unchanged. Robust to the
% nEvents == 1 case (where SQUEEZE would misorient the result).

if ndims(x) == 3
    m = reshape(x, size(x, 2), size(x, 3));
else
    m = x;
end
end
