function ripp_plotSpks(rippSpks, spkPeth, varargin)
% RIPP_PLOTSPKS Visualizes ripple-modulated spiking activity.
%
%   RIPP_PLOTSPKS(rippSpks, spkPeth, varargin)
%
%   SUMMARY:
%       Generates a comprehensive summary figure:
%       Row 1 (MUA):
%           - PETH Heatmap (Events x Time).
%           - Mean Rates (Ripple vs Control) with SEM shading.
%           - Spike Count Probability Distribution.
%       Row 2 (SUA):
%           - Normalized PETH Heatmap (Units x Time) sorted by peak latency.
%           - Population Mean Response (Split by Cell Type if available).
%           - Modulation Scatter Plot (Baseline FR vs Ripple FR).
%
%   INPUTS:
%       rippSpks    - (Struct) Stats from ripp_spks.m (.frRipp, .frRand, etc).
%       spkPeth     - (Struct) Maps from ripp_spkPeth.m (.ripp, .ctrl).
%       varargin    - Parameter/Value pairs:
%           'basepath'   - (Char) Save location.
%           'flgSaveFig' - (Log)  Save generated figure? (Default: true).
%
%   OUTPUTS:
%       None. Generates and saves a figure.
%
%   DEPENDENCIES:
%       PlotColorMap, plot_stdShade.
%
%   HISTORY:
%       Updated: 23 Jan 2026
%

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'rippSpks', @isstruct);
addRequired(p, 'spkPeth', @isstruct);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgSaveFig', true, @islogical);
parse(p, rippSpks, spkPeth, varargin{:});

basepath = p.Results.basepath;
flgSaveFig = p.Results.flgSaveFig;

%% ========================================================================
%  SETUP
%  ========================================================================

cd(basepath);
[~, basename] = fileparts(basepath);

% Load Unit Classification (if available)
unitsFile = fullfile(basepath, [basename, '.units.mat']);
flgUnits = false;
idxRS = [];
idxFS = [];

try
    load(unitsFile, 'units');
    if isfield(units, 'type')
        flgUnits = true;
        idxRS = units.type == 'RS';
        idxFS = units.type == 'FS';
        fprintf('Loaded unit types: % d RS, %d FS.\n', sum(idxRS), sum(idxFS));
    end
catch ME
    warning('Failed to load units file: %s');
end


%% ========================================================================
%  PREP FIG
%  ========================================================================

% Colors
clrRipp = [0 0 0];          % Black
clrCtrl = [0.5 0.5 0.5];    % Gray
clrRS   = [0 0 1];          % Blue
clrFS   = [1 0 0];          % Red

hFig = figure('Name', [basename '_rippSpks'], 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.1 0.1 0.6 0.8]);

tl = tiledlayout(hFig, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
title(tl, [basename ' - Ripple Modulation'], 'Interpreter', 'none');

%% ========================================================================
%  MUA (Row 1)
%  ========================================================================

muMaps = spkPeth.mu;
muStats = rippSpks.mu;
tVec = muMaps.tstamps;

% PETH Map (Ripples)
% ---------------------
nexttile;
% Subsample if too many events to visualize clearly
nEvents = size(muMaps.ripp, 2);
maxPlot = 200;
if nEvents > maxPlot
    showIdx = sort(randperm(nEvents, maxPlot));
else
    showIdx = 1:nEvents;
end

% Data: [1 x nEvents x nBins] -> squeeze -> [nEvents x nBins]
mapData = squeeze(muMaps.ripp(1, showIdx, :));

PlotColorMap(mapData, 'x', tVec, 'bar', 'on');
title(sprintf('MUA PETH (n=%d)', nEvents));
ylabel('Event #');
xlabel('Time (s)');

% Mean Rates (Ripple vs Control)
% ---------------------------------
nexttile; hold on;

% Ripple (Mean +/- SEM)
% Transpose to [nBins x nEvents] for plot_stdShade validity if needed,
% but plot_stdShade expects [Units/Events x Time].
% Data: [1 x nEvents x nBins] -> squeeze -> [nEvents x nBins]
rData = squeeze(muMaps.ripp);
cData = squeeze(muMaps.ctrl);

plot_stdShade('hAx', gca, 'dataMat', rData, 'xVal', tVec, ...
    'clr', clrRipp, 'alpha', 0.2);

plot_stdShade('hAx', gca, 'dataMat', cData, 'xVal', tVec, ...
    'clr', clrCtrl, 'alpha', 0.2);

axis tight;
xlabel('Time (s)');
ylabel('Spikes / Bin');
title('MUA Mean Response');
legend({'Ripple', 'Control'}, 'Location', 'best', 'Box', 'off');

% Spiking Probability Distribution
% -----------------------------------
nexttile; hold on;

% Sum spikes per event (integrating over time window)
rCounts = sum(rData, 2);
cCounts = sum(cData, 2);

histogram(rCounts, 'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', 'EdgeColor', clrRipp, 'LineWidth', 2);

histogram(cCounts, 'Normalization', 'probability', ...
    'DisplayStyle', 'stairs', 'EdgeColor', clrCtrl, 'LineWidth', 2);

xlabel('Spike Count');
ylabel('Probability');
title('MUA Spike Counts');
axis tight;

%% ========================================================================
%  SUA (Row 2)
%  ========================================================================

suMaps = spkPeth.su;
suStats = rippSpks.su;
tVec = suMaps.tstamps;
nUnits = size(suMaps.ripp, 1);

% Prepare Normalized PETHs
% Average across ripples: [nUnits x nEvents x nBins] -> [nUnits x nBins]
meanPeth = squeeze(mean(suMaps.ripp, 2, 'omitnan'));

% Normalize by Peak
peakRates = max(meanPeth, [], 2);
peakRates(peakRates == 0) = 1; % Avoid div/0
normPeth = meanPeth ./ peakRates;

% Sort by peak latency or magnitude? Let's sort by peak latency.
[~, peakIdx] = max(normPeth, [], 2);
[~, sortOrd] = sort(peakIdx);

% Normalized PETH Map (Sorted Units)
% -------------------------------------
nexttile;
PlotColorMap(normPeth(sortOrd, :), 'x', tVec, 'bar', 'on');
title(sprintf('SU Norm PETH (n=%d)', nUnits));
ylabel('Unit # (Sorted)');
xlabel('Time (s)');

% Population Mean by Type (RS vs FS)
% -------------------------------------
nexttile; hold on;

if flgUnits
    if any(idxRS)
        plot_stdShade('hAx', gca, 'dataMat', normPeth(idxRS, :), 'xVal', tVec, ...
            'clr', clrRS, 'alpha', 0.2);
    end
    if any(idxFS)
        plot_stdShade('hAx', gca, 'dataMat', normPeth(idxFS, :), 'xVal', tVec, ...
            'clr', clrFS, 'alpha', 0.2);
    end
    legend({sprintf('RS (n=%d)', sum(idxRS)), sprintf('FS (n=%d)', sum(idxFS))}, ...
        'Location', 'best', 'Box', 'off');
else
    % Plot all if no types
    plot_stdShade('hAx', gca, 'dataMat', normPeth, 'xVal', tVec, ...
        'clr', 'k', 'alpha', 0.2, 'plotMean', true);
    legend({sprintf('All Units (n=%d)', nUnits)}, 'Location', 'best', 'Box', 'off');
end

axis tight;
xlabel('Time (s)');
ylabel('Norm. Firing Rate');
title('SUA Mean Response');

% Modulation Scatter (Ripple vs Baseline FR)
% ---------------------------------------------
nexttile; hold on;

frRipp = suStats.frRipp;
frRand = suStats.frRand;

% Unity Line
minVal = min([frRipp; frRand], [], 'all', 'omitnan');
maxVal = max([frRipp; frRand], [], 'all', 'omitnan');
% Handle pure zero case for log plot
if minVal <= 0, minVal = 0.01; end

plot([minVal maxVal], [minVal maxVal], 'k--');

if flgUnits
    if any(idxRS)
        scatter(frRand(idxRS), frRipp(idxRS), 20, clrRS, 'filled', ...
            'MarkerFaceAlpha', 0.6);
    end
    if any(idxFS)
        scatter(frRand(idxFS), frRipp(idxFS), 20, clrFS, 'filled', ...
            'MarkerFaceAlpha', 0.6);
    end
else
    scatter(frRand, frRipp, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
end

set(gca, 'XScale', 'log', 'YScale', 'log');
xlim([minVal maxVal]); ylim([minVal maxVal]);
xlabel('Baseline FR (Hz)');
ylabel('Ripple FR (Hz)');
title('Rate Modulation');



%% ========================================================================
%  SAVING
%  ========================================================================

if flgSaveFig
    figDir = fullfile(basepath, 'graphics');
    if ~exist(figDir, 'dir')
        mkdir(figDir);
    end

    figFile = fullfile(figDir, [basename, '_ripp_spks.png']);
    saveas(hFig, figFile);

    figFile = fullfile(figDir, [basename, '_ripp_spks.png']);
    saveas(hFig, figFile);
end
end

