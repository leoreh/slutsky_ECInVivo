function ripp_plotSpks(rippSpks, varargin)
% RIPP_PLOTSPKS Generates summary plots for ripple-related spiking activity.
%
% SUMMARY:
% Visualizes MU and SU PETHs, gain, and rate modulation.
%
% INPUT:
%   rippSpks        Structure from ripp_spks.m.
%                   (Fields: .su, .mu, .info.mapDur, .info.nBinsMap)
%   basepath        (Optional) Path to session {pwd}.
%   flgSaveFig      (Optional) Save figure {true}.
%
% OUTPUT:
%   None. Generates figure.
%
% DEPENDENCIES: setMatlabGraphics, PlotColorMap, plot_stdShade.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'rippSpks', @isstruct);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgSaveFig', true, @islogical);
parse(p, rippSpks, varargin{:});

rippSpks = p.Results.rippSpks;
basepath = p.Results.basepath;
flgSaveFig = p.Results.flgSaveFig;

%% ========================================================================
%  INITIALIZATION
%  ========================================================================
cd(basepath);
[~, basename] = fileparts(basepath);
unitsfile = fullfile(basepath, [basename, '.units.mat']);

% Load Units (RS/FS)
if exist(unitsfile, 'file')
    load(unitsfile, 'units');
    flgUnits = true;
    idxRS = units.type == 'RS';
    idxFS = units.type == 'FS';
else
    flgUnits = false;
end

% Params
mapDur = rippSpks.info.mapDur;
nbinsMap = rippSpks.info.nBinsMap;
xMap = linspace(mapDur(1), mapDur(2), nbinsMap);
nUnits = length(rippSpks.su.pVal);
nRipples = size(rippSpks.mu.rippMap, 1); % Estimate directly from map

unitClr = [0 0 1; 1 0 0];

%% ========================================================================
%  PLOTTING
%  ========================================================================
fh = figure('Name', [basename '_rippSpks'], 'NumberTitle', 'off');

% --- MU Plots ---
if ~isempty(rippSpks.mu) && ~isempty(rippSpks.mu.rippMap)

    % PETH Map
    subplot(4, 3, [1, 4]);
    nPlot = min([100, nRipples]);
    if nRipples > 0
        randIdx = randperm(nRipples, nPlot);
        mapData = rippSpks.mu.rippMap(randIdx, :);
        PlotColorMap(mapData, 'x', xMap, 'bar', 'on');
    end
    xlabel('Time relative to peak [s]');
    ylabel('Ripple #');
    title('MU Spike Counts (PETH)');

    % Mean PETH (Ripple)
    sb2 = subplot(4, 3, 2); hold on;
    ydata = rippSpks.mu.rippMap;
    plot_stdShade('hAx', sb2, 'dataMat', ydata', 'alpha', 0.3, 'clr', [0 0 0], 'xVal', xMap);
    axis tight;
    yLimit = ylim; yLimit(1) = 0; ylim(yLimit);
    xlabel('Time [s]'); ylabel('MU Spikes'); title('Ripple'); box off;

    % Mean PETH (Control)
    sb3 = subplot(4, 3, 5); hold on;
    ydata = rippSpks.mu.ctrlMap;
    plot_stdShade('hAx', sb3, 'dataMat', ydata', 'alpha', 0.3, 'clr', [0.5 0.5 0.5], 'xVal', xMap);
    ylim(yLimit); xlim([xMap(1), xMap(end)]);
    xlabel('Time [s]'); ylabel('MU Spikes'); title('Control'); box off;

    % Histogram
    subplot(4, 3, [3, 6]); hold on;
    rCnts = sum(rippSpks.mu.rippMap, 2);
    cCnts = sum(rippSpks.mu.ctrlMap, 2);
    histogram(rCnts, 50, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'LineWidth', 1.5);
    histogram(cCnts, 50, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 1.5);
    legend({'Ripple','Control'},'Location','best');
    ylabel('Probability'); xlabel('Spikes per Event'); title('MU Count Dist.'); box off;
end

% --- SU Plots ---
if ~isempty(rippSpks.su) && ~isempty(rippSpks.su.rippMap)

    % Normalized Map
    subplot(4, 3, [7, 10]);
    ydata = squeeze(mean(rippSpks.su.rippMap, 2, 'omitnan'));
    maxPerUnit = max(ydata, [], 2);
    maxPerUnit(maxPerUnit == 0) = 1;
    normSpking = ydata ./ maxPerUnit;

    PlotColorMap(normSpking, 'x', xMap, 'bar', 'on');
    xlabel('Time [s]'); ylabel('Unit #'); title('Norm. Mean SU PETH');

    % RS Units
    sb6 = subplot(4, 3, 8); hold on;
    if flgUnits && ~isempty(idxRS)
        ydata = normSpking(idxRS, :);
        plot_stdShade('hAx', sb6, 'dataMat', ydata', 'alpha', 0.3, 'clr', unitClr(1,:), 'xVal', xMap);
        title(sprintf('RS Units (n=%d)', sum(idxRS)));
        ylim([0 1]);
    end
    xlabel('Time [s]'); ylabel('Norm. FR'); box off; axis tight;

    % FS Units
    sb7 = subplot(4, 3, 11); hold on;
    if flgUnits && ~isempty(idxFS)
        ydata = normSpking(idxFS, :);
        plot_stdShade('hAx', sb7, 'dataMat', ydata', 'alpha', 0.3, 'clr', unitClr(2,:), 'xVal', xMap);
        title(sprintf('FS Units (n=%d)', sum(idxFS)));
        ylim([0 1]);
    end
    xlabel('Time [s]'); ylabel('Norm. FR'); box off; axis tight;

    % Scatter
    subplot(4, 3, [9, 12]);
    rFR = rippSpks.su.frRipp;
    cFR = rippSpks.su.frRand;

    if flgUnits
        plot(cFR(idxRS), rFR(idxRS), '.b', 'MarkerSize', 10); hold on;
        plot(cFR(idxFS), rFR(idxFS), '.r', 'MarkerSize', 10);
        legend({'RS','FS'},'Location','best');
    else
        plot(cFR, rFR, '.k', 'MarkerSize', 10); hold on;
    end

    set(gca, 'yscale', 'log', 'xscale', 'log');
    eqLim = [0.01 100];
    plot(eqLim, eqLim, '--k');
    xlim(eqLim); ylim(eqLim);
    xlabel('Control Rate [Hz]'); ylabel('Ripple Rate [Hz]');
    title('Firing Rate Modulation'); box off;
end

sgtitle([basename ' - Ripple Spikes'], 'Interpreter', 'none');

%% ========================================================================
%  SAVE
%  ========================================================================
if flgSaveFig
    figDir = fullfile(basepath, 'graphics');
    if ~exist(figDir, 'dir'), mkdir(figDir); end
    saveas(fh, fullfile(figDir, [basename '_ripp_plotSpks.png']));
end

end