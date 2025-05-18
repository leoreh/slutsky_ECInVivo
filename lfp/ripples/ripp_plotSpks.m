function ripp_plotSpks(ripp, varargin)
% RIPP_PLOTSPKS Generates summary plots for ripple-related spiking activity.
%
% SUMMARY:
% This function takes the output structure from ripp_spks.m (which itself
% takes output from ripp_detect.m) and creates a figure summarizing key
% spiking patterns relative to ripple events. It visualizes:
%   - Multi-unit (MU) spike PETHs (Peri-Event Time Histograms) around ripple
%     peaks and control event centers.
%   - Comparison of MU spike counts during ripples vs. control events.
%   - Single-unit (SU) spike PETHs (normalized) for all units.
%   - Average normalized SU PETHs for putative regular-spiking (RS) and
%     fast-spiking (FS) units.
%   - Comparison of mean firing rates during ripples vs. control events per unit.
%
% INPUT:
%   ripp          	Structure containing ripple detection and spike analysis
%                   results from ripp_detect.m and ripp_spks.m. Must include
%                   fields under ripp.spks (e.g., mu.rippMap, su.rippMap,
%                   su.ctrlRates, etc.) and ripp.maps (e.g., durWin).
%   basepath        (Optional) Path to recording session directory {pwd}.
%   flgSaveFig      (Optional) Logical flag to save the figure {true}.
%
% OUTPUT:
%   None. Generates and optionally saves a figure.
%
% DEPENDENCIES:
%   setMatlabGraphics (custom), PlotColorMap (custom), export_fig (external),
%   Requires units.mat file in basepath for RS/FS separation.
%
% HISTORY:
% Aug 2024 LH - Refactored from plot_rippleSpks.m to align with ripp_* style.
% 11 Jan 23 LH - Original version (plot_rippleSpks).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING & INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addRequired(p, 'ripp', @isstruct); % Require the ripple structure
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgSaveFig', true, @islogical);

parse(p, ripp, varargin{:}); % Pass 'ripp' directly
ripp        = p.Results.ripp; % Retrieve ripp from parsed results
basepath    = p.Results.basepath;
flgSaveFig  = p.Results.flgSaveFig;

% Files
cd(basepath);
[~, basename] = fileparts(basepath);
rippfile = fullfile(basepath, [basename, '.ripp.mat']); % For loading if needed
unitsfile = fullfile(basepath, [basename, '.units.mat']); % For RS/FS classification

% Load required data if 'ripp' is empty or units are missing
if isempty(ripp)
    if exist(rippfile, 'file')
        load(rippfile, 'ripp');
    else
        error('Input ''ripp'' structure is empty and %s not found.', rippfile);
    end
end

if exist(unitsfile, 'file')
    load(unitsfile, 'units'); % Contains 'units' struct with 'clean' field for RS/FS
else
    warning('units.mat file not found at %s. Cannot separate RS/FS units.', unitsfile);
    units = []; % Set units to empty to handle downstream plotting gracefully
end

% Validate necessary fields in ripp structure
if ~isfield(ripp, 'spks') || ~isfield(ripp.spks, 'info') || ~isfield(ripp.spks.info, 'mapDur')
    error('Input ''ripp'' structure is missing required fields');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARATIONS & PARAMETER EXTRACTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nRipples = size(ripp.times, 1);
nUnits = size(ripp.spks.su.rippMap, 1);             % Get nUnits from the map dimension

mapDur = ripp.maps.durWin;                          % Use map window from ripp struct [s]
nbinsMap = ripp.spks.info.nBinsMap;                 % Use bin count from ripp.spks.info
xMap = linspace(mapDur(1), mapDur(2), nbinsMap);    % Time axis for maps [s]

% Determine indices for RS and FS units if 'units' struct loaded successfully
if ~isempty(units) && size(units.clean, 1) >= 2
    flgUnits = true;
    idxRS = units.clean(1, :);
    idxFS = units.clean(2, :);
else
    flgUnits = false;
    idxRS = [];
    idxFS = [];
end

% graphics params
unitClr = [0 0 1; 1 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAPHICS GENERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMatlabGraphics(true); % Apply consistent graphic settings
fh = figure('Name', [basename '_rippleSpks'], 'NumberTitle', 'off');

% --- Multi-Unit (MU) Plots ---
if ~isempty(ripp.spks.mu) 

    % PETH map for a subset of ripples
    sb1 = subplot(4, 3, [1, 4]);
    nRipples_plot = min([100, nRipples]);
    repIdx = randperm(nRipples, nRipples_plot);
    mapData = ripp.spks.mu.rippMap(repIdx, :);
    PlotColorMap(mapData, 'x', xMap, 'bar','on');
    xlabel('Time relative to peak [s]');
    ylabel(sprintf('Ripple #'));
    title('MU Spike Counts (PETH)');

    % Mean MU PETH during ripples
    sb2 = subplot(4, 3, 2); hold on; cla;
    ydata = ripp.spks.mu.rippMap;
    plot_stdShade('axh', sb2, 'dataMat', ydata', 'alpha', 0.3, 'clr', [0 0 0], 'xVal', xMap);
    axis tight;
    yLimit = ylim; % Store y-limit for consistent scaling with control plot
    yLimit(1) = 0;
    ylim(yLimit);
    xlabel('Time s]');
    ylabel('MU Spikes');
    title('Ripple');
    box off;

    % Mean MU PETH during control events
    sb3 = subplot(4, 3, 5); hold on; cla;
    ydata = ripp.spks.mu.ctrlMap;
    plot_stdShade('axh', sb3, 'dataMat', ydata', 'alpha', 0.3, 'clr', [0.5 0.5 0.5], 'xVal', xMap);
    ylim(yLimit) % Apply same y-limit
    xlim([xMap(1), xMap(end)])
    xlabel('Time [s]');
    ylabel('MU Spikes');
    title('Control');
    box off;

    % Histogram of total MU spikes per event (ripple vs. control)
    sb4 = subplot(4, 3, [3, 6]); hold on;
    rippleCnts = sum(ripp.spks.mu.rippMap, 2);
    histogram(rippleCnts, 50, 'Normalization', 'probability', ...
        'DisplayStyle', 'stairs', 'EdgeColor', 'k', 'LineWidth', 1.5);
    ctrlCnts = sum(ripp.spks.mu.ctrlMap, 2);
    histogram(ctrlCnts, 50, 'Normalization', 'probability', ...
        'DisplayStyle', 'stairs', 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 1.5); % Gray for control
    legend({'Ripple', 'Control'}, 'Location', 'best');
    ylabel('Probability');
    xlabel('MU Spikes per Event');
    title('MU Spike Count Distribution');
    box off;
end

% --- Single-Unit (SU) Plots ---
if ~isempty(ripp.spks.su) 

    % Normalized SU PETH map (mean across ripples, normalized per unit)
    sb5 = subplot(4, 3, [7, 10]);
    ydata = squeeze(mean(ripp.spks.su.rippMap, 2, 'omitnan'));  % [nUnits x nBins]
    maxPerUnit = max(ydata, [], 2);
    maxPerUnit(maxPerUnit == 0) = 1;                            % Avoid division by zero for silent units
    normSpking = ydata ./ maxPerUnit;
    PlotColorMap(normSpking, 'x', xMap, 'bar','on');
    xlabel('Time[s]');
    ylabel('Unit #');
    title('Norm. Mean SU PETH (Ripple)');

    % Mean normalized SU PETH for RS units
    sb6 = subplot(4, 3, 8); hold on; cla
    if flgUnits
        ydata = normSpking(idxRS, :);
        plot_stdShade('axh', sb6, 'dataMat', ydata', 'alpha', 0.3, 'clr', unitClr(1, :), 'xVal', xMap);
        title(sprintf('RS Units (n=%d)', sum(idxRS)));
        ylim([0 1]); 
    end
    xlabel('Time[s]');
    ylabel('Norm. FR');
    box off; axis tight;

    % Mean normalized SU PETH for FS units
    sb7 = subplot(4, 3, 11); hold on;
    if flgUnits
        ydata = normSpking(idxFS, :);
        plot_stdShade('axh', sb7, 'dataMat', ydata', 'alpha', 0.3, 'clr', unitClr(2, :), 'xVal', xMap);
        title(sprintf('FS Units (n=%d)', sum(idxFS)));
        ylim([0 1]); 
    end
    xlabel('Time [s]');
    ylabel('Norm. FR');
    box off; axis tight;

    % Scatter plot of mean ripple rate vs. mean control rate per unit
    sb8 = subplot(4, 3, [9, 12]); cla;
    rippFR = mean(ripp.spks.su.rippRates, 2, 'omitnan');
    ctrlFR = mean(ripp.spks.su.ctrlRates, 2, 'omitnan');

    if flgUnits
        plot(ctrlFR(idxRS), rippFR(idxRS), '.b', 'MarkerSize', 10); hold on;
        plot(ctrlFR(idxFS), rippFR(idxFS), '.r', 'MarkerSize', 10);
        legend({'RS', 'FS'}, 'Location', 'best');
    else
        plot(ctrlFR, rippFR, '.k', 'MarkerSize', 10); hold on;
    end

    set(gca, 'yscale', 'log', 'xscale', 'log');
    maxVal = max([ctrlFR; rippFR], [], 'all', 'omitnan');
    minVal = min(nonzeros([ctrlFR; rippFR]), [], 'all', 'omitnan');
    eqLim = [minVal, maxVal];
    eqLim = [0.01, 100];
    plot(eqLim, eqLim, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlim(eqLim); ylim(eqLim);
    xlabel('MFR in Control [Hz]');
    ylabel('MFR in Ripple [Hz]');
    title('Unit Firing Rates');
    box off;
end

% --- Overall Title ---
sgtitle([basename ' - Ripple Spike Analysis'], 'Interpreter', 'none');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flgSaveFig
    figpath = fullfile(basepath, 'graphics');
    if ~exist(figpath, 'dir')
        mkdir(figpath);
    end
    figname = fullfile(figpath, sprintf('%s_ripp_plotSpks', basename)); % Use new function name
    saveas(fh, [figname '.png']); % Fallback to MATLAB's saveas
end

end % END FUNCTION ripp_plotSpks

% EOF