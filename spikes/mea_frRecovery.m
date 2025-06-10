function frr = mea_frRecovery(spktimes, varargin)
% MEAFRRECOVERY Calculates firing rate recovery parameters for all neurons.
%
% SUMMARY:
% This function analyzes firing rate recovery after a network-wide perturbation.
% It uses a robust workflow to calculate metrics for each neuron:
%   1. Firing rates are calculated and smoothed using a Savitzky-Golay filter
%      to preserve recovery shape features (e.g., overshoots).
%   2. The perturbation time is automatically detected from the smoothed
%      population mean firing rate (MFR).
%   3. Each unit's post-perturbation firing rate is fitted with a
%      double-exponential model to capture both recovery and overshoot dynamics.
%   4. Three key, largely independent metrics are derived from the model fit:
%      a) Baseline Firing Rate: Pre-perturbation average rate from robust linear fit.
%      b) Recovery Fidelity: Ratio of steady-state to baseline rate.
%      c) Recovery Kinetics: Normalized slope and time to 90% recovery,
%         derived from the model's primary recovery time constant.
%
% INPUT (Required):
%   spktimes      - Cell array. spktimes{i} contains spike times (s) for neuron i.
%
% INPUT (Optional Key-Value Pairs):
%   basepath      - Path to recording session directory {pwd}.
%   flgSave       - Logical flag to save results to .frr.mat file {true}.
%   flgPlot       - Logical flag to generate all analysis plots {true}.
%   winLim        - [start end] time window to analyze [s] {[0 Inf]}.
%   spkThr        - Minimum number of spikes required per unit {300}.
%   binSize       - Bin size for firing rate calculation [s] {30}.
%
% OUTPUT:
%   frr           - Structure containing recovery analysis results:
%     .bslFr         - Baseline firing rate [Hz]. NaN for bad units.
%     .hCapacity     - Homeostatic capacity score (PCA-based composite metric).
%     .nif           - Network Impact Factor [Hz*s]. Leave-one-out analysis
%                      of each neuron's unique contribution to network recovery.
%     .recovFr       - Post-recovery steady-state firing rate from model [Hz].
%     .recovError    - Recovery error (absolute log2 fold change from baseline).
%     .recovChange   - Recovery change (log2 fold change from baseline).
%     .recovTime     - Time to 90% recovery [s].
%     .recovSlope    - Initial recovery slope [Hz/min].
%     .normSlope     - Normalized slope (% of recovery range per min).
%     .pertOnset     - Index of detected perturbation onset (steepest decline).
%     .recovOnset    - Index of detected recovery onset (trough).
%     .fr            - Smoothed firing rate curves for each neuron.
%     .frFit         - Fitted model curves for good units.
%     .t             - Time vector for rate curves [s].
%     .info          - Analysis parameters and metadata.
%
% DEPENDENCIES:
%   times2rate.m
%   Signal Processing Toolbox (for sgolayfilt)
%   Optimization Toolbox (for lsqcurvefit)
%
% HISTORY:
%   Sep 2024 - Refactored to use robust model-fitting approach.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING & INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addRequired(p, 'spktimes', @iscell);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgSave', true, @islogical);
addParameter(p, 'flgPlot', true, @islogical);
addParameter(p, 'winLim', [0 Inf], @(x) isnumeric(x) && numel(x)==2);
addParameter(p, 'spkThr', 300, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'binSize', 30, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'sgPolyOrder', 3, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'sgFrameSec', 600, @(x) isnumeric(x) && isscalar(x) && x > 0);

parse(p, spktimes, varargin{:});
basepath = p.Results.basepath;
flgSave = p.Results.flgSave;
flgPlot = p.Results.flgPlot;
winLim = p.Results.winLim;
spkThr = p.Results.spkThr;
binSize = p.Results.binSize;
sgPolyOrder = p.Results.sgPolyOrder;
sgFrameSec = p.Results.sgFrameSec;

cd(basepath);
[~, basename] = fileparts(basepath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIRING RATE CALCULATION & SMOOTHING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes and smooths firing rates for all units:
% 1. Calculate raw rates using fixed-width bins (times2rate)
% 2. Apply Savitzky-Golay filter to preserve recovery features
% 3. Filter: polynomial order=3, frame length=600s

% Set analysis window if not fully specified
if isinf(winLim(2))
    winLim(2) = max(cellfun(@(x) max(x, [], 'omitnan'), spktimes, 'UniformOutput', true));
end
spktimes = cellfun(@(x) x(x >= winLim(1) & x <= winLim(2)), spktimes, 'UniformOutput', false);
nUnits = length(spktimes);

% Calculate raw firing rates. Units are rows
[frOrig, ~, t] = times2rate(spktimes, 'binsize', binSize, 'c2r', false);

% Denoise FRs for all units to get clean traces
fr = mea_frDenoise(frOrig, t);

% Unit selection 
nSpks = cellfun(@length, spktimes)';
uOutlier = isoutlier(range(fr, 2), 'median', 'ThresholdFactor', 10);
uGood = find(nSpks >= spkThr & ~uOutlier);
nGood = length(uGood);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PERTURBATION DETECTION FROM MFR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detects perturbation and recovery times from population activity:
% 1. Calculate mean firing rate (MFR) across all units
% 2. Find perturbation onset as steepest decline in MFR
% 3. Identify recovery onset as the start of a significant upward trend

% Calculate MFR from the smoothed individual unit traces
mfr = mean(fr, 1, 'omitnan');
dmfr = diff(mfr);

% Detect perturbation onset time from the MFR's steepest decline. Include
% buffer from recording start
pertBuffer = 10;
[~, pertOnset] = min(dmfr(pertBuffer : end));

% To robustly detect recovery onset, look for a sustained period
% of MFR increase, defined as the first time the MFR derivative is
% positive for nUp consecutive bins.
nUp = 5;
runsUp = conv(dmfr(pertOnset:end) > 0, ones(1, nUp), 'valid');
longUp = find(runsUp >= nUp, 1, 'first');

if ~isempty(longUp)
    % Found a sustained rise. Onset is the start of this period.
    recovDelay = longUp;
    recovOnset = pertOnset + recovDelay - 1;
else
    % Fallback: No sustained rise. Look for the *first* non-negative
    % derivative, which indicates the activity has stopped declining.
    recovDelay = find(dmfr(pertOnset:end) >= 0, 1, 'first');
    recovOnset = pertOnset + recovDelay - 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KINETIC MODEL FITTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fits recovery model to each unit's post-perturbation data:
% Model: y(t) = frSS - (frSS - frTrough)*exp(-kRecov*t) + aOver*t*exp(-kOver*t)
% Parameters:
%   frSS: steady-state firing rate [Hz]
%   kRecov: primary recovery rate constant [1/s]
%   aOver: overshoot amplitude [Hz/s]
%   kOver: overshoot rate constant [1/s]

% Fit model to each unit
for iUnit = 1:nUnits
    frFit(iUnit) = mea_frFit(fr(iUnit, :), recovOnset, pertOnset);
end
frFit = catfields(frFit, 'addim', true, [3, 2, 1]);
frFit.fitCurve = squeeze(frFit.fitCurve);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PER-UNIT RECOVERY METRICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates recovery metrics for good units:
% 1. Recovery fidelity: steady-state / baseline rate
% 2. Recovery time: time to 90% of baseline
% 3. Recovery slope: initial rate of recovery
% 4. Bad units marked with NaN in output arrays
% Note that although we modeled the overshoot, focusing on kRecov means
% that only the homeostatic recovery is included in the rate estimates.

% Initialize output arrays with NaNs for all units
recovFr = nan(nUnits, 1);
recovError = nan(nUnits, 1);
recovChange    = nan(nUnits, 1);
recovTime     = nan(nUnits, 1);
recovSlope    = nan(nUnits, 1);
normSlope     = nan(nUnits, 1);
bslFr         = nan(nUnits, 1);

% Loop through good units to calculate and store metrics
for iGood = 1:nGood
    idxRef = uGood(iGood);

    if frFit.rsquare(idxRef) > 0
        % Extract fitted parameters for clarity
        pFit = frFit.pFit(idxRef, :);
        frSS = pFit(1);
        kRecov = pFit(2);
        frTrough = frFit.troughFr(idxRef);
        bslFr(idxRef) = frFit.bslFr(idxRef);

        % Store steady-state FR from model
        recovFr(idxRef) = frSS;

        % Recovery Error; fold change relative to to baseline
        recovChange(idxRef) = log2(frSS / bslFr(idxRef));
        recovError(idxRef) = abs(recovChange(idxRef));

        % Metric 2: Recovery Kinetics (Time and Slope)
        if kRecov > eps
            % Time to 90% recovery (based on the primary exponential term)
            recovTime(idxRef) = -log(0.1) / kRecov; % in seconds

            % Normalized slope [% of recovery range per min]
            normSlope(idxRef) = kRecov * 60 * 100;

            % Initial recovery slope [Hz/min]. Note this depends on the
            % petrubation effect size.
            recovSlope(idxRef) = (frSS - frTrough) * kRecov * 60;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IDENTIFY BAD UNITS BASED ON FIT QUALITY AND STABILITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find bad units based on criteria

% Define thresholds
thrRsquare = 0.1;
thrBslSlope = 1;    % Hz/min (equivalent to 45 degrees)
thrPert = 0.2;      % Percent reduction from baseline
thrRecov = 0.1;
thrBsl = 0.001;

tauRecov = 1 ./ frFit.pFit(:, 2);
badFast = tauRecov < binSize;                
badPert = (frFit.bslFr - frFit.troughFr) ./ frFit.bslFr < thrPert;    
badBsl = frFit.bslFr < thrBsl;
badRecov = recovError < thrRecov;    
badSlope = abs(squeeze(frFit.bslFit(:, 2))) > thrBslSlope;
badRsquare = frFit.rsquare < thrRsquare;
badSpks = nSpks < spkThr;

uBad = badFast | badPert | badBsl | badRecov | badSlope | badRsquare | badSpks | uOutlier;
uGood = ~uBad;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HOMEOSTATIC CAPACITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates a composite score that represents a neuron's homeostatic
% capacitiy, a combination of recovery speed and accuracy. 

% Invert recovFidelity so that higher values indicate better recovery
% (lower error = better fidelity)
dataVars = [zscore(-recovError(uGood)), zscore(normSlope(uGood))];

% Run PCA
[loadings, pc, ~] = pca(dataVars);
pc = pc(:, 1);

% Ensure direction is correct (high score = good recovery)
% Both loadings should be positive for higher values to indicate better recovery
if loadings(1, 1) < 0 || loadings(2, 1) < 0
    pc = -pc;
end

hCapacity = nan(nUnits, 1);
hCapacity(uGood) = pc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NETWORK IMPACT FACTOR (NIF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates each neuron's unique contribution to network MFR recovery
% using leave-one-out analysis. NIF answers: "What is the unique,
% indispensable contribution of each neuron to the network's MFR recovery?"
% 1. Calculate baseline network recovery using fitted steady-state values
% 2. For each unit, recalculate network recovery excluding that unit
% 3. NIF = difference in recovery error (baseline - leave-one-out)
% 4. Higher NIF = more critical unit for network recovery

% Calculate baseline network recovery (simple mean across units)
recovMfr = mean(recovFr(uGood));
bslMfr = mean(bslFr(uGood));
mfrError = log2(recovMfr / bslMfr);

% Initialize
nif = nan(nUnits, 1);
looBsl = nan(nUnits, 1);

% Loop through each good unit for leave-one-out analysis
for iUnit = 1:nUnits
    if ~uGood(iUnit)
        continue
    end
    
    % Remove this unit from the analysis
    uLOO = uGood;
    uLOO(iUnit) = false;
    
    % Calculate network recovery without this unit
    recovLOO = mean(recovFr(uLOO));
    looBsl(iUnit) = mean(bslFr(uLOO));
    looError = log2(recovLOO / looBsl(iUnit));
    
    % NIF = difference in recovery error
    nif(iUnit) = mfrError - looError;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORGANIZE OUTPUT & SAVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Packages results and generates outputs:
% 1. Store all metrics in output structure
% 2. Save to .frr.mat file if requested
% 3. Generate summary plots if requested

% Organize output structure
frr.fr = fr;
frr.t = t;
frr.bslFr = bslFr(:);
frr.hCapacity = hCapacity;
frr.nif = nif;
frr.recovFr = recovFr;
frr.looBsl = looBsl;
frr.recovError = recovError;
frr.recovChange = recovChange;
frr.recovTime = recovTime;
frr.recovSlope = recovSlope;
frr.normSlope = normSlope;
frr.frFit = frFit;

frr.info.basename = basename;
frr.info.runtime = datetime("now");
frr.info.pertOnset = pertOnset;
frr.info.recovOnset = recovOnset;
frr.info.nSpks = nSpks;
frr.info.uGood = uGood;
frr.info.winLim = winLim;
frr.info.spkThr = spkThr;
frr.info.binSize = binSize;
frr.info.sgPolyOrder = sgPolyOrder;
frr.info.sgFrameSec = sgFrameSec;

% Save results
if flgSave
    save(fullfile(basepath, [basename, '.frr.mat']), 'frr', '-v7.3');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flgPlot
    plot_frr(frr);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION: PLOT RECOVERY RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_frr(frr, unitIdx)
% PLOT_FRR Generates a 4-panel summary figure of firing rate recovery analysis.
%
% INPUT:
%   frr         - Structure containing recovery analysis results
%   unitIdx     - Optional array of unit indices to plot. If empty or not
%                 provided, uses frr.info.uGood.
%
% DEPENDENCIES:
%   plot_scatterCorr.m

% Extract data from frr structure
t = frr.t;
fr = frr.fr;  % Units are rows
frFit = frr.frFit;
uGood = find(frr.info.uGood);

% Use provided unitIdx if available, otherwise use uGood
if nargin < 2 || isempty(unitIdx)
    unitIdx = uGood;
end

pertOnset = frr.info.pertOnset;
recovOnset = frr.info.recovOnset;
pertOnsetTime = t(frr.info.pertOnset);
recovOnsetTime = t(frr.info.recovOnset);
bslFr = frr.bslFr;
recovError = frr.recovError;
recovTime = frr.recovTime;
recovSlope = frr.recovSlope;
normSlope = frr.normSlope;

nGood = length(unitIdx);

hndFig = figure('NumberTitle', 'off', ...
    'Position', [100 100 1400 800], 'Color', 'w');
txtTtl = sprintf('Firing Rate Recovery - %s', frr.info.basename);

% --- Panel 1: Example fits ---
ax1 = subplot(2, 3, [1,2]);
hold(ax1, 'on');
nSmpl = min(5, nGood);
rng(1); % for reproducibility
smpl_locs = randperm(nGood, nSmpl);
colors = lines(nSmpl);

for i = 1:nSmpl
    i_into_good_list = smpl_locs(i);
    idxRef = unitIdx(i_into_good_list);
    plot(ax1, t/60, fr(idxRef, :), 'Color', [colors(i,:), 0.4], 'LineWidth', 1.5);
    plot(ax1, t/60, frFit.fitCurve(idxRef, :), 'Color', colors(i,:), 'LineWidth', 2.5, 'DisplayName', sprintf('Unit %d', idxRef));
end
xline(ax1, pertOnsetTime/60, 'k--', {'Pert. Onset'}, 'LineWidth', 1.5, 'LabelOrientation', 'horizontal');
xline(ax1, recovOnsetTime/60, 'r--', {'Recov. Onset'}, 'LineWidth', 1.5, 'LabelOrientation', 'horizontal');
title(ax1, txtTtl);
xlabel(ax1, 'Time (min)');
ylabel(ax1, 'Smoothed FR (Hz)');
legend(ax1, 'show', 'Location', 'best');
grid(ax1, 'on');
box(ax1, 'on');
xlim(ax1, [t(1)/60, t(end)/60]);

% --- Panel 2: MFR and Model Fit ---
ax2 = subplot(2, 3, 3);
hold(ax2, 'on');

% Calculate MFR from good units only
mfr_good = mean(fr(unitIdx, :), 1, 'omitnan');

% Fit model directly to MFR
mfr_fit = mea_frFit(mfr_good, recovOnset, pertOnset);

plot(ax2, t/60, mfr_good, 'b-', 'LineWidth', 2, 'DisplayName', 'MFR (smoothed)');
plot(ax2, t/60, mfr_fit.fitCurve, 'r-', 'LineWidth', 2, 'DisplayName', 'MFR (model fit)');
xline(ax2, pertOnsetTime/60, 'k--', {'Pert. Onset'}, 'LineWidth', 1.5, 'LabelOrientation', 'horizontal');
xline(ax2, recovOnsetTime/60, 'r--', {'Recov. Onset'}, 'LineWidth', 1.5, 'LabelOrientation', 'horizontal');
title(ax2, 'Population Mean Firing Rate');
xlabel(ax2, 'Time (min)');
ylabel(ax2, 'MFR (Hz)');
legend(ax2, 'show', 'Location', 'best');
grid(ax2, 'on');
box(ax2, 'on');
xlim(ax2, [t(1)/60, t(end)/60]);

% --- Panel 3: Baseline FR vs Fidelity ---
hndAx = subplot(2, 3, 4);
plot_scatterCorr(bslFr(unitIdx), recovError(unitIdx), ...
    'xLbl', 'Baseline FR (Hz)', ...
    'yLbl', 'Recovery Error (log2 fold change from baseline)', ...
    'hndAx', hndAx);

% --- Panel 4: Baseline FR vs Recovery Time ---
hndAx = subplot(2, 3, 5);
plot_scatterCorr(bslFr(unitIdx), recovTime(unitIdx)/60, ...
    'xLbl', 'Baseline FR (Hz)', ...
    'yLbl', 'Time to 90% Recovery (min)', ...
    'hndAx', hndAx);

% --- Panel 5: Recovery Error vs Recovery Time ---
hndAx = subplot(2, 3, 6);
plot_scatterCorr(recovError(unitIdx), recovTime(unitIdx)/60, ...
    'xLbl', 'Recovery Error (log2 fold change)', ...
    'yLbl', 'Time to 90% Recovery (min)', ...
    'hndAx', hndAx);

end




