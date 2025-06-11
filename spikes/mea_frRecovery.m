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
%     .frBsl         - Baseline firing rate [Hz]. NaN for bad units.
%     .frRecov       - Post-recovery steady-state firing rate from model [Hz].
%     .frTrough      - Trough firing rate [Hz].
%     .hCapacity     - Homeostatic capacity score (PCA-based composite metric).
%     .nif           - Network Impact Factor [Hz*s]. Leave-one-out analysis
%                      of each neuron's unique contribution to network recovery.
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
uOtl = isoutlier(range(fr, 2), 'median', 'ThresholdFactor', 7);
uBad = nSpks < spkThr & uOtl;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PERTURBATION DETECTION FROM MFR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detects perturbation and recovery times from population activity:
% 1. Calculate mean firing rate (MFR) across all units
% 2. Find perturbation onset as steepest decline in MFR
% 3. Identify recovery onset as the start of a significant upward trend

% Calculate MFR from the smoothed individual unit traces
mfr = mean(fr(~uBad, :), 1, 'omitnan');
dmfr = diff(mfr);

% Detect perturbation onset time. Find the N largest MFR drops (most
% negative derivative) and select the one that occurs at the lowest MFR,
% which is characteristic of a major network state transition.
pertWin = round([10, 100] * 60 / binSize);
pertWin = pertWin(1) : pertWin(2);
dmfrWin = dmfr(pertWin);
mfrWin = mfr(pertWin);

% Find the top N most negative derivatives
nCandidates = 10;
[sortedDrops, sortIdx] = sort(dmfrWin, 'ascend');

nTop = min(nCandidates, length(sortedDrops));
topIdx = sortIdx(1:nTop);

% From these candidates, find the one with the lowest MFR value
candidateMFRs = mfrWin(topIdx);
[~, minMfrIdx] = min(candidateMFRs);

% The perturbation onset is the index of this best candidate
bestCandidateIdx = topIdx(minMfrIdx);
pertOnset = pertWin(bestCandidateIdx);

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
    frFit(iUnit) = mea_frFit(fr(iUnit, :), pertOnset,...
        'flgPlot', false, 'binSize', binSize);
end
frFit = catfields(frFit, 'addim', true, [3, 2, 1]);
frFit.fitCurve = squeeze(frFit.fitCurve);

% Good unit fits
uGood = frFit.goodFit & ~uBad;
nGood = sum(uGood);

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
frRecov = nan(nUnits, 1);
recovError = nan(nUnits, 1);
recovChange    = nan(nUnits, 1);
recovTime     = nan(nUnits, 1);
recovSlope    = nan(nUnits, 1);
normSlope     = nan(nUnits, 1);
frBsl         = nan(nUnits, 1);
frTrough      = nan(nUnits, 1);

% Loop through good units to calculate and store metrics
for iUnit = 1:nUnits
    if frFit.exitflag(iUnit) > 0
        % Extract fitted parameters for clarity
        frRecov(iUnit) = frFit.frRecov(iUnit);
        frTrough(iUnit) = frFit.frTrough(iUnit);
        frBsl(iUnit) = frFit.frBsl(iUnit);

        % Recovery Error; fold change relative to to baseline
        if frBsl(iUnit) > 0 && frRecov(iUnit) > 0
            recovChange(iUnit) = log2(frRecov(iUnit) / frBsl(iUnit));
            recovError(iUnit) = abs(recovChange(iUnit));
        end

        % Metric 2: Recovery Kinetics (Time and Slope)
        % Calculated numerically from the full fitted curve for model-agnostic results
        recovRange = frRecov(iUnit) - frTrough(iUnit);
        
        % Time to threshold (50%) recovery (numerically)
        recovThr = 0.5;
        recovVal = frTrough(iUnit) + recovThr * recovRange;
        
        % Use individual unit's recovery onset
        recovOnsetUnit = frFit.recovOnset(iUnit);
        recovCurve = frFit.fitCurve(iUnit, recovOnsetUnit:end);
        recovIdx = find(recovCurve >= recovVal, 1, 'first');

        if ~isempty(recovIdx)
            % Convert to time relative to perturbation onset for comparison between units
            recovTime(iUnit) = (recovOnsetUnit - pertOnset + recovIdx - 1) * binSize; % in seconds
            
            % Find max slope within the window leading up to 50% recovery time.
            slopeWindow = recovCurve(1:recovIdx);
            if length(slopeWindow) > 1
                maxSlopePerBin = max(diff(slopeWindow));
            else
                maxSlopePerBin = 0;
            end
            
            % Maximum recovery slope [Hz/min]
            recovSlope(iUnit) = (maxSlopePerBin / binSize) * 60;
            
            % Normalized slope [% of recovery range per min]
            if recovRange > 0
                normSlope(iUnit) = (maxSlopePerBin / recovRange) / binSize * 60 * 100;
            else
                normSlope(iUnit) = 0;
            end

        else
            recovTime(iUnit) = length(recovCurve) * binSize; % maximum value (end of recording)

            % For units that don't recover, use average slope over the recovery period.
            if recovTime(iUnit) > 0 && recovRange > 0
                avgSlopePerSec = recovRange / recovTime(iUnit);
                recovSlope(iUnit) = avgSlopePerSec * 60; % Convert to Hz/min
                normSlope(iUnit) = (avgSlopePerSec / recovRange) * 60 * 100;
            else
                recovSlope(iUnit) = 0;
                normSlope(iUnit) = 0;
            end
        end

    end
end

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
recovMfr = mean(frRecov(uGood));
bslMfr = mean(frBsl(uGood));
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
    recovLOO = mean(frRecov(uLOO));
    looBsl(iUnit) = mean(frBsl(uLOO));
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
frr.frBsl = frBsl(:);
frr.frRecov = frRecov(:);
frr.frTrough = frTrough(:);
frr.hCapacity = hCapacity;
frr.nif = nif;
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
frr.info.nSpks = nSpks;
frr.info.uGood = uGood;
frr.info.uOutlier = uOtl;
frr.info.winLim = winLim;
frr.info.spkThr = spkThr;
frr.info.binSize = binSize;
frr.info.sgPolyOrder = sgPolyOrder;
frr.info.sgFrameSec = sgFrameSec;
frr.info.qualityChecks = frFit.qualityChecks;

% Save results
if flgSave
    save(fullfile(basepath, [basename, '.frr.mat']), 'frr', '-v7.3');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section generates a 4-panel summary figure of the analysis.

if flgPlot
    
    % Extract data from frr structure for plotting
    t = frr.t;
    fr = frr.fr;
    frFit = frr.frFit;
    unitIdx = find(frr.info.uGood);
    nGood = length(unitIdx);

    pertOnsetTime = t(frr.info.pertOnset);
    frBsl = frr.frBsl;
    recovError = frr.recovError;
    recovTime = frr.recovTime;
    
    hndFig = figure('NumberTitle', 'off', ...
        'Position', [100 100 1400 800], 'Color', 'w');
    txtTtl = sprintf('Firing Rate Recovery - %s', frr.info.basename);

    % --- Panel 1: Example fits ---
    ax1 = subplot(2, 3, [1,2]);
    hold(ax1, 'on');
    nSmpl = min(5, nGood);
    if nSmpl > 0
        rng(1); % for reproducibility
        smpl_locs = randperm(nGood, nSmpl);
        colors = lines(nSmpl);

        for i = 1:nSmpl
            i_into_good_list = smpl_locs(i);
            idxRef = unitIdx(i_into_good_list);
            plot(ax1, t/60, fr(idxRef, :), 'Color', [colors(i,:), 0.4], 'LineWidth', 1.5, 'DisplayName', sprintf('Unit %d', idxRef));
            plot(ax1, t/60, frFit.fitCurve(idxRef, :), 'Color', colors(i,:), 'LineWidth', 2.5, 'HandleVisibility', 'off');
        end
    end
    xline(ax1, pertOnsetTime/60, 'k--', {'Pert. Onset'}, 'LineWidth', 1.5, 'LabelOrientation', 'horizontal', 'HandleVisibility', 'off');
    title(ax1, txtTtl);
    xlabel(ax1, 'Time (min)');
    ylabel(ax1, 'Smoothed FR (Hz)');
    if nSmpl > 0
        legend(ax1, 'show', 'Location', 'best');
    end
    grid(ax1, 'on');
    box(ax1, 'on');
    xlim(ax1, [t(1)/60, t(end)/60]);

    % --- Panel 2: MFR and Model Fit ---
    ax2 = subplot(2, 3, 3);
    hold(ax2, 'on');

    % Calculate MFR from good units only
    mfr_good = mean(fr(unitIdx, :), 1, 'omitnan');

    % Fit model directly to MFR
    mfr_fit = mea_frFit(mfr_good, frr.info.pertOnset);

    plot(ax2, t/60, mfr_good, 'b-', 'LineWidth', 2, 'DisplayName', 'MFR (smoothed)');
    plot(ax2, t/60, mfr_fit.fitCurve, 'r-', 'LineWidth', 2, 'DisplayName', 'MFR (model fit)');
    xline(ax2, pertOnsetTime/60, 'k--', {'Pert. Onset'}, 'LineWidth', 1.5, 'LabelOrientation', 'horizontal', 'HandleVisibility', 'off');
    title(ax2, 'Population Mean Firing Rate');
    xlabel(ax2, 'Time (min)');
    ylabel(ax2, 'MFR (Hz)');
    legend(ax2, 'show', 'Location', 'best');
    grid(ax2, 'on');
    box(ax2, 'on');
    xlim(ax2, [t(1)/60, t(end)/60]);

    % --- Panel 3: Baseline FR vs Fidelity ---
    hndAx = subplot(2, 3, 4);
    plot_scatterCorr(frBsl(unitIdx), recovError(unitIdx), ...
        'xLbl', 'Baseline FR (Hz)', ...
        'yLbl', 'Recovery Error (log2 fold change from baseline)', ...
        'hndAx', hndAx);

    % --- Panel 4: Baseline FR vs Recovery Time ---
    hndAx = subplot(2, 3, 5);
    plot_scatterCorr(frBsl(unitIdx), recovTime(unitIdx)/60, ...
        'xLbl', 'Baseline FR (Hz)', ...
        'yLbl', 'Time to 50% Recovery (min)', ...
        'hndAx', hndAx);

    % --- Panel 5: Recovery Error vs Recovery Time ---
    hndAx = subplot(2, 3, 6);
    plot_scatterCorr(recovError(unitIdx), recovTime(unitIdx)/60, ...
        'xLbl', 'Recovery Error (log2 fold change)', ...
        'yLbl', 'Time to 50% Recovery (min)', ...
        'hndAx', hndAx);
end

end     % EOF
