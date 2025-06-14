function frr = mea_frr(spktimes, varargin)
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
%   flgForce      - Logical flag to force reanalysis even if frFit.mat exists {false}.
%   winLim        - [start end] time window to analyze [s] {[0 Inf]}.
%   spkThr        - Minimum number of spikes required per unit {300}.
%   binSize       - Bin size for firing rate calculation [s] {30}.
%
% OUTPUT:
%   frr           - Structure containing recovery analysis results:
%     .frBsl         - Baseline firing rate [Hz]. NaN for bad units.
%     .frRcv         - Post-recovery steady-state firing rate from model [Hz].
%     .frTrough      - Trough firing rate [Hz].
%     .hCapacity     - Homeostatic capacity score (PCA-based composite metric).
%     .nif           - Network Impact Factor [Hz*s]. Leave-one-out analysis
%                      of each neuron's unique contribution to network recovery.
%     .rcvErr        - Recovery error (absolute log2 fold change from baseline).
%     .rcvGain       - Recovery gain (log2 fold change from baseline).
%     .rcvTime       - Time to 90% recovery [s].
%     .rcvSlope      - Initial recovery slope [Hz/min].
%     .normSlope     - Normalized slope (% of recovery range per min).
%     .pertOnset     - Index of detected perturbation onset (steepest decline).
%     .rcvOnset      - Index of detected recovery onset (trough).
%     .fr            - Smoothed firing rate curves for each neuron.
%     .fitCurve      - Fitted model curves for good units.
%     .model         - Best-fitting model name for each unit.
%     .pFit          - Fitted parameters for each unit.
%     .rsquare       - R-squared goodness of fit for each unit.
%     .exitflag      - Optimization exit flag for each unit.
%     .goodFit       - Logical flag indicating good fit quality for each unit.
%     .t             - Time vector for rate curves [s].
%     .info          - Analysis parameters and metadata.
%     .frFit         - Complete frFit structure for all units.
%
% DEPENDENCIES:
%   times2rate.m
%   Signal Processing Toolbox (for sgolayfilt)
%   Optimization Toolbox (for lsqcurvefit)
%
% HISTORY:
%   Sep 2024 - Refactored to use robust model-fitting approach.
%   Dec 2024 - Refactored to use frFit output as basis for entire output structure,
%              eliminating duplication and simplifying data organization.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING & INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addRequired(p, 'spktimes', @iscell);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgSave', true, @islogical);
addParameter(p, 'flgPlot', true, @islogical);
addParameter(p, 'flgForce', false, @islogical);
addParameter(p, 'winLim', [0 Inf], @(x) isnumeric(x) && numel(x)==2);
addParameter(p, 'spkThr', 300, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'binSize', 30, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'sgPolyOrder', 3, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'sgFrameSec', 600, @(x) isnumeric(x) && isscalar(x) && x > 0);

parse(p, spktimes, varargin{:});
basepath = p.Results.basepath;
flgSave = p.Results.flgSave;
flgPlot = p.Results.flgPlot;
flgForce = p.Results.flgForce;
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
% Detect perturbation onset time. Find the N largest MFR drops (most
% negative derivative) and select the one that occurs at the lowest MFR,
% which is characteristic of a major network state transition.

% Calculate MFR from the smoothed individual unit traces
mfr = mean(fr(~uBad, :), 1, 'omitnan');
dmfr = diff(mfr);

% Limit the search to a pre-defined window
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

% Check if frFit.mat exists and load if not forcing reanalysis
fitFile = fullfile(basepath, [basename, '.frFit.mat']);
if ~flgForce && exist(fitFile, 'file')
    fprintf('Loading existing frFit from %s\n', fitFile);
    load(fitFile, 'frFit');
else
    % Get frFit structure template by calling mea_frFit with dummy data and
    % pre-allocate the array
    frFit = mea_frFit(zeros(100, 1), pertOnset, 'flgPlot', false);
    frFit = repmat(frFit, nUnits, 1);

    for iUnit = 1:nUnits
        frFit(iUnit) = mea_frFit(fr(iUnit, :), pertOnset,...
            'flgPlot', false);        
    end
    frFit = catfields(frFit, 'addim', true, [3, 2, 1]);
    frFit.fitCurve = squeeze(frFit.fitCurve);
    frFit.mdlName = squeeze(frFit.mdlName);

    % Save frFit if requested
    if flgSave
        save(fitFile, 'frFit', '-v7.3');
        fprintf('Saved frFit to %s\n', fitFile);
    end
end

% Good unit fits
uGood = frFit.goodFit & ~uBad;

% Extract specific fields from frFit to the main frr structure
frFlds = {'frBsl', 'frRcv', 'frTrough', 'fitCurve', 'rcvOnset', 'rcvErr', 'rcvGain', ...
    'pertDepth', 'mdlName', 'idxMdl'};

% Extract the specified fields to the top level
for iFld = 1:length(frFlds)
    fldName = frFlds{iFld};
    if isfield(frFit, fldName)
        frr.(fldName) = frFit.(fldName);
        frFit = rmfield(frFit, fldName);
    end
end
frr.frFit = frFit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PER-UNIT RECOVERY METRICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates recovery metrics for good units:
% 1. Recovery fidelity: steady-state / baseline rate
% 2. Recovery time: time to 90% of baseline
% 3. Recovery slope: initial rate of recovery
% 4. Bad units marked with NaN in output arrays
% Note that although we modeled the overshoot, focusing on kRcv means
% that only the homeostatic recovery is included in the rate estimates.

% Initialize additional metrics arrays with NaNs for all units
rcvTime = nan(nUnits, 1);
rcvSlope = nan(nUnits, 1);
normSlope = nan(nUnits, 1);
uRcv = false(nUnits, 1);

% Loop through good units to calculate recovery Kinetics (Time and Slope).
% Calculated numerically from the full fitted curve for model-agnostic
% results
for iUnit = 1:nUnits
    if ~uGood(iUnit)
        continue
    end

    % Time to threshold recovery
    thrRcv = 0.5;
    rcvRange = frr.frBsl(iUnit) - frr.frTrough(iUnit);
    rcvVal = frr.frTrough(iUnit) + thrRcv * rcvRange;

    % Use individual unit's recovery onset
    rcvCurve = frr.fitCurve(iUnit, frr.rcvOnset(iUnit):end);
    rcvIdx = find(rcvCurve >= rcvVal, 1, 'first');

    % Convert to time relative to perturbation onset. For units that did
    % not recover, use maximum recording length
    if ~isempty(rcvIdx)
        rcvTime(iUnit) = (frr.rcvOnset(iUnit) + rcvIdx - 1 - pertOnset) * binSize;
        uRcv(iUnit) = true;
    else
        rcvTime(iUnit) = t(end); % maximum value + penality (end of recording)
    end
    
    % Find max slope during recovery [Hz/min]
    rcvSlope(iUnit) = (max(diff(rcvCurve)) / binSize) * 60;

    % Normalized slope [% of recovery range per min]
    normSlope(iUnit) = (rcvSlope(iUnit) / rcvRange) * 100;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL-FREE RECOVERY METRICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates recovery metrics without relying on model fits.
% 1. Recovery time is calculated from the smoothed firing rate.
% 2. Spike deficit is the integrated difference between pre-perturbation
%    baseline firing and the actual firing rate post-perturbation.

% Unit baseline FR in window before perturbation (with margin)
mrgnBins = round(600/binSize);     % 5 min
bslWin = 1:max(1, pertOnset - mrgnBins);
mf.frBsl = mean(fr(:, bslWin), 2, 'omitnan');

% Calculate trough FR as minimum in first 30 mins post-perturbation. Use it
% as the recovery onset.
troughWin = pertOnset : pertOnset + round(1800/binSize);
[mf.frTrough, troughIdx] = min(fr(:, troughWin), [], 2);
mf.rcvOnset = pertOnset + troughIdx - 1;

% Calculate "steady state" firing rate as the last 20% of the recording
ssWin = [length(t) - round(length(t) * 0.2) : length(t)];
mf.frSs = mean(fr(:, ssWin), 2, 'omitnan');

% Determine if units recovered based on this metric
mf.uRcv = mf.frSs > mf.frBsl * thrRcv;

% Initialize model-free metrics and post-perturbation time vector
mf.rcvTime = nan(nUnits, 1);
mf.spkDfct = nan(nUnits, 1);
tPost = t(pertOnset:end);

for iUnit = 1:nUnits
    if ~uGood(iUnit)
        continue
    end
        
    % Time to 50% recovery from model-free baseline
    rcvRange = mf.frBsl(iUnit) - mf.frTrough(iUnit);
    if rcvRange > 0
        rcvVal = mf.frTrough(iUnit) + 0.5 * rcvRange;
        
        rcvCurve = fr(iUnit, mf.rcvOnset(iUnit):end);
        rcvIdx = find(rcvCurve >= rcvVal, 1, 'first');
        
        if ~isempty(rcvIdx)
            mf.rcvTime(iUnit) = (mf.rcvOnset(iUnit) + rcvIdx - 1 - pertOnset) * binSize;
        else
            mf.rcvTime(iUnit) = t(end); % Max value if not recovered
        end
    end
    
    % Spike deficit calculation
    % Area under the curve for baseline and actual firing rate
    aucExp = trapz(tPost, ones(size(tPost)) * mf.frBsl(iUnit));
    aucObs = trapz(tPost, fr(iUnit, pertOnset:end));
    % mf.spkDeficit(iUnit) = (aucExp - aucObs) / aucExp;
    mf.spkDfct(iUnit) = log2(aucObs / aucExp);
    if mf.spkDfct(iUnit) == -Inf
        uGood(iUnit) = false;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORGANIZE OUTPUT & SAVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Packages results and generates outputs:
% 1. Store all metrics in output structure
% 2. Save to .frr.mat file if requested
% 3. Generate summary plots if requested

% Calculate MFR from good units and fit model
mfr = mean(fr(~uBad, :), 1, 'omitnan');
mfrFit = mea_frFit(mfr, pertOnset, 'flgPlot', false);

% Add additional data and metrics to the frr structure
frr.fr = fr;
frr.frOrig = frOrig;
frr.t = t;
frr.rcvTime = rcvTime;
frr.rcvSlope = rcvSlope;
frr.normSlope = normSlope;
frr.uRcv = uRcv;
frr.uGood = uGood;
frr.mf = mf;
frr.mfr = mfr;
frr.mfrFit = mfrFit;
frr = orderfields(frr);

frr.info.basename = basename;
frr.info.runtime = datetime("now");
frr.info.pertOnset = pertOnset;
frr.info.nSpks = nSpks;
frr.info.uOutlier = uOtl;
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

% -------------------------------------------------------------------------
% Prepare table for plotting correlations
clear varMap
varMap.uGood      = 'uGood';
varMap.frBsl      = 'frBsl';
varMap.rcvTime    = 'rcvTime';
varMap.rcvGain    = 'rcvGain';
varMap.spkDfct    = 'mf.spkDfct';
varMap.pertDepth  = 'pertDepth';
tbl = v2tbl('v', frr, 'varMap', varMap);
dataTbl = tbl(tbl.uGood, :);
dataTbl(:, {'uGood', 'UnitID'}) = [];

if flgPlot
    
    % Extract data from frr structure for plotting
    t = frr.t / 60;  % Convert to minutes once
    fr = frr.fr;
    unitIdx = find(frr.uGood);

    pertOnsetTime = t(frr.info.pertOnset);
    
    hFig = figure('NumberTitle', 'off', ...
        'Position', [100 100 1400 800], 'Color', 'w');
    txtTtl = sprintf('Firing Rate Recovery - %s', frr.info.basename);

    % --- Example fits ---
    hAx = subplot(2, 2, 1);
    hold(hAx, 'on');
    nSmpl = min(7, length(unitIdx));
    rng(1); % for reproducibility
    uSmpl = randperm(length(unitIdx), nSmpl);
    colors = lines(nSmpl);

    for iSmpl = 1:nSmpl
        iUnit = uSmpl(iSmpl);
        idxRef = unitIdx(iUnit);
        plot(hAx, t, fr(idxRef, :), 'Color',...
            [colors(iSmpl,:), 0.4], 'LineWidth', 1.5,...
            'DisplayName', sprintf('Unit %d', idxRef));
        plot(hAx, t, frr.fitCurve(idxRef, :),...
            'Color', colors(iSmpl,:), 'LineWidth', 2.5,...
            'HandleVisibility', 'off');
    end

    xline(hAx, pertOnsetTime, 'k--', {'Pert. Onset'},...
        'LineWidth', 1.5, 'LabelOrientation', 'horizontal',...
        'HandleVisibility', 'off');
    title(hAx, txtTtl);
    xlabel(hAx, 'Time (min)');
    ylabel(hAx, 'Smoothed FR (Hz)');
    legend(hAx, 'show', 'Location', 'best');
    grid(hAx, 'on');
    box(hAx, 'on');
    xlim(hAx, [t(1), t(end)]);
    
    % --- MFR and Model Fit ---
    hAx = subplot(2, 2, 3);
    hold(hAx, 'on');
    plot(hAx, t, frr.mfr, 'b-', 'LineWidth', 2,...
        'DisplayName', 'MFR (smoothed)');
    plot(hAx, t, frr.mfrFit.fitCurve, 'r-', 'LineWidth', 2,...
        'DisplayName', 'MFR (model fit)');
    xline(hAx, pertOnsetTime, 'k--', {'Pert. Onset'},...
        'LineWidth', 1.5, 'LabelOrientation', 'horizontal',...
        'HandleVisibility', 'off');
    title(hAx, 'Population Mean Firing Rate');
    xlabel(hAx, 'Time (min)');
    ylabel(hAx, 'MFR (Hz)');
    legend(hAx, 'show', 'Location', 'best');
    grid(hAx, 'on');
    box(hAx, 'on');
    xlim(hAx, [t(1), t(end)]);

    % --- Correlations ---
    hAx = subplot(2, 2, [2, 4]);
    hold(hAx, 'on');
    corrplot(hAx, dataTbl, 'Type', 'Spearman', 'TestR', 'on');


    % -------------------------------------------------------------------------
    % Comparison of model and model free parameters as a sanity check
    flgSanity = false;
    if flgSanity
        clear varMap
        varMap.uGood    = 'info.uGood';
        varMap.Bsl      = 'frBsl';
        varMap.BslMf    = 'mf.frBsl';
        varMap.Rcv      = 'frRcv';
        varMap.Ss       = 'mf.frSs';
        varMap.RcvTime  = 'rcvTime';
        varMap.MfTime   = 'mf.rcvTime';
        tbl = v2tbl('v', frr, 'varMap', varMap);
        dataTbl = tbl(tbl.uGood, :);
        dataTbl(:, 'uGood') = [];
        dataTbl(:, 'UnitID') = [];

        hFig = figure;
        corrplot(gca, dataTbl, 'Type', 'Spearman', 'TestR', 'on');
    end

end

end     % EOF
