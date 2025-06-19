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
%   binSize       - Bin size for firing rate calculation [s] {60}.
%
% OUTPUT:
%   frr           - Structure containing recovery analysis results:
%     .frBsl         - Baseline firing rate [Hz]. NaN for bad units.
%     .frSs          - Post-recovery steady-state firing rate from model [Hz].
%     .frTrough      - Trough firing rate [Hz].
%     .hCapacity     - Homeostatic capacity score (PCA-based composite metric).
%     .nif           - Network Impact Factor [Hz*s]. Leave-one-out analysis
%                      of each neuron's unique contribution to network recovery.
%     .rcvErr        - Recovery error (absolute log2 fold change from baseline).
%     .rcvGain       - Recovery gain (log2 fold change from baseline).
%     .rcvTime       - Time to 90% recovery [s].
%     .rcvSlope      - Initial recovery slope [Hz/min].
%     .normSlope     - Normalized slope (% of recovery range per min).
%     .idxPert       - Index of detected perturbation onset (steepest decline).
%     .idxTrough     - Index of detected recovery onset (trough).
%     .fr            - Smoothed firing rate curves for each neuron.
%     .frMdl      - Fitted model curves for good units.
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
addParameter(p, 'spkThr', 4, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'binSize', 60, @(x) isnumeric(x) && isscalar(x) && x>0);
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

% Set analysis window, limit and threshold spikes
if isinf(winLim(2))
    winLim(2) = max(cellfun(@(x) max(x, [], 'omitnan'), spktimes, 'UniformOutput', true));
end
spktimes = cellfun(@(x) x(x >= winLim(1) & x <= winLim(2)), spktimes, 'UniformOutput', false);
nUnits = length(spktimes);
nSpks = cellfun(@length, spktimes)';
uBad = nSpks < spkThr;     

% Calculate raw firing rates. Units are rows
[frOrig, ~, t] = times2rate(spktimes, 'binsize', binSize, 'c2r', true);

% Denoise FRs for all units to get clean traces
fr = mea_frDenoise(frOrig, t, 'flgPlot', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PERTURBATION DETECTION FROM MFR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detect perturbation onset time. Find the N largest MFR drops (most
% negative derivative) and select the one that occurs at the lowest MFR,
% which is characteristic of a major network state transition.

% Crude outlier detection, only for cleaning the MFR before perturbation
% detection
uOtl = isoutlier(range(fr, 2), 'median', 'ThresholdFactor', 7);
mfr = mean(fr(~uBad | ~uOtl, :), 1, 'omitnan');
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
idxPert = pertWin(bestCandidateIdx);

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
    frFit = mea_frFit(zeros(100, 1), idxPert, 'flgPlot', false);
    frFit = repmat(frFit, nUnits, 1);

    for iUnit = 1:nUnits
        frFit(iUnit) = mea_frFit(fr(iUnit, :), idxPert,...
            'flgPlot', false);
    end
    frFit = catfields(frFit, 'addim', true, [3, 2, 1]);
    frFit.frMdl = squeeze(frFit.frMdl);
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
% Note: frBsl, frTrough, and frSs (steady-state) from the model fit are
% preserved in frFit, as they represent the outcome of the modeling process.
% The frr_metrics function will later calculate purely numerical equivalents.
frFlds = {'frMdl', 'mdlName', 'idxMdl'};

% Extract the specified fields to the top level
for iFld = 1:length(frFlds)
    fldName = frFlds{iFld};
    if isfield(frFit, fldName)
        frr.(fldName) = frFit.(fldName);
    end
end
frr.frFit = frFit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PER-UNIT RECOVERY METRICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: although we previously extracted model-based frBsl, frSs, etc.
% into the frr struct, we now calculate all metrics numerically from the
% fitted curves for consistency. The original model parameters are retained
% within the frr.frFit structure.

% Calculate model-based metrics
frr.mdl = frr_metrics(frr.frMdl, t, idxPert, uGood, binSize, true);

% Calculates recovery metrics without relying on model fits, using the
% smoothed firing rate data directly.
frr.mdlF = frr_metrics(fr, t, idxPert, uGood, binSize, false);

% Hardcoded exclusion of units that weren't perturbed. Confirmed by manual
% inspection.
uBad = uBad | ~frr.mdlF.uPert;
uGood(uBad) = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORGANIZE OUTPUT & SAVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Packages results and generates outputs:
% 1. Store all metrics in output structure
% 2. Save to .frr.mat file if requested
% 3. Generate summary plots if requested

% Calculate MFR from good units and fit model
mfr = mean(fr(uGood, :), 1, 'omitnan');
mfrFit = mea_frFit(mfr, idxPert, 'flgPlot', false);

% Add additional data and metrics to the frr structure
frr.fr = fr;
frr.frOrig = frOrig;
frr.t = t;
frr.uGood = uGood;
frr.uBad = uBad;
frr.mfr = mfr;
frr.mfrFit = mfrFit;
frr = orderfields(frr);

frr.info.basename = basename;
frr.info.runtime = datetime("now");
frr.info.idxPert = idxPert;
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

% Plotting
if flgPlot
    frr_plot(frr);
end

end     % EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION: FRR_METRICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function metrics = frr_metrics(fr, t, idxPert, uGood, binSize, flgMdl)
%FRR_METRICS Calculates firing rate recovery metrics.
%
% This helper function calculates a set of recovery metrics from firing rate
% data, which can be either model-based fitted curves or model-free
% smoothed data. It numerically derives baseline, trough, and steady-state
% values from the provided curves.
%
% INPUTS:
%   frCurves - Firing rate curves (nUnits x nTimebins).
%   t        - Time vector for frCurves.
%   idxPert  - Perturbation onset index (scalar).
%   uGood    - Logical vector for good units (nUnits x 1).
%   binSize  - Bin size in seconds.
%
% OUTPUT:
%   metrics  - Structure with all calculated metric fields.

% -------------------------------------------------------------------------
% NUMERICAL FR VALUES
nUnits = size(fr, 1);
winMarg = 5 * 60 / binSize;
winDur = 60 * 60 / binSize;

% FR Baseline: 1 hour window starting 5 min after recording onset, ensuring
% we don't go past perturbation
bslStart = round(winMarg);  
bslEnd = bslStart + round(winDur);  
bslWin = bslStart:min(bslEnd, idxPert-1);  
frBsl = mean(fr(:, bslWin), 2, 'omitnan');

% FR Steady state: 1 hour window ending 5 min before recording end, ensuring
% we don't go before start
ssEnd = length(t) - round(winMarg);  
ssStart = ssEnd - round(winDur);  
ssWin = max(1, ssStart):ssEnd;  % Ensure 
frSs = mean(fr(:, ssWin), 2, 'omitnan');

% FR Trough: minimum in 20 min window after perturbation
troughWin = idxPert : idxPert + round(15 * 60 / binSize);
[frTrough, idxMin] = min(fr(:, troughWin), [], 2);
idxTrough = idxPert + idxMin - 1;
if ~flgMdl
    troughWin = mode(idxTrough) : mode(idxTrough) + winMarg;
    frTrough = mean(fr(:, troughWin), 2, 'omitnan');
end



% -------------------------------------------------------------------------
% RECOVERY METRICS

% Constanst for numerical stability
c(1) = min(frTrough(frTrough > 0));
c(2) = min(frSs(frSs > 0));
if isempty(c(1)) || isinf(c(1)) || isnan(c(1))
    c(1) = 1e-5;
end
if isempty(c(2)) || isinf(c(2)) || isnan(c(2))
    c(2) = 1e-4;
end
% Hardcode c based on inspection of entire data set
c(1) = 5e-5;
c(2) = 3e-4;
% Hardcode c based on 1 spike in an hour
c(1) = 1 / 3600;
c(2) = c(1);

% Cap FR values used in ratios
frBslC = max([frBsl, repmat(c(2), nUnits, 1)], [], 2);
frTroughC = max([frTrough, repmat(c(1), nUnits, 1)], [], 2);
frSsC = max([frSs, repmat(c(2), nUnits, 1)], [], 2);

% Perturbation Depth
pertDepth = log2(frBslC ./ frTroughC);
if any(pertDepth < 0)
    idxBad = find(pertDepth < 0)';
    warning('Problem with perturbation: %s', num2str(idxBad));
    fprintf('\nBSL FR: %s\n', num2str(frBsl(idxBad)))
end

% Recovery Fidelity
rcvErr = abs(log2(frSsC ./ frBslC));
rcvDiff = abs(frSsC - frBslC);

% Recovery Gain and Work. RcvGain is capped at zero for numerical stability
% in regression models, and because the meaning of further loosing fr is
% not of interest.
rcvGain = log2(frSsC ./ frTroughC);
rcvGain(rcvGain < 0) = 0;
rcvWork = rcvGain ./ pertDepth;

% Successful Recovery
% A unit is considered to have recovered if 1. It was significantly
% perturbed and 2. its recovery is at least thrRcv% of its perturbation depth.
thrPert = 1; 
thrRcv = 0.67;  
uPert = pertDepth >= thrPert;
uRcv = uPert & (frSs >= frTrough + thrRcv * (frBsl - frTrough));
uRcv(~uGood) = false; % Only good units can be considered recovered

% uRcv = uPert & rcvWork > thrRcv;
% frTrgt = frTrough .* (frBsl ./ frTrough) .^ thrRcv


% -------------------------------------------------------------------------
% RECOVERY KINETICS

% Initialize 
bslTime = nan(nUnits, 1);
rcvTime = nan(nUnits, 1);
rcvSlope = nan(nUnits, 1);
normSlope = nan(nUnits, 1);
spkDfct = nan(nUnits, 1);
tPost = t(idxPert:end);

for iUnit = 1:nUnits
    if ~uGood(iUnit)
        continue;
    end
    
    % SPIKE DEFICIT  
    aucExp = trapz(tPost, ones(size(tPost)) * frBsl(iUnit));
    aucObs = trapz(tPost, fr(iUnit, idxPert:end));
    spkDfct(iUnit) = log2(aucExp / (aucObs + 1));

    % KINETICS
    rcvCurve = fr(iUnit, idxTrough(iUnit):end);

    % To baseline (thrRcv% of perturbation depth)
    bslTrgt = frTrough(iUnit) + thrRcv * (frBsl(iUnit) - frTrough(iUnit));
    bslPnt = find(rcvCurve >= bslTrgt, 1, 'first');
    if ~isempty(bslPnt)
        bslTime(iUnit) = (idxTrough(iUnit) + bslPnt - 1 - idxPert) * binSize;
    else
        bslTime(iUnit) = t(end); % Max value if not recovered
    end    

    % To steady state 
    rcvRange = frSs(iUnit) - frTrough(iUnit);
    rcvTrgt = frTrough(iUnit) + thrRcv * rcvRange;
    rcvPnt = find(rcvCurve >= rcvTrgt, 1, 'first');
    if ~isempty(rcvPnt)
        rcvTime(iUnit) = (idxTrough(iUnit) + rcvPnt - 1 - idxPert) * binSize;
    else
        rcvTime(iUnit) = t(end); % Max value if not recovered
    end

    % Find max slope during recovery [Hz/min] and normalize
    if ~isempty(rcvCurve) && length(rcvCurve) > 1
        rcvSlope(iUnit) = (max(diff(rcvCurve)) / binSize) * 60;
        normSlope(iUnit) = (rcvSlope(iUnit) / rcvRange) * 100;
    end
end

% -------------------------------------------------------------------------
% OUTPUT
metrics.frBsl = frBsl;
metrics.frTrough = frTrough;
metrics.frSs = frSs;
metrics.idxTrough = idxTrough;
metrics.pertDepth = pertDepth;
metrics.bslTime = bslTime;         
metrics.rcvTime = rcvTime;
metrics.rcvSlope = rcvSlope;
metrics.normSlope = normSlope;
metrics.rcvGain = rcvGain;
metrics.rcvWork = rcvWork;
metrics.rcvErr = rcvErr;
metrics.rcvDiff = rcvDiff;
metrics.spkDfct = spkDfct;
metrics.uRcv = uRcv;
metrics.uPert = uPert;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION: FRR_PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function frr_plot(frr)
%FRR_PLOT Generates summary plots for firing rate recovery analysis.
%
% INPUTS:
%   frr - The main output structure from mea_frr.

% Prepare table for plotting correlations using model-based metrics
clear varMap
varMap.uGood      = 'uGood';
varMap.frBsl      = 'mdl.frBsl';
varMap.rcvTime    = 'mdl.rcvTime';
varMap.rcvGain    = 'mdl.rcvGain';
varMap.spkDfct    = 'mdl.spkDfct';
varMap.pertDepth  = 'mdl.pertDepth';
tbl = v2tbl('v', frr, 'varMap', varMap);
dataTbl = tbl(tbl.uGood, :);
dataTbl(:, {'uGood', 'UnitID'}) = [];

% Extract data from frr structure for plotting
t = frr.t / 60;  % Convert to minutes once
fr = frr.fr;
unitIdx = find(frr.uGood);

idxPertTime = t(frr.info.idxPert);

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
    plot(hAx, t, frr.frMdl(idxRef, :),...
        'Color', colors(iSmpl,:), 'LineWidth', 2.5,...
        'HandleVisibility', 'off');
end

xline(hAx, idxPertTime, 'k--', {'Pert. Onset'},...
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
plot(hAx, t, frr.mfrFit.frMdl, 'r-', 'LineWidth', 2,...
    'DisplayName', 'MFR (model fit)');
xline(hAx, idxPertTime, 'k--', {'Pert. Onset'},...
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
% Comparison of model and model-free parameters as a sanity check
flgSanity = false;  % Enable this to compare model vs model-free metrics
if flgSanity
    clear varMap
    varMap.uGood    = 'uGood';
    varMap.Bsl      = 'mdl.frBsl';
    varMap.BslMf    = 'mdlF.frBsl';
    varMap.Rcv      = 'mdl.frSs';
    varMap.Ss       = 'mdlF.frSs';
    varMap.RcvTime  = 'mdl.rcvTime';
    varMap.MfTime   = 'mdlF.rcvTime';
    tbl = v2tbl('v', frr, 'varMap', varMap);
    dataTbl = tbl(tbl.uGood, :);
    dataTbl(:, 'uGood') = [];
    dataTbl(:, 'UnitID') = [];

    hFig2 = figure('NumberTitle', 'off', ...
        'Position', [200 200 800 600], 'Color', 'w');
    corrplot(gca, dataTbl, 'Type', 'Spearman', 'TestR', 'on');
    title('Model vs Model-Free Metrics Comparison');
end
end
