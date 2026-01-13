function rcv = mea_frRcv(t, frMat, varargin)
% MEA_FRRCV Calculates firing rate recovery metrics.
%
%   RCV = MEA_FRRCV(T, FRMAT) calculates recovery metrics based on firing rate
%   traces, evaluating the depth of perturbation and the fidelity of recovery.
%   It compares baseline, trough (perturbation), and steady-state (recovery)
%   activity levels.
%
%   INPUTS:
%       t           - (double) Time vector in seconds. The perturbation is
%                     assumed to occur at t=0.
%       frMat       - (double) Firing rate matrix [Units x Time].
%
%   OPTIONAL (Key-Value Pairs):
%       idxTrough   - (double) Index of the trough center around which the
%                     trough window is defined.
%       uGood       - (logical) Vector indicating valid "good" units. {Default: all true}
%       binSize     - (double) Bin size of the firing rate in seconds. {Default: 60}
%       flgSave     - (logical) Whether to save the output struct to disk. {Default: false}
%       basepath    - (char) Base path for saving output files. {Default: pwd}
%       flgPlot     - (logical) Whether to plot summary figures. {Default: false}
%
%   OUTPUTS:
%       rcv         - (struct) Recovery metrics structure containing:
%                       .frBsl      - Baseline mean firing rate [Hz].
%                       .frTrough   - Trough mean firing rate [Hz].
%                       .frSs       - Steady-state (Recovery) mean firing rate [Hz].
%                       .pertDepth  - Perturbation Depth: Log2 ratio of Baseline / Trough.
%                                     (Positive value indicates rate reduction).
%                       .rcvErr     - Recovery Error: Absolute log2 ratio of Steady-State / Baseline.
%                                     (0 indicates perfect recovery).
%                       .rcvBsl     - Recovery Percentage: Ratio of Steady-State / Baseline in percent.
%                                     (100% indicates full recovery).
%                       .rcvDiff    - Recovery Difference: Absolute difference between Steady-State and Baseline rates [Hz].
%                       .rcvGain    - Recovery Gain: Log2 ratio of Steady-State / Trough.
%                                     (Indicates amount of recovery from trough).
%                       .rcvWork    - Recovery Work: Ratio of Recovery Gain to Perturbation Depth.
%                       .bslTime    - Baseline Recovery Time [s]: Time from trough to reach 67% of the distance to Baseline.
%                       .rcvTime    - Steady-State Recovery Time [s]: Time from trough to reach 67% of the distance to Steady-State.
%                       .rcvSlope   - Recovery Slope [Hz/min]: Maximum slope during recovery.
%                       .normSlope  - Normalized Slope [%/min]: Slope normalized by recovery range.
%                       .spkDfct    - Spike Deficit: Log2 ratio of Expected vs Observed AUC.
%                       .uRcv       - (logical) Units that successfully recovered.
%                       .uPert      - (logical) Units that were significantly perturbed.
%
%   NOTES:
%       - Bias Correction: To prevent undefined log ratios when firing rates are
%         zero, a small constant 'c' (1 spike/hour) is added to ALL units if
%         ANY 'good' unit (uGood) has a rate lower than 'c'. This check is
%         performed independently for Baseline, Trough, and Steady-state periods.
%         This ensures unbiased shifting of rates.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 't', @isnumeric);
addRequired(p, 'frMat', @isnumeric);
addParameter(p, 'idxTrough', [], @isnumeric);
addParameter(p, 'uGood', [], @(x) islogical(x) || isnumeric(x));
addParameter(p, 'binSize', 60, @isnumeric);
addParameter(p, 'flgSave', false, @islogical);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgPlot', false, @islogical);

parse(p, t, frMat, varargin{:});

idxTrough = p.Results.idxTrough;
uGood = p.Results.uGood;
binSize = p.Results.binSize;
flgSave = p.Results.flgSave;
basepath = p.Results.basepath;
flgPlot = p.Results.flgPlot;

% Infer defaults
nUnits = size(frMat, 1);
if isempty(uGood)
    uGood = true(nUnits, 1);
end
uGood = logical(uGood);

% Ensure time is a row vector
if iscolumn(t)
    t = t';
end


%% ========================================================================
%  INITIALIZATION
%  ========================================================================

% Constants
c = 1 / 3600;           % 1 spike per hour (~0.0002 Hz, min FR recorded ~0.004 Hz)
thrPert = 1.00;         % Threshold for significant perturbation (50% reduction)
thrRcv = 0.5;           % Threshold for recovery (linear fraction of range)

% Find perturbation index (t=0)
idxPert = find(t >= -1e-9, 1, 'first');

% Initialize output vectors
bslTime   = nan(nUnits, 1);
rcvTime   = nan(nUnits, 1);
bslIdx    = nan(nUnits, 1);
rcvIdx    = nan(nUnits, 1);
rcvSlope  = nan(nUnits, 1);
normSlope = nan(nUnits, 1);
spkDfct   = nan(nUnits, 1);


%% ========================================================================
%  CALCULATE RATES
%  ========================================================================

% Define Windows
winMarg = round(5 * 60 / binSize);
winDur  = round(60 * 60 / binSize);

% Baseline Rate
% -------------------------------------------------------------------------
% Window: 60 bins (3hr abs) ending 5 bins (30min) before perturbation
bslEnd = idxPert - winMarg;
bslStart = max(1, bslEnd - winDur);
bslWin = bslStart : bslEnd;
frBsl = mean(frMat(:, bslWin), 2, 'omitnan');

% Conditional Bias Correction (Baseline)
if any(frBsl(uGood) < c)
    warning('[MEA_FRRCV] %d units have BSL FR < c', sum(frBsl(uGood) < c))
    frBsl = frBsl + c;
end

% Steady State Rate
% -------------------------------------------------------------------------
% Window: 60 bins (6hr abs) ending 5 bins (30min) before end of recording
ssEnd = length(t) - winMarg;
ssStart = ssEnd - winDur;
ssWin = ssStart : ssEnd;
frSs = mean(frMat(:, ssWin), 2, 'omitnan');

% Conditional Bias Correction (Steady State)
if any(frSs(uGood) < c)
    frSs = frSs + c;
end

% Trough Rate
% -------------------------------------------------------------------------
% Window: 20 bins (2hr abs) starting 5 bins (30 min abs) before idxTrough
troughStart = max(idxPert, idxTrough - winMarg);
troughEnd = idxTrough + 15;
troughWin = troughStart : troughEnd;
frTrough = mean(frMat(:, troughWin), 2, 'omitnan');

% Conditional Bias Correction (Trough)
% Censored Regression
% -------------------------------------------------------------------------
% Calculate Theoretical Minimum (Limit of Detection)
% We clamp all smoothed values below this limit to 'frMin' and flag them
% as censored. This ensures the Tobit model treats them as "At or Below Limit".

troughDur = length(troughWin) * binSize;
frMin = 1 / troughDur;

% Create Censoring Mask & Clamp
censMask = frTrough < (frMin + 1e-9); % Epsilon for float tolerance
frTrough(censMask) = frMin;

tbl = table();
tbl.frBsl = frBsl(uGood);
tbl.frTrough = frTrough(uGood);
tbl.censMask = censMask(uGood); 
frml = 'frTrough ~ frBsl';

% Run Tobit Imputation
[frTroughRcv, ~] = mea_lmCens(tbl, frml, 'censVar', 'censMask');
frTrough(uGood) = frTroughRcv;


%% ========================================================================
%  CALCULATE METRICS
%  ========================================================================

% Perturbation Depth
pertDepth = log2(frBsl ./ frTrough);

% Recovery Fidelity
rcvErr  = abs(log2(frSs ./ frBsl));
rcvBsl  = (frSs ./ frBsl) * 100;
rcvDiff = abs(frSs - frBsl);

% Recovery Gain and Work
rcvGain = log2(frSs ./ frTrough);
rcvGain(rcvGain < 0) = 0; % Clamp negative gain to 0
rcvWork = rcvGain ./ pertDepth;

% Perturbed Units
uPert = pertDepth >= thrPert;

% Recoverd Units
% Significant perturbation AND reached thr% of recovery range
uRcv = uPert & (frSs >= frTrough + thrRcv * (frBsl - frTrough));

% Assert uGood on unit indices
uRcv(~uGood) = false;
uPert(~uGood) = false;


%% ========================================================================
%  KINETICS ANALYSIS
%  ========================================================================

tPost = t(idxPert:end); % Post-perturbation time vector

for iUnit = 1:nUnits
    if ~uGood(iUnit)
        continue;
    end

    % Spike Deficit (AUC Metric)
    % ---------------------------------------------------------------------
    if ~isnan(frBsl(iUnit))
        aucExp = trapz(tPost, ones(size(tPost)) * frBsl(iUnit));
        aucObs = trapz(tPost, frMat(iUnit, idxPert:end));
        spkDfct(iUnit) = log2(aucExp / (aucObs + 1));
    end

    % Recovery Time & Indices
    % ---------------------------------------------------------------------
    rcvCurve = frMat(iUnit, idxTrough:end);

    % Range for recovery calculation
    rangeVal = frBsl(iUnit) - frTrough(iUnit);

    % Time to Recover to Baseline (Approximate)
    % Target: Trough + 67% of (Baseline - Trough)
    bslTrgt = frTrough(iUnit) + thrRcv * rangeVal;
    bslPnt = find(rcvCurve >= bslTrgt, 1, 'first');

    if ~isempty(bslPnt)
        absIdx = idxTrough + bslPnt - 1;
        bslIdx(iUnit) = absIdx;
        bslTime(iUnit) = t(absIdx) - t(idxTrough);
    else
        bslIdx(iUnit) = length(t);
        bslTime(iUnit) = t(end) - t(idxTrough);
    end

    % Time to Recover to Steady State (Actual Recovery)
    % Target: Trough + 67% of (SteadyState - Trough)
    rcvRange = frSs(iUnit) - frTrough(iUnit);
    rcvTrgt = frTrough(iUnit) + thrRcv * rcvRange;
    rcvPnt = find(rcvCurve >= rcvTrgt, 1, 'first');

    if ~isempty(rcvPnt)
        absIdx = idxTrough + rcvPnt - 1;
        rcvIdx(iUnit) = absIdx;
        rcvTime(iUnit) = t(absIdx) - t(idxTrough);
    else
        rcvIdx(iUnit) = length(t);
        rcvTime(iUnit) = t(end) - t(idxTrough);
    end

    % Slope Analysis
    % ---------------------------------------------------------------------
    if ~isempty(rcvCurve) && length(rcvCurve) > 1
        % Max instantaneous slope (Hz/min)
        currSlope = (max(diff(rcvCurve)) / binSize) * 60;
        rcvSlope(iUnit) = currSlope;

        % Normalized slope (%/min)
        if rcvRange ~= 0
            normSlope(iUnit) = (currSlope / rcvRange) * 100;
        end
    end
end


%% ========================================================================
%  OUTPUT CONSTRUCTION
%  ========================================================================

rcv.frBsl       = frBsl;
rcv.frTrough    = frTrough;
rcv.frSs        = frSs;
rcv.pertDepth   = pertDepth;
rcv.rcvErr      = rcvErr;
rcv.rcvBsl      = rcvBsl;
rcv.rcvGain     = rcvGain;
rcv.rcvWork     = rcvWork;
rcv.rcvDiff     = rcvDiff;
rcv.bslTime     = bslTime;
rcv.rcvTime     = rcvTime;
rcv.rcvSlope    = rcvSlope;
rcv.normSlope   = normSlope;
rcv.spkDfct     = spkDfct;
rcv.uRcv        = uRcv;
rcv.uPert       = uPert;

rcv.info.idxTrough  = idxTrough;
rcv.info.winTrough  = [troughStart, troughEnd] * 60;
rcv.info.winBsl     = [bslStart, bslEnd] * 60;
rcv.info.winSs      = [ssStart, ssEnd] * 60;
rcv.info.binSize    = binSize;

if flgSave
    [~, basename] = fileparts(basepath);
    save(fullfile(basepath, [basename, '.frRcv.mat']), 'rcv', '-v7.3');
end


%% ========================================================================
%  PLOTTING
%  ========================================================================

if flgPlot
    % Select n random good units
    goodIdx = find(uGood);
    nPlot = min(6, length(goodIdx));

    % Use current random stream to shuffle
    rng('shuffle');
    shuffledIdx = goodIdx(randperm(length(goodIdx)));
    plotIdx = shuffledIdx(1:nPlot);

    % Overwrite random for comparison of raw with model
    plotIdx = [1 : 6];

    figure('Position', [100, 100, 1200, 800], 'Name', 'MEA Recovery Metrics');
    nRows = ceil(nPlot/2);

    % Convert time to Hours for plotting
    tHr = t / 3600;

    for iPlot = 1:nPlot
        uIdx = plotIdx(iPlot);
        subplot(nRows, 2, iPlot);
        hold on;

        % TRACES
        % -----------------------------------------------------------------
        % Raw FR
        plot(tHr, frMat(uIdx, :), 'Color', [0.7 0.7 0.7], 'LineWidth', 1, ...
            'DisplayName', 'Raw FR');

        % Baseline Level
        plot(tHr(bslWin), repmat(frBsl(uIdx), size(bslWin)), ...
            'Color', [0 0.4470 0.7410], 'LineWidth', 2, 'DisplayName', 'Baseline');

        % Trough Level
        plot(tHr(troughWin), repmat(frTrough(uIdx), size(troughWin)), ...
            'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2, 'DisplayName', 'Trough');

        % Steady State Level
        plot(tHr(ssWin), repmat(frSs(uIdx), size(ssWin)), ...
            'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2, 'DisplayName', 'Steady State');


        % METRIC INDICATORS
        % -----------------------------------------------------------------
        % Perturbation Depth Line
        plot([tHr(idxTrough), tHr(idxTrough)], [frTrough(uIdx), frBsl(uIdx)], ...
            '--', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2, ...
            'DisplayName', 'Pert. Depth');

        % Recovery Gain Line
        idxSSCenter = round(mean(ssWin));
        plot([tHr(idxSSCenter), tHr(idxSSCenter)], [frTrough(uIdx), frSs(uIdx)], ...
            '--', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2, ...
            'DisplayName', 'Rcv. Gain');

        % Recovery Time Point
        if ~isnan(rcvIdx(uIdx))
            tRcvPlot = tHr(rcvIdx(uIdx));
            tRcvValHr = rcvTime(uIdx) / 3600;

            plot(tRcvPlot, frMat(uIdx, rcvIdx(uIdx)), 'p', 'MarkerSize', 10, ...
                'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', ...
                'DisplayName', sprintf('RcvTime: %.2f hr', tRcvValHr));
        end

        % Spike Deficit Area
        idxPost = idxPert:length(t);
        if ~isempty(idxPost)
            tPostHr = tHr(idxPost);
            frObs = frMat(uIdx, idxPost);
            frExp = ones(size(frObs)) * frBsl(uIdx);

            xFill = [tPostHr, fliplr(tPostHr)];
            yFill = [frExp, fliplr(frObs)];
            fill(xFill, yFill, [0.5 0.5 0.5], 'FaceAlpha', 0.1, 'EdgeColor', 'none', ...
                'DisplayName', sprintf('SpkDfct: %.2f', spkDfct(uIdx)));
        end

        % Text Labels
        text(tHr(idxTrough), mean([frTrough(uIdx), frBsl(uIdx)]), ...
            sprintf(' %.2f', pertDepth(uIdx)), ...
            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', ...
            'FontSize', 8, 'Color', [0.8500 0.3250 0.0980], 'FontWeight', 'bold');

        text(tHr(idxSSCenter), mean([frTrough(uIdx), frSs(uIdx)]), ...
            sprintf(' %.2f', rcvGain(uIdx)), ...
            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', ...
            'FontSize', 8, 'Color', [0.4660 0.6740 0.1880], 'FontWeight', 'bold');

        % Formatting
        xlabel('Time (Hr)');
        ylabel('Rate (Hz)');
        title(sprintf('Unit %d', uIdx));
        grid on;
        axis tight;
        xline(0, '-', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');

        if iPlot == 1
            legend('Location', 'best', 'NumColumns', 2, 'FontSize', 7);
        end
    end
end

end     % EOF



%% ========================================================================
%  NOTE: THEORETICAL VS. OBSERVED MINIMUM & CENSORING STRATEGY
%  ========================================================================
%  OBSERVATION
%   The calculated Theoretical Minimum firing rate (1 spike / window duration)
%   is often orders of magnitude higher than the absolute minimum non-zero
%   value observed in the data trace.
%
%  EXPLANATION: SMOOTHING ARTIFACTS
%   This discrepancy confirms that the input firing rate matrix (frMat) has
%   undergone temporal smoothing. While biological spiking is discrete,
%   smoothing filters redistribute this energy into continuous curves. The
%   extremely low values observed in the data (e.g., 1e-6 Hz) represent the
%   asymptotic tails of these smoothing kernels ("fractional spikes") rather
%   than physical firing events.
%
%  STRATEGY: WHERE TO APPLY THRESHOLDS
%   Crucially, one must NOT apply a hard threshold (cap) during the pre-
%   processing or denoising stages (e.g., mea_frPrep or fr_denoise). Doing
%   so would be destructive. Because smoothing spreads a single spike across
%   multiple bins, the peak value in any single bin is often lower than
%   1/binSize. A pre-processing cap would erroneously erase these valid,
%   isolated spiking events.
%
%  IMPLICATION FOR ANALYSIS
%   The correct approach is to handle these values only during the aggregation
%   stage (mea_frRcv). When averaging over a specific window (e.g., Trough),
%   any mean rate falling below the resolution limit of that specific window
%   (1 / WindowDuration) should be treated as "Censored" (Left-Censored at
%   Limit L) for the purpose of Tobit regression, rather than deleted or
%   clamped continuously in the raw trace.
%  ========================================================================


