function isiVal = brst_isiValley(spktimes, varargin)
% BRST_ISIVALLEY Detects the ISI valley separating bursts from non-bursts.
%
%   isiValley = BRST_ISIVALLEY(SPKTIMES, ...) calculates the inter-spike
%   intervals (ISIs) for all cells, aggregates them, and finds the valley
%   in the log-log histogram (or KDE) that typically separates burst ISIs
%   from tonic ISIs.
%
%   INPUTS:
%       spktimes    - (cell) Spike times per unit (e.g., {unit1, unit2}).
%                     Times should be in seconds.
%       varargin    - (param/value) Optional parameters:
%                     'nSpks'        : (num) N spikes for ISI calculation {3}
%                     'rngVal'       : (vec) [min, max] ISI (s) to search for valley {[0.005, 0.25]}
%                     'bandwidth'    : (num) Bandwidth for KDE (log10 scale) {0.05}
%                     'flgPlot'      : (log) Plot the PDF and detected valley {true}
%
%   OUTPUTS:
%       isiVal      - (num) The ISI value (in seconds) at the valley.
%
%   EXAMPLE:
%       val = brst_isiValley(spktimes, 'flgPlot', true);
%

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addParameter(p, 'nSpks', 3, @isnumeric);
addParameter(p, 'rngVal', [0.01, 0.25], @isnumeric);
addParameter(p, 'bandwidth', 0.05, @isnumeric);
addParameter(p, 'flgPlot', true, @islogical);

parse(p, varargin{:});
nSpks        = p.Results.nSpks;
rngVal       = p.Results.rngVal;
bw           = p.Results.bandwidth;
flgPlot      = p.Results.flgPlot;


%% ========================================================================
%  PREPARATIONS
%  ========================================================================

% Flatten spktimes if needed and handle single vector input
if isnumeric(spktimes)
    spktimes = {spktimes};
end

% Collect all ISIs
% Calculate ISI_N: Distance between spike i and spike i + (nSpks - 1)
isiVec = [];
nISI   = nSpks - 1;

for iUnit = 1:length(spktimes)
    st = spktimes{iUnit}(:);
    if isempty(st) || length(st) < nSpks
        continue;
    end

    % Simple diff for nSpks=2, higher order diff for nSpks > 2
    isiVec = [isiVec; st(1+nISI:end) - st(1:end-nISI)]; %#ok<AGROW>
end

% Remove invalid ISIs (refractory period violations)
isiVec(isiVec <= 0.0005) = [];

% Remove extremely long ISIs that compromise range, essentially limites
% minimum firing rate to 0.0001 Hz
isiVec(isiVec >= 10000) = [];

% Validation
if isempty(isiVec)
    warning('brst_isiValley:NoISIs', 'No valid ISIs found.');
    isiVal = nan;
    return;
end


%% ========================================================================
%  CALCULATION
%  ========================================================================

% Log transform
isiLog = log10(isiVec);

% Kernel Density Estimation
% Define evaluation grid: from min to max, with sufficient resolution
gridMin = min(isiLog);
gridMax = max(isiLog);
nPoints = 1000;
xIdx = linspace(gridMin, gridMax, nPoints);

% Use ksdensity for robust smoothing
[isiPdf, xIdx] = ksdensity(isiLog, xIdx, 'Bandwidth', bw);

% Limit search to the specified range
rngLog = log10(rngVal);
rngIdx = xIdx >= rngLog(1) & xIdx <= rngLog(2);

% Segment of PDF/Grid to search
pdfSeg = isiPdf(rngIdx);
xiSeg  = xIdx(rngIdx);

% Find local minima (valleys) in this segment
[~, locsMin] = findpeaks(-pdfSeg);

if ~isempty(locsMin)
    % If multiple valleys, pick the deepest one (lowest density)
    [~, bestMinIdx] = min(pdfSeg(locsMin));
    idxInSeg = locsMin(bestMinIdx);
    isiVal = 10^xiSeg(idxInSeg);
else
    % Fallback: Absolute minimum in the range
    [~, idxMin] = min(pdfSeg);
    isiVal = 10^xiSeg(idxMin);
end


%% ========================================================================
%  PLOTTING
%  ========================================================================

if flgPlot

    [hFig, hAx] = plot_axSize('szOnly', false, 'flgFullscreen', true, ...
        'flgPos', true);

    % Plot Histogram
    hHist = histogram(isiLog, 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    hold on;

    % Plot KDE
    hPlt = plot(xIdx, isiPdf, 'k-', 'LineWidth', 2);

    % Plot Valley
    if ~isnan(isiVal)
        xline(log10(isiVal), 'r--', 'LineWidth', 1.5, 'Label', sprintf('Valley: %.1f ms', isiVal*1000));
        plot(log10(isiVal), interp1(xIdx, isiPdf, log10(isiVal)), 'ro', 'MarkerFaceColor', 'r');
    end

    % Convert X-axis to Frequency (Hz)
    % log10(ISI) = -log10(Freq)
    xLims = xlim(hAx);
    minLogFreq = -xLims(2);
    maxLogFreq = -xLims(1);
    tickPowers = ceil(minLogFreq):floor(maxLogFreq);
    xticks(hAx, tickPowers);
    xticklabels(hAx, arrayfun(@(p) sprintf('%.0d', 10^-p), tickPowers, 'UniformOutput', false));

    % Formatting
    xlabel(hAx, 'Frequency (Hz)');
    ylabel('PDF')
    set(hAx, 'TickLabelInterpreter', 'tex');
end

end     % EOF


%% ========================================================================
%  NOTE: ISI VALLEY METHOD
%  ========================================================================
% The ISI valley method identifies a data-driven temporal threshold by
% finding the minimum (trough) between the two peaks of a bimodal
% log-ISI distribution. This threshold distinguishes between the two
% primary firing states of a neuron: high-frequency clusters (bursts)
% and low-frequency background activity (tonic firing).
%
% This approach explicitly separates the concept of burst "speed" (ISI
% threshold) from burst "validity" (Spike count).
% The threshold defines the maximum allowable speed for an interval to
% be considered part of a burst. However, a single fast interval only
% suggests a potential burst; the event is only validated if it
% sustains this high-frequency state across multiple consecutive
% intervals as defined by the minimum spike criterion.
%  ========================================================================

%% ========================================================================
%  NOTE: LOG-ISI DISTRIBUTION ANALYSIS
%  ========================================================================
% Visualizing the distribution of inter-spike intervals in logarithmic
% space is the gold standard for verifying burst detection parameters.
% Unlike linear histograms, which compress all high-frequency data
% into the first few bins, the log-space histogram clearly displays
% the multi-scale nature of neuronal timing.
%
% A well-chosen burst threshold should fall exactly in the "valley" or
% lowest point between the short-interval peak (bursts) and the long-
% interval peak (tonic firing). If a chosen threshold falls directly on
% a peak, it will be highly unstable, as minor fluctuations in firing
% rate will cause massive changes in the detected burst metrics.
%  ========================================================================

%% ========================================================================
%  NOTE: HIGH-FREQUENCY DISCRETIZATION
%  ========================================================================
% The "comb-like" spikes and fine valleys visible in the high-frequency
% (low ISI) region of the log-ISI histogram are primarily non-biological
% artifacts of temporal discretization and digital sampling.
%
% * Sampling Rate Limits: At high frequencies, the time between discrete
%   samples (e.g., 50 s at 20 kHz) is a large percentage of the total
%   ISI. Because spikes are assigned to a fixed temporal grid, intervals
%   can only exist as integer multiples of the sampling period.
%
% * Log-Transform Artifacts: Mapping these fixed integer steps into the
%   continuous log-domain stretches the mathematical gaps between them.
%   This creates artificial "peaks" where data exists and "empty bins"
%   where no interval is mathematically possible, resulting in the
%   discretized comb appearance.
%
% * Biological Resonance: While the fine comb is a sampling artifact,
%   larger distinct peaks in this region (e.g., 100–200 Hz) represent the
%   natural "stereotypical" firing frequency of hippocampal neurons
%   during a sustained bursting state.
%
% To avoid mistaking these sampling artifacts for the true burst-tonic
% valley, a Kernel Density Estimation (KDE) with a slightly higher
% bandwidth is typically required to smooth over the "comb" and reveal
% the underlying biological distribution.
%  ========================================================================

%% ========================================================================
%  NOTE: THRESHOLD STRATEGY (FIXED VS. ADAPTIVE)
%  ========================================================================
% Selecting between a global (fixed) and individual (adaptive) threshold
% is a critical decision that impacts the statistical validity of
% cross-group comparisons and pharmacological tracking. While adaptive
% methods maximize detection sensitivity for heterogeneous units, they
% introduce specific technical risks that can confound biological
% interpretations.
%
% THE NORMALIZATION BIAS (MASKING EFFECTS). A significant risk of
% unit-specific adaptive thresholds is the unintended "normalization" of
% the very biological effects under investigation. If a genotype, such as
% MCU-KO, or a drug treatment causes a global slowing of network activity,
% an adaptive algorithm will shift its threshold toward longer intervals to
% find the new statistical valley. Consequently, the metric for
% "burstiness" (e.g., fraction of spikes in bursts) may appear identical
% across groups because the definition of a burst has been expanded for the
% slower preparation, effectively masking a primary biological phenotype.
%
% PHARMACOLOGICAL TRACKING AND BIMODALITY. Adaptive thresholds rely on the
% presence of a bimodal log-ISI distribution to function accurately. During
% pharmacological perturbations like Baclofen application, firing rates may
% drop significantly, causing the distinct "burst" and "tonic" peaks to
% merge into a single unimodal distribution. In these states,
% valley-detection algorithms often return unstable or arbitrary values,
% making longitudinal tracking and recovery analysis technically
% unreliable.
%
% CORRELATION AND MATHEMATICAL INDEPENDENCE. To study the correlation
% between baseline burstiness and drug response, the baseline metric must
% serve as a stable predictor. Using an individual threshold derived from
% the baseline recording creates a mathematical dependency between the
% predictor and the outcome. A fixed threshold provides an absolute,
% independent frame of reference that ensures the definition of a "burst"
% remains constant throughout the baseline, perturbation, and recovery
% phases.
%
% FIXED THRESHOLDS AS PHYSIOLOGICAL LIMITS. Fixed thresholds are often
% interpreted as "physiological speed limits" defined by the stereotypical
% firing frequencies of specific neuronal types. For hippocampal CA1
% preparations, a threshold of ~13ms identifies the high-frequency core of
% a burst regardless of the unit's background firing rate.
%
% STRATEGIC COMPROMISE: GLOBAL ADAPTATION
% A robust intermediate strategy involves calculating a single "Global
% Adaptive Threshold" by aggregating all intervals from baseline
% Wild-Type recordings. This approach uses the data-driven valley of
% the specific biological system but applies it as a fixed reference
% across all genotypes and conditions, preserving the ability to
% detect absolute shifts in network state.
%  ========================================================================