function dyn = brst_dynamics(brst, spktimes, varargin)
% BRST_DYNAMICS Computes time-varying burst dynamics (Population Ready).
%
%   dyn = BRST_DYNAMICS(BRST, SPKTIMES, ...) converts discrete burst events
%   (from brst_detect) into continuous time-series using density estimation
%   (for rates) and interpolation with adaptive masking (for properties).
%
%   INPUTS:
%       brst        - (struct) Output from brst_detect.m
%       spktimes    - (cell) Spike times per unit.
%       varargin    - (param/value) Optional parameters:
%                     'binSize'   : (num) Time bin size {60} (s)
%                     'ksd'       : (num) Gaussian kernel sigma {300} (s)
%                     'winSm'     : (num) Smoothing window {10} (events)
%                     'ibiPct'    : (num) Percentile for masking gap {99}
%                     'flgPlot'   : (log) Plot dynamics {true}
%                     'flgSave'   : (log) Save result as brstDyn {false}
%                     'basepath'  : (char) Base path for saving {pwd}
%
%   OUTPUTS:
%       dyn         - (struct) Dynamics structure.
%                     .time       : (1 x nTime) Time vector.
%                     .eventRate  : (nUnits x nTime) Burst event rate (Hz).
%                     .frTot      : (nUnits x nTime) Total firing rate (Hz).
%                     .frBspk     : (nUnits x nTime) Burst spike firing rate (Hz).
%                     .frSspk     : (nUnits x nTime) Single spike firing rate (Hz).
%                     .pBspk      : (nUnits x nTime) Burst fraction (0-1).
%                     .dur        : (nUnits x nTime) Duration (s).
%                     .nBspk      : (nUnits x nTime) Spikes per burst.
%                     .freq       : (nUnits x nTime) Intra-burst freq (Hz).
%                     .ibi        : (nUnits x nTime) Inter-burst interval (s).
%
%   NOTE:
%   Properties (.dur, .nBspk, etc) are interpolated to the global time
%   grid. Periods of silence exceeding the 'ibiPct' percentile of a unit's
%   IBI distribution are masked as NaN to distinguish "undefined" states
%   from zero values.
%
%   See also: BRST_DETECT, BRST_STATS

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'brst', @isstruct);
addRequired(p, 'spktimes', @iscell);
addParameter(p, 'binSize', 60, @isnumeric);
addParameter(p, 'ksd', 300, @isnumeric);
addParameter(p, 'winSm', 10, @isnumeric);
addParameter(p, 'ibiPct', 95, @isnumeric);
addParameter(p, 'flgPlot', true, @islogical);
addParameter(p, 'flgSave', false, @islogical);
addParameter(p, 'basepath', pwd, @ischar);

parse(p, brst, spktimes, varargin{:});
binSize   = p.Results.binSize;
ksd       = p.Results.ksd;
winSm     = p.Results.winSm;
ibiPct    = p.Results.ibiPct;
flgPlot   = p.Results.flgPlot;
flgSave   = p.Results.flgSave;
basepath  = p.Results.basepath;


%% ========================================================================
%  INITIALIZE
%  ========================================================================

nUnits = length(brst.times);

% Time vector (excluding tail)
maxTime = max(cellfun(@max, spktimes, 'uni', true));
chunks = n2chunks('n', maxTime, 'chunksize', binSize, 'lastChunk', 'exclude');
t = chunks(:, 1)' - 1;            % Bin starts
tEdges = [t, chunks(end, 2)];     % Bin edges
nBins = length(t);

% Initialize Matrices (nUnits x nTime)
dyn.time        = t;
dyn.eventRate   = zeros(nUnits, nBins);
dyn.frTot       = zeros(nUnits, nBins);
dyn.frBspk      = zeros(nUnits, nBins);
dyn.frSspk      = zeros(nUnits, nBins);
dyn.pBspk       = zeros(nUnits, nBins);
dyn.nBspk       = nan(nUnits, nBins);

dyn.dur      = nan(nUnits, nBins);
dyn.freq     = nan(nUnits, nBins);
dyn.ibi      = nan(nUnits, nBins);

% Kernel
kRng = -3*ksd : binSize : 3*ksd;
kd = normpdf(kRng, 0, ksd);
kd = kd / sum(kd);

% Edge Effect Correction Vector
% Normalized kernel convolution with unity vector reveals boundary loss
kCorr = conv(ones(1, nBins), kd, 'same');

%% ========================================================================
%  COMPUTE DYNAMICS
%  ========================================================================

for iUnit = 1:nUnits

    st = spktimes{iUnit};
    bTimes = brst.times{iUnit};

    if isempty(st) || isempty(bTimes)
        continue;
    end

    % Access properties
    nBspk = brst.nBspk{iUnit};
    dur   = brst.dur{iUnit};
    freq  = brst.freq{iUnit};
    ibi   = brst.ibi{iUnit};
    bst   = brst.spktimes{iUnit};

    % ---------------------------------------------------------------------
    % Density Estimation (Continuous)
    % ---------------------------------------------------------------------

    % Burst Event Rate
    bStart = bTimes(:, 1);
    bCounts = histcounts(bStart, tEdges);
    dyn.eventRate(iUnit, :) = (conv(bCounts, kd, 'same') ./ kCorr) / binSize;

    % Burst Densities (Firing Rates)
    spkCounts = histcounts(st, tEdges);
    frTot = (conv(spkCounts, kd, 'same') ./ kCorr) / binSize;

    brstCounts = histcounts(bst, tEdges);
    frBspk = (conv(brstCounts, kd, 'same') ./ kCorr) / binSize;

    % Single Spike Density
    frSspk = frTot - frBspk;

    % Store
    dyn.frTot(iUnit, :)  = frTot;
    dyn.frBspk(iUnit, :) = frBspk;
    dyn.frSspk(iUnit, :) = frSspk;

    % Probability of spikes in bursts
    pBspk = zeros(size(frTot));
    mask = frTot > 1e-6;
    pBspk(mask) = frBspk(mask) ./ frTot(mask);
    dyn.pBspk(iUnit, :) = pBspk;

    % Validation
    % frSspk + frBspk = frTot
    if max(abs((frSspk + frBspk) - frTot)) > 1e-6
        warning('brst_dynamics:ValidationFailed', ...
            'Unit %d: frSspk + frBspk does not equal frTot', iUnit);
    end

    % frTot * (1 - pBspk) = frSspk
    % (Only checking on mask where frTot > 0)
    err = abs((frTot(mask) .* (1 - pBspk(mask))) - frSspk(mask));
    if max(err) > 1e-6
        warning('brst_dynamics:Validation2Failed', ...
            'Unit %d: frTot * (1 - pBspk) != frSspk', iUnit);
    end


    % ---------------------------------------------------------------------
    % Property Tracking (Interpolated & Masked)
    % ---------------------------------------------------------------------

    % Calculate Mask Threshold (Adaptive) and ensure isn't too small
    maskThr = prctile(ibi, ibiPct);
    maskThr = max(maskThr, 2 * binSize);

    % Apply to all props
    bStart = bTimes(:, 1);

    if length(bStart) < 2
        continue
    end

    % Process Properties
    dyn.dur(iUnit, :)   = proc_prop(dur, bStart, t, winSm, maskThr);
    dyn.nBspk(iUnit, :) = proc_prop(nBspk, bStart, t, winSm, maskThr);
    dyn.freq(iUnit, :)  = proc_prop(freq, bStart, t, winSm, maskThr);
    dyn.ibi(iUnit, :)   = proc_prop(ibi, bStart, t, winSm, maskThr);

end


%% ========================================================================
%  PLOTTING
%  ========================================================================

if flgPlot
    tbl = struct2table(rmfield(dyn, 'time'));
    tblGUI_xy(dyn.time, tbl, 'yVar', 'eventRate');
end


%% ========================================================================
%  SAVE
%  ========================================================================

if flgSave
    [~, basename] = fileparts(basepath);
    dynfile = fullfile(basepath, [basename, '.brstDyn.mat']);
    brstDyn = dyn;
    save(dynfile, 'brstDyn');
end

end     % EOF


%% ========================================================================
%  NESTED FUNCTIONS
%  ========================================================================

function bOut = proc_prop(bVal, bStart, t, winSm, maskThr)
% Processes a property vector: Smooth -> Interpolate -> Mask

if length(bVal) < winSm
    % Not enough points to smooth well, just use mean or raw
    smVal = bVal;
else
    smVal = movmedian(bVal, winSm, 'omitnan');
end

% Interpolate to global time grid
% 'linear' provides smooth transitions between events
% 'extrap' with NaN prevents wild guesses outside range
bOut = interp1(bStart, smVal, t, 'linear', NaN);

% Identify Gaps (Masking)

% Gap before first burst
bOut(t < bStart(1) - maskThr) = NaN;

% Gap after last burst
bOut(t > bStart(end) + maskThr) = NaN;

% Interior Gaps
diffs = diff(bStart);
gaps = find(diffs > maskThr);
for iGap = 1:length(gaps)
    idx1 = gaps(iGap);
    idx2 = gaps(iGap) + 1;
    t_start = bStart(idx1) + maskThr;
    t_end   = bStart(idx2) - maskThr;

    if t_end > t_start
        bOut(t > t_start & t < t_end) = NaN;
    end
end

% Ensure row vector for output matrix
bOut = bOut(:)';

end


%% ========================================================================
%  NOTE: BINNED STATS VS. CONTINUOUS DYNAMICS (KINETICS)
%  ========================================================================
%  When analyzing temporal changes, one must choose between sliding-window
%  statistics (brst_stats) and continuous density estimation (brst_dynamics).
%  While both track changes over time, they serve different analytical goals:
%
%  1. TEMPORAL RESOLUTION VS. STABILITY
%     Binned statistics (e.g., 20-min windows) require a minimum number of
%     events per bin to be statistically valid. If a unit fires only two
%     bursts in 20 minutes, the mean is highly volatile. Density estimation
%     (brst_dynamics) uses Gaussian smoothing and interpolation to provide
%     a higher-resolution "state estimate" that is less sensitive to the
%     exact timing of individual triggers.
%
%  2. THE DURATION BIAS (TIME-WEIGHTING)
%     brst_stats treats every burst as one N. A 10s burst and a 1s burst
%     average to 5.5s. In brst_dynamics, the 10s burst occupies 10x more
%     temporal bins than the 1s burst. Consequently, the "average" value
%     of a dynamic trace is weighted by duration. Use stats for
%     measuring the *physical shape* of events and dynamics for measuring
%     the *network's state*.
%
%  3. HANDLING DROPOUTS (SILENCE)
%     In a 20-minute bin with zero bursts, a structural metric (like dur)
%     is simply NaN. In dynamics, we use adaptive masking (ibiPct) to
%     distinguish between "the neuron is still in a bursting state but
%     hasn't fired yet" and "the neuron has dropped out of the bursting
%     regime entirely".
%
%  SUMMARY:
%  - Use brst_stats for: Formal genotype comparisons (WT vs KO) and
%    calculating physical burst geometry across experimental phases.
%  - Use brst_dynamics for: Visualizing the kinetic "flow" of recovery,
%    drug onset, or high-resolution population transitions.
%  ========================================================================

%% ========================================================================
%  NOTE: THE NECESSITY OF INTERPOLATION FOR POPULATION ANALYSIS
%  ========================================================================
%  This function uses interpolation to transform asynchronous unit data
%  into a synchronized matrix format suitable for genotype comparisons.
%  Because Unit A might burst at t=10s and Unit B at t=15s, simply
%  averaging their properties without a shared temporal grid is impossible.
%  Interpolation projects the property value of a specific burst onto the
%  global time vector T. This creates a continuous "estimate" of what the
%  burst duration would be if a burst occurred at that exact moment, based
%  on the nearest real observations.
%  ========================================================================

%% ========================================================================
%  NOTE: DISTINGUISHING ZERO FROM UNDEFINED STATES
%  ========================================================================
%  If we do not mask our interpolated data, a unit that bursts at the start
%  of a recording and then remains silent for the remainder will appear to
%  have a constant, "flat" duration across the entire experiment. This is
%  biologically misleading because the property (duration) only exists
%  during the event itself. While a burst rate of zero is a valid
%  measurement of inactivity, a burst duration of zero is a physical
%  impossibility in this detection framework (as the algorithm requires
%  minimum spikes and duration). Therefore, periods of silence are treated
%  as NaN (undefined) rather than zero.
%  ========================================================================

%% ========================================================================
%  NOTE: DEFINING THE MASKING WINDOW NON-ARBITRARILY
%  ========================================================================
%  To avoid arbitrary thresholds, we determine the masking window based on
%  the statistical behavior of the Inter-Burst Interval (IBI). The function
%  calculates the 95th (or user-specified) percentile of the IBI for each
%  unit. If the time elapsed since the last detected burst exceeds this
%  calculated "Expected IBI," the property is masked as NaN. This ensures
%  that the population average is only calculated from units that are
%  actively participating in the network's bursting state at that specific
%  time.
%  ========================================================================

%% ========================================================================
%  NOTE: MOVING MEDIANS AND MISSING DATA
%  ========================================================================
%  Applying a movmedian across the entire length of the recording is only
%  possible if you have a value for every time bin, which effectively
%  forces you to assume a value (like zero) for bins without bursts. This
%  is incorrect for structural metrics like brstDur or nBspk. By performing
%  the movmedian on the event index first - meaning we smooth across the
%  last N bursts regardless of when they happened - we preserve the
%  internal structural integrity of the burst properties. The interpolation
%  and subsequent masking then allow us to align these "active state"
%  values across the population while correctly identifying periods where
%  the unit has dropped out of the bursting state.
%  ========================================================================
%
%  %% ========================================================================
%  NOTE: EDGE EFFECT CORRECTION
%  ========================================================================
%  Convolution with a Gaussian kernel assumes zero-padding outside the
%  signal boundaries. This causes rate estimates to artificially drop at
%  the start and end of the recording because the kernel "sees" zeros where
%  there is simply no data. To correct this, we calculate a normalization
%  vector (kCorr) by convolving the kernel with a vector of ones. Dividing
%  the raw convolution by this vector compensates for the missing mass at
%  the edges, ensuring the rate estimate remains accurate throughout the
%  entire timeline without altering the smoothed shape in the center.
%  ========================================================================

%% ========================================================================
%  NOTE: STRUCTURAL (EVENT-BASED) VS. KINETIC (DENSITY) ANALYSIS
%  ========================================================================
%  Analyzing burst data requires choosing between two distinct methods:
%  structural analysis, which focuses on individual events, and kinetic
%  analysis, which focuses on continuous time-varying density. Choosing
%  between them determines whether the resulting data reflects the physical
%  shape of a burst or the prevailing state of the neuronal network over
%  time.
%
%  STRUCTURAL ANALYSIS (EVENT-BASED)
%  Structural analysis treats every detected burst as a single,
%  equivalent data point regardless of when it occurs. This
%  method is used to calculate the physical geometry of bursts, such as
%  average duration, spikes per burst, and internal frequency.
%  Because each event carries the same weight, this approach provides
%  the most accurate measurement of how a burst is built.
%
%  * Measures physical burst shape.
%  * Weights every burst equally.
%  * Best for baseline characterization.
%  * Ignores time between events.
%
%  KINETIC ANALYSIS (DENSITY-BASED)
%  Kinetic analysis converts discrete bursts into a continuous time-series
%  using Gaussian kernels and interpolation. This transforms
%  burst events into a "state estimate" that exists for every time bin
%  of the recording. This approach is necessary for comparing
%  bursting behavior directly to firing rate recovery or other
%  time-varying signals on a synchronized global grid.
%
%  * Measures time-varying network state.
%  * Weights events by duration.
%  * Synchronizes units across time.
%  * Best for recovery modeling.
%
%  THE DURATION BIAS IN KINETIC METRICS
%  A significant difference between these methods is how they handle
%  burst length. Kinetic analysis samples property values from a
%  continuous trace. A long burst occupies more time bins than
%  a short burst, meaning it contributes more to the final average
%  value. While this is helpful for measuring the "total
%  impact" of bursting on a network, it is mathematically biased if the
%  goal is to measure the average physical size of individual events.
%
%  * Long bursts dominate density.
%  * Short bursts are diluted.
%  * Creates time-weighted state estimates.
%
%  STATE MASKING AND SILENCE
%  Kinetic analysis uses adaptive masking based on the Inter-Burst Interval
%  (IBI) to handle periods where a unit stops bursting. If a unit is silent
%  for longer than its expected interval, the density trace becomes NaN.
%  This allows statistical models to identify that the bursting "state" has
%  ended, whereas structural analysis simply runs out of events to average.
%
%  * Identifies burst-state dropout.
%  * Prevents zero-duration artifacts.
%  * Uses adaptive IBI percentiles.
%  ========================================================================
