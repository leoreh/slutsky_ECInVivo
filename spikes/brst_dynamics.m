function dyn = brst_dynamics(brst, spktimes, varargin)
% BRST_DYNAMICS Computes time-varying burst dynamics (Population Ready).
%
%   dyn = BRST_DYNAMICS(BRST, SPKTIMES, ...) converts discrete burst events
%   (from brst_maxInt) into continuous time-series using density estimation
%   (for rates) and interpolation with adaptive masking (for properties).
%
%   INPUTS:
%       brst        - (struct) Output from brst_maxInt.m
%       spktimes    - (cell) Spike times per unit.
%       varargin    - (param/value) Optional parameters:
%                     'binSize'   : (num) Time bin size {60} (s)
%                     'kernelSD'  : (num) Gaussian kernel sigma {300} (s)
%                     'smoothWin' : (num) Smoothing window {10} (events)
%                     'ibiPct'    : (num) Percentile for masking gap {99}
%                     'flgPlot'   : (log) Plot dynamics {true}
%                     'flgSave'   : (log) Save result as brstDyn {false}
%                     'basepath'  : (char) Base path for saving {pwd}
%
%   OUTPUTS:
%       dyn         - (struct) Dynamics structure.
%                     .time     : (1 x nTime) Time vector.
%                     .rate     : (nUnits x nTime) Burst rate (Hz).
%                     .bfrac    : (nUnits x nTime) Burst fraction.
%                     .dur      : (nUnits x nTime) Duration (s).
%                     .nspks    : (nUnits x nTime) Spikes per burst.
%                     .freq     : (nUnits x nTime) Intra-burst freq (Hz).
%                     .ibi      : (nUnits x nTime) Inter-burst interval (s).
%
%   NOTE:
%   Properties (.dur, .nspks, etc) are interpolated to the global time
%   grid. Periods of silence exceeding the 'ibiPct' percentile of a unit's
%   IBI distribution are masked as NaN to distinguish "undefined" states
%   from zero values.
%

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'brst', @isstruct);
addRequired(p, 'spktimes', @iscell);
addParameter(p, 'binSize', 60, @isnumeric);
addParameter(p, 'kernelSD', 300, @isnumeric);
addParameter(p, 'smoothWin', 10, @isnumeric);
addParameter(p, 'ibiPct', 95, @isnumeric);
addParameter(p, 'flgPlot', true, @islogical);
addParameter(p, 'flgSave', false, @islogical);
addParameter(p, 'basepath', pwd, @ischar);

parse(p, brst, spktimes, varargin{:});
binSize   = p.Results.binSize;
kernelSD  = p.Results.kernelSD;
smoothWin = p.Results.smoothWin;
ibiPct    = p.Results.ibiPct;
flgPlot   = p.Results.flgPlot;
flgSave   = p.Results.flgSave;
basepath  = p.Results.basepath;


%% ========================================================================
%  INITIALIZE
%  ========================================================================

nUnits = length(brst.all);

% Determine global time range
maxTime = 0;
for iUnit = 1:nUnits
    if ~isempty(spktimes{iUnit})
        maxTime = max(maxTime, max(spktimes{iUnit}));
    end
end

% Time vector
t = 0 : binSize : ceil(maxTime);
nTime = length(t);

% Initialize Matrices (nUnits x nTime)
dyn.time  = t;
dyn.rate  = zeros(nUnits, nTime);
dyn.bfrac = zeros(nUnits, nTime);
dyn.dur   = nan(nUnits, nTime);
dyn.nspks = nan(nUnits, nTime);
dyn.freq  = nan(nUnits, nTime);
dyn.ibi   = nan(nUnits, nTime);


%% ========================================================================
%  COMPUTE DYNAMICS
%  ========================================================================

% Create Gaussian Kernel
kRange = -3*kernelSD : binSize : 3*kernelSD;
kernel = normpdf(kRange, 0, kernelSD);
kernel = kernel / sum(kernel);

for iUnit = 1:nUnits

    st = spktimes{iUnit};
    b_struct = brst.all{iUnit};

    if isempty(st) || isempty(b_struct) || isempty(b_struct.times)
        continue;
    end

    % ---------------------------------------------------------------------
    % Density Estimation (Continuous)
    % ---------------------------------------------------------------------

    % Burst Rate
    onsets = b_struct.times(:, 1);
    bCounts = histcounts(onsets, [t, inf]);
    % conv produces row vector if bCounts is row
    rawRate = conv(bCounts, kernel, 'same') / binSize;
    dyn.rate(iUnit, :) = rawRate;

    % Burst Fraction
    % Identify burst spikes
    isBurstSpike = false(size(st));
    for b = 1:size(b_struct.times, 1)
        isBurstSpike = isBurstSpike | ...
            (st >= b_struct.times(b, 1) & st <= b_struct.times(b, 2));
    end
    st_burst = st(isBurstSpike);

    allCounts = histcounts(st, [t, inf]);
    densAll = conv(allCounts, kernel, 'same') / binSize;
    brstCounts = histcounts(st_burst, [t, inf]);
    densBrst = conv(brstCounts, kernel, 'same') / binSize;

    frac = zeros(size(densAll));
    mask = densAll > 1e-6;
    frac(mask) = densBrst(mask) ./ densAll(mask);
    dyn.bfrac(iUnit, :) = frac;


    % ---------------------------------------------------------------------
    % Property Tracking (Interpolated & Masked)
    % ---------------------------------------------------------------------

    % Calculate Mask Threshold (Adaptive)
    unit_ibis = b_struct.ibi(~isnan(b_struct.ibi));
    if isempty(unit_ibis)
        maskThresh = 60;
    else
        maskThresh = prctile(unit_ibis, ibiPct);
    end

    % Ensure threshold isn't too small (min 2*binSize)
    maskThresh = max(maskThresh, 2 * binSize);

    % Apply to all props
    tTimes = b_struct.times(:, 1);

    if length(tTimes) < 2
        continue
    end

    % Process Properties
    dyn.dur(iUnit, :)   = proc_prop(b_struct.dur, tTimes, t, smoothWin, maskThresh);
    dyn.nspks(iUnit, :) = proc_prop(b_struct.nspks, tTimes, t, smoothWin, maskThresh);
    dyn.freq(iUnit, :)  = proc_prop(b_struct.freq, tTimes, t, smoothWin, maskThresh);
    dyn.ibi(iUnit, :)   = proc_prop(b_struct.ibi, tTimes, t, smoothWin, maskThresh);

end


%% ========================================================================
%  PLOTTING
%  ========================================================================

if flgPlot
    % Convert to table for tblGUI_xy
    % rmfield removes 'time'
    % struct2table creates table where fields are variables.
    % If fields are matrices (nUnits x nTime), table variables will be same.
    tbl = struct2table(rmfield(dyn, 'time'));

    % Plot
    % If tblGUI_xy accepts yVar optional, we can default to 'rate' or let user pick.
    % We launch it with 'rate' as default.
    tblGUI_xy(dyn.time, tbl, 'yVar', 'rate');
end




%% ========================================================================
%  SAVE
%  ========================================================================

if flgSave
    [~, basename] = fileparts(basepath);
    dynfile = fullfile(basepath, [basename, '.brstDyn.mat']);
    brstDyn = dyn; %#ok<NASGU>
    save(dynfile, 'brstDyn');
end

end     % EOF


%% ========================================================================
%  NESTED FUNCTIONS
%  ========================================================================

function vec = proc_prop(rawVal, timePoints, tVec, smoothWin, maskThresh)
% Processes a property vector: Smooth -> Interpolate -> Mask

if length(rawVal) < smoothWin
    % Not enough points to smooth well, just use mean or raw
    smVal = rawVal;
else
    smVal = movmedian(rawVal, smoothWin, 'omitnan');
end

% Interpolate to global time grid
% 'linear' provides smooth transitions between events
% 'extrap' with NaN prevents wild guesses outside range
vec = interp1(timePoints, smVal, tVec, 'linear', NaN);

% Identify Gaps (Masking)
% Manual Gap Filling approach:

% Gap before first burst
vec(tVec < timePoints(1) - maskThresh) = NaN;

% Gap after last burst
vec(tVec > timePoints(end) + maskThresh) = NaN;

% Interior Gaps
diffs = diff(timePoints);
gaps = find(diffs > maskThresh);
for g = 1:length(gaps)
    idx1 = gaps(g);
    idx2 = gaps(g) + 1;
    t_start = timePoints(idx1) + maskThresh;
    t_end   = timePoints(idx2) - maskThresh;

    if t_end > t_start
        vec(tVec > t_start & tVec < t_end) = NaN;
    end
end

% Ensure row vector for output matrix
vec = vec(:)';
end


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
%  To avoid arbitrary thresholds (e.g., "2 minutes silence"), we determine
%  the masking window based on the statistical behavior of the Inter-Burst
%  Interval (IBI). The function calculates the 95th (or user-specified)
%  percentile of the IBI for each unit. If the time elapsed since the last
%  detected burst exceeds this calculated "Expected IBI," the property is
%  masked as NaN. This ensures that the population average is only
%  calculated from units that are actively participating in the network's
%  bursting state at that specific time.
%  ========================================================================

%% ========================================================================
%  NOTE: BEHAVIOR OF MOVING MEDIANS AND MISSING DATA
%  ========================================================================
%  Applying a movmedian across the entire length of the recording is only
%  possible if you have a value for every time bin, which effectively
%  forces you to assume a value (like zero) for bins without bursts. This
%  is incorrect for structural metrics like brstDur or nspks. By performing
%  the movmedian on the event index first—meaning we smooth across the last
%  N bursts regardless of when they happened—we preserve the internal
%  structural integrity of the burst properties. The interpolation and
%  subsequent masking then allow us to align these "active state" values
%  across the population while correctly identifying periods where the unit
%  has dropped out of the bursting state.
%  ========================================================================
