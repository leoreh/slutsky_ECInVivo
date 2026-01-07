function stats = brst_stats(brst, spktimes, varargin)
% BRST_STATS Calculates burst statistics per unit for each time window.
%
%   stats = BRST_STATS(BRST, SPKTIMES, ...) calculates summary statistics
%   (e.g., burst rate, duration, spikes per burst) for bursts detected by
%   brst_detect, within specified time windows.
%
%   INPUTS:
%       brst        - (struct) Output from brst_detect.m (must contain .times, etc.)
%       spktimes    - (cell) Spike times per unit (e.g., {unit1, unit2}).
%       varargin    - (param/value) Optional parameters:
%                     'winCalc'  : (num) [M x 2] matrix of time windows.
%                                  Default: [0, max(spktimes)].
%                     'basepath' : (char) Base path for saving {pwd}
%                     'flgSave'  : (log) Save result as stats struct {false}
%
%   OUTPUTS:
%       stats       - (struct) Burst statistics structure.
%                     Fields are matrices of size [nUnits x nWin]:
%                     .nb       : Number of bursts
%                     .rate     : Burst rate (Hz) (Count / Window Duration)
%                     .dur      : Mean burst duration (s)
%                     .freq     : Mean intra-burst frequency (Hz)
%                     .ibi      : Mean inter-burst interval (s)
%                     .nBspk    : Mean spikes per burst
%                     .pBspk    : Porbability of spikes in bursts (0-1)
%                     .winCalc  : The time windows used [nWin x 2]
%
%   NOTES:
%       - Bursts are assigned to a window based on their START time.
%       - IBI statistics for a window are the mean of the IBIs of bursts
%         starting in that window. (IBI is the interval preceding the burst).
%       - If no bursts occur in a window, count/rate/pBspk are 0, others NaN.
%
%   See also: BRST_DETECT, BRST_DYNAMICS

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'brst', @isstruct);
addRequired(p, 'spktimes', @iscell);
addParameter(p, 'winCalc', [], @isnumeric);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgSave', false, @islogical);

parse(p, brst, spktimes, varargin{:});
winCalc  = p.Results.winCalc;
basepath = p.Results.basepath;
flgSave  = p.Results.flgSave;


%% ========================================================================
%  INITIALIZE
%  ========================================================================

nUnits = length(brst.times);

% Handle winCalc
if isempty(winCalc)
    % Determine max time from spktimes
    maxTime = max(cellfun(@(x) max([0; x(:)]), spktimes));
    winCalc = [0, maxTime];
end

nWin = size(winCalc, 1);

% Initialize Output Matrices [nUnits x nWin]
stats.nb      = zeros(nUnits, nWin);
stats.rate    = zeros(nUnits, nWin);
stats.nBspk   = nan(nUnits, nWin);
stats.dur = nan(nUnits, nWin);
stats.freq    = nan(nUnits, nWin);
stats.ibi     = nan(nUnits, nWin);
stats.pBspk   = zeros(nUnits, nWin);

% Info
stats.info.input   = p.Results;
stats.info.winCalc = winCalc;


%% ========================================================================
%  COMPUTE LOOP
%  ========================================================================

for iUnit = 1:nUnits

    % Access burst properties
    times = brst.times{iUnit};
    nBspk = brst.nBspk{iUnit};
    dur   = brst.dur{iUnit};
    freq  = brst.freq{iUnit};
    ibi   = brst.ibi{iUnit};

    % Check if unit has spikes
    st = spktimes{iUnit};
    if isempty(st) || isempty(times)
        continue;
    end

    bStart = times(:, 1);

    for iWin = 1:nWin
        wStart = winCalc(iWin, 1);
        wEnd   = winCalc(iWin, 2);
        wDur   = wEnd - wStart;

        bIdx = (bStart >= wStart) & (bStart <= wEnd);
        if ~any(bIdx)
            continue;
        end

        % Count & Rate
        nb = sum(bIdx);
        stats.nb(iUnit, iWin) = nb;
        stats.rate(iUnit, iWin) = nb / wDur;

        % Structural Means
        stats.nBspk(iUnit, iWin)   = mean(nBspk(bIdx));
        stats.dur(iUnit, iWin)     = mean(dur(bIdx));
        stats.freq(iUnit, iWin)    = mean(freq(bIdx));
        stats.ibi(iUnit, iWin)     = mean(ibi(bIdx), 'omitnan');

        % Burst Spike Fraction
        % Spikes in window
        stIdx = (st >= wStart) & (st <= wEnd);
        nst = sum(stIdx);

        if nst > 0
            bSpks = sum(nBspk(bIdx));
            stats.pBspk(iUnit, iWin) = bSpks / nst;
        else
            stats.pBspk(iUnit, iWin) = NaN;
        end

    end
end


%% ========================================================================
%  SAVE
%  ========================================================================

if flgSave
    [~, basename] = fileparts(basepath);
    fname = fullfile(basepath, [basename, '.brstStats.mat']);
    save(fname, 'stats');
end

end     % EOF


%% ========================================================================
%  NOTE: PREDICTED VALUE FOR RATE NORMALIZATION
%  ========================================================================
%  Burst metrics are intrinsically linked to the underlying firing rate
%  because the probability of inter-spike intervals (ISIs) falling below a
%  threshold increases with spike density. Simple division (ratio
%  normalization) is often insufficient because the relationship between
%  rate and burstiness is rarely proportional; it typically features
%  non-zero intercepts and non-linearities. The Predicted Value approach,
%  as detailed by Eisenman et al. (2015), provides a robust empirical
%  framework to isolate physiological spike- pattern changes from these
%  activity-driven artifacts.
%
%  THEORY AND IMPLEMENTATION
%  The method establishes a "physiological map" by fitting regression
%  models to baseline (WT/Control) data across a wide range of activity
%  levels. This regression defines the "expected" burst
%  parameter for any given firing rate. During experimental phases
%  where rates may shift, the "Predicted Value" is calculated for each
%  unit based on its current rate. The final metric is reported as a
%  percentage of this prediction, where 100% signifies a spike pattern
%  identical to a WT unit at that specific rate.
%
%  By utilizing this method, a researcher can conclude that a change in
%  bursting is a fundamental alteration of the neuron's signaling
%  strategy rather than a trivial consequence of an increased or
%  decreased firing rate.
%  ========================================================================

%% ========================================================================
%  NOTE: NETWORK DRIVE VS. UNIT ACTIVITY
%  ========================================================================
%  A critical decision in burst normalization is selecting the reference
%  rate for the x-axis: the individual unit's firing rate or the Array-
%  Wide Spike Detection Rate (ASDR). In highly synchronized preparations
%  like hippocampal cultures, bursts are seldom isolated events; they are
%  network-wide phenomena driven by the collective excitatory tone of the
%  population.
%
%  ADVANTAGES OF ASDR (NETWORK RATE)
%  The ASDR serves as a proxy for the total "network drive" or synaptic
%  pressure experienced by every neuron in the culture.
%  Using ASDR as the independent variable for normalization is
%  statistically more stable than individual rates because it averages
%  out the high-frequency noise and sorting artifacts inherent to
%  single-unit detection. It ensures that the "map"
%  represents the network's state rather than a single cell's volatility.
%
%  * Reflects global excitatory drive.
%  * Reduces single-unit signal noise.
%  * Accounts for network synchrony.
%  * Stabilizes normalization reference frames.
%
%  While individual FR is useful for analyzing autonomous firing
%  properties, the ASDR provides a more accurate physiological context
%  for understanding how a neuron's burstiness scales within a
%  communicating network.
%  ========================================================================

%% ========================================================================
%  NOTE: STATISTICAL CONTROL VIA FIRING RATE COVARIATES IN LME
%  ========================================================================
%  Linear Mixed-Effects (LME) models allow for the integration of
%  activity-level corrections directly into the primary statistical
%  analysis. By including Firing Rate (FR) as a continuous covariate
%  in the model formula (e.g., BSpks ~ Group + FR + (1|Name)), the
%  model mathematically partitions the variance in burstiness.
%  This process "partials out" the effect of activity, effectively
%  comparing the genotypes as if they were firing at the same rate.
%
%  INTERPRETATION AND INTERACTION
%  This approach is particularly valuable for the MCU-KO vs. WT
%  comparison. If the genotype effect remains significant after
%  including FR as a covariate, the difference in burstiness is
%  statistically independent of activity levels. Furthermore,
%  implementing an interaction term (Group * FR) allows the researcher
%  to test if the relationship between rate and burstiness itself has
%  changed. For example, if KO neurons become "disproportionately
%  bursty" as they recover activity.
%
%  * Partials out rate-driven variance.
%  * Adjusts means to common rates.
%  * Detects genotype-rate interaction effects.
%  * Preserves unit-level variability.
%
%  This method is an essential adjunct to the Predicted Value method,
%  offering a formalized p-value for genotype differences that is
%  rigorously controlled for the confounding influence of firing
%  frequency.
%  ========================================================================