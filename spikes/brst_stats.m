function stats = brst_stats(brst, spktimes, varargin)
% BRST_STATS Calculates burst statistics per unit for each time window.
%
%   stats = BRST_STATS(BRST, SPKTIMES, ...) calculates summary statistics
%   (e.g., burst rate, duration, spikes per burst) for bursts detected by
%   brst_detect, within specified time windows.
%
%   INPUTS:
%       brst        - (struct) Output from brst_detect.m (must contain .all)
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
%                     .detect   : Number of bursts
%                     .rate     : Burst rate (Hz) (Count / Window Duration)
%                     .nspks    : Mean spikes per burst
%                     .brstDur  : Mean burst duration (s)
%                     .freq     : Mean intra-burst frequency (Hz)
%                     .ibi      : Mean inter-burst interval (s)
%                     .bspks    : Fraction of spikes in bursts (0-1)
%                     .winCalc  : The time windows used [nWin x 2]
%
%   NOTES:
%       - Bursts are assigned to a window based on their START time.
%       - IBI statistics for a window are the mean of the IBIs of bursts
%         starting in that window. (IBI is the interval preceding the burst).
%       - If no bursts occur in a window, count/rate/bspks are 0, others NaN.
%

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
%  PREPARATIONS
%  ========================================================================

nUnits = length(brst.all);

% Handle winCalc
if isempty(winCalc)
    % Determine max time from spktimes
    maxTime = 0;
    for i = 1:length(spktimes)
        if ~isempty(spktimes{i})
            maxTime = max(maxTime, max(spktimes{i}));
        end
    end
    winCalc = [0, maxTime];
end

nWin = size(winCalc, 1);

% Initialize Output Matrices [nUnits x nWin]
stats.detect  = zeros(nUnits, nWin);
stats.rate    = zeros(nUnits, nWin);
stats.nspks   = nan(nUnits, nWin);
stats.brstDur = nan(nUnits, nWin);
stats.freq    = nan(nUnits, nWin);
stats.ibi     = nan(nUnits, nWin);
stats.bspks   = zeros(nUnits, nWin);

% Info
stats.info.runtime = datetime("now");
stats.info.input   = p.Results;
stats.info.winCalc = winCalc;


%% ========================================================================
%  CALCULATE STATISTICS
%  ========================================================================

for iUnit = 1:nUnits

    b_struct = brst.all{iUnit};
    st = spktimes{iUnit};

    if isempty(st)
        continue;
        % stats initialized to 0/NaN already
    end

    % If no bursts detected for unit, b_struct might be empty or fields empty
    if isempty(b_struct) || isempty(b_struct.times)
        % No bursts detected at all
        % detect=0, rate=0, bspks=0, others NaN. (Already set)
        continue;
    end

    % Prepare burst properties
    bStarts = b_struct.times(:, 1);
    bNspks  = b_struct.nspks;
    bDur    = b_struct.dur;
    bFreq   = b_struct.freq;
    bIbi    = b_struct.ibi;

    for iWin = 1:nWin
        tStart = winCalc(iWin, 1);
        tEnd   = winCalc(iWin, 2);
        wDur   = tEnd - tStart;

        % Filter bursts in this window (based on Start Time)
        % Using >= tStart and < tEnd (standard binning)
        % Check strict inequality handling? User didn't specify.
        % Using [Start, End] inclusive logic usually safer for single windows
        % but standard is [ ).
        % Defaulting to: Start >= tStart & Start <= tEnd

        idx = (bStarts >= tStart) & (bStarts <= tEnd);

        if ~any(idx)
            % No bursts in this window
            continue;
        end

        % Count & Rate
        nb = sum(idx);
        stats.detect(iUnit, iWin) = nb;
        if wDur > 0
            stats.rate(iUnit, iWin) = nb / wDur;
        else
            stats.rate(iUnit, iWin) = NaN;
        end

        % Structural Means
        stats.nspks(iUnit, iWin)   = mean(bNspks(idx));
        stats.brstDur(iUnit, iWin) = mean(bDur(idx));
        stats.freq(iUnit, iWin)    = mean(bFreq(idx));
        stats.ibi(iUnit, iWin)     = mean(bIbi(idx), 'omitnan');

        % Burst Spike Fraction
        % 1. Count total spikes in window
        % 2. Count burst spikes in window (approximate or exact?)
        % bspks definition in brst_detect was "total_spikes_in_bursts / length(st)".
        % Here it should be "spikes in bursts in window" / "total spikes in window".

        % Spikes in window
        st_idx = (st >= tStart) & (st <= tEnd);
        nStWin = sum(st_idx);

        if nStWin > 0
            % Bursts contributing to this window:
            % If we define "bursts in window" by start time, we should probably
            % sum the spikes OF THOSE BURSTS.
            % But technically some spikes of a burst starting in window might be outside?
            % OR: Count spikes that fall in window AND belong to a burst?
            % Code in brst_detect: sum(b_struct.nspks) / length(st).
            % This implies using the spike count of the bursts *identified*.
            % I will stick to "filtered bursts" logic for consistency with other stats.

            total_bSpks = sum(bNspks(idx));
            stats.bspks(iUnit, iWin) = total_bSpks / nStWin;
        else
            stats.bspks(iUnit, iWin) = 0; % Or NaN? If no spikes, fraction is undefined.
            % But detect=0 implies 0 fraction.
            % If nStWin=0, then nb must be 0 (since bursts require spikes).
            % So 0/0 -> 0 is reasonable or NaN.
            stats.bspks(iUnit, iWin) = NaN;
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
%  * Maps expected physiological states.
%  * Corrects for non-proportional scaling.
%  * Isolates rate-independent genotype effects.
%  * Benchmarks against WT population.
%  * Validates drug-induced pattern shifts.
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