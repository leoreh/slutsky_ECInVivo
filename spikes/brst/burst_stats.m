function stats = burst_stats(burst, spktimes, varargin)
% BURST_STATS Calculates burst statistics per unit for each time window.
%
%   stats = BURST_STATS(BRST, SPKTIMES, ...) calculates summary statistics
%   (e.g., burst rate, duration, spikes per burst) for bursts detected by
%   burst_detect, within specified time windows.
%
%   INPUTS:
%       burst        - (struct) Output from burst_detect.m (must contain .times, etc.)
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
%                     .bN         : Number of bursts
%                     .br  : Burst event rate (Hz) (Count / Window Duration)
%                     .fr      : Total firing rate (Hz)
%                     .frBurst     : Burst spike firing rate (Hz)
%                     .frSingle     : Single spike firing rate (Hz)
%                     .dur        : Mean burst duration (s)
%                     .freq       : Mean intra-burst frequency (Hz)
%                     .ibi        : Mean inter-burst interval (s)
%                     .bSize      : Mean spikes per burst
%                     .pBurst      : Probability of spikes in bursts (0-1)
%                     .winCalc    : The time windows used [nWin x 2]
%
%   NOTES:
%       - Bursts are assigned to a window based on their START time.
%       - IBI statistics for a window are the mean of the IBIs of bursts
%         starting in that window. (IBI is the interval preceding the burst).
%       - If no bursts occur in a window, count/rate/pBspk are 0, others NaN.
%
%   See also: BURST_DETECT, BURST_DYNAMICS

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'burst', @isstruct);
addRequired(p, 'spktimes', @iscell);
addParameter(p, 'winCalc', [], @isnumeric);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgSave', false, @islogical);

parse(p, burst, spktimes, varargin{:});
winCalc  = p.Results.winCalc;
basepath = p.Results.basepath;
flgSave  = p.Results.flgSave;


%% ========================================================================
%  INITIALIZE
%  ========================================================================

nUnits = length(burst.times);

% Handle winCalc
if isempty(winCalc)
    % Determine max time from spktimes
    maxTime = max(cellfun(@(x) max([0; x(:)]), spktimes));
    winCalc = [0, maxTime];
end

nWin = size(winCalc, 1);

% Initialize Output Matrices [nUnits x nWin]
stats.bN        = zeros(nUnits, nWin);
stats.br = zeros(nUnits, nWin);
stats.fr     = zeros(nUnits, nWin);
stats.frBurst    = zeros(nUnits, nWin);
stats.frSingle    = zeros(nUnits, nWin);
stats.pBurst     = zeros(nUnits, nWin);
stats.bSize     = nan(nUnits, nWin);
stats.dur       = nan(nUnits, nWin);
stats.freq      = nan(nUnits, nWin);
stats.ibi       = nan(nUnits, nWin);

% Info
stats.info.input   = p.Results;
stats.info.winCalc = winCalc;


%% ========================================================================
%  COMPUTE LOOP
%  ========================================================================

for iUnit = 1:nUnits

    % Access burst properties
    times = burst.times{iUnit};
    nBspk = burst.size{iUnit};
    dur   = burst.dur{iUnit};
    freq  = burst.freq{iUnit};
    ibi   = burst.ibi{iUnit};
    st = spktimes{iUnit};

    % Count bursts fully contained in window.
    if isempty(times)
        bStart = [];
        bEnd   = [];
    else
        bStart = times(:, 1);
        bEnd   = times(:, 2);
    end

    for iWin = 1:nWin
        wStart = winCalc(iWin, 1);
        wEnd   = winCalc(iWin, 2);
        wDur   = wEnd - wStart;

        % Count Bursts Fully Contained in Window
        bIdx = (bStart >= wStart) & (bEnd <= wEnd);

        % Floor of detection (1 event per window)
        c = 1 / wDur;

        % Count & Event Rate
        nb = sum(bIdx);
        stats.bN(iUnit, iWin) = nb;
        stats.br(iUnit, iWin) = (nb / wDur) + c;

        % Structural Means
        if nb > 0
            stats.bSize(iUnit, iWin)   = mean(nBspk(bIdx));
            stats.dur(iUnit, iWin)     = mean(dur(bIdx));
            stats.freq(iUnit, iWin)    = mean(freq(bIdx));
            stats.ibi(iUnit, iWin)     = mean(ibi(bIdx), 'omitnan');
        end

        % Firing Rates & Partitioning
        % Spikes in window
        stIdx = (st >= wStart) & (st <= wEnd);
        nst = sum(stIdx);

        frRawTot = nst / wDur;
        stats.fr(iUnit, iWin) = frRawTot + c;

        if nst > 0
            % Count spikes from bursts fully contained in window
            % Since burst is fully contained, all its spikes are in window.
            bSpks = sum(nBspk(bIdx));
            frBspk = bSpks / wDur;

            stats.frBurst(iUnit, iWin) = frBspk + c;
            stats.pBurst(iUnit, iWin)  = bSpks / nst;
        else
            frBspk = 0;
            stats.frBurst(iUnit, iWin) = 0 + c;
            stats.pBurst(iUnit, iWin)  = NaN;
        end

        % Single Spike Rate
        stats.frSingle(iUnit, iWin) = (frRawTot - frBspk) + c;

    end
end


%% ========================================================================
%  SAVE
%  ========================================================================

if flgSave
    [~, basename] = fileparts(basepath);
    fname = fullfile(basepath, [basename, '.burstStats.mat']);
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