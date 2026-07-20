function stats = burst_stats(burst, spktimes, varargin)
% BURST_STATS Calculates burst statistics per unit for each condition.
%
%   stats = BURST_STATS(burst, spktimes, ...)
%
%   SUMMARY:
%       Summary statistics (burst rate, duration, spikes per burst, the
%       burst / single decomposition of firing rate) for bursts detected by
%       burst_detect, restricted to time windows.
%
%       Each output column is one condition. By default every row of winCalc
%       is a condition of its own, which is how the MEA pipeline gets its
%       BSL / Acute / SS columns. With flgPool the rows are instead the bouts
%       of a single condition (e.g. every NREM bout) and collapse to one
%       column - that is the form spk_byCond calls, since it supplies the
%       label and stacks the conditions as rows.
%
%   INPUTS:
%       burst    - (Struct) Output of burst_detect (.times, .size, .dur,
%                           .freq, .ibi).
%       spktimes - (Cell)   {nUnits x 1} of spike times [s].
%       varargin - Parameter/Value pairs:
%           'winCalc'  - (Mat)  [nWin x 2] time windows [s].
%                               {[0, max(spktimes)]}
%           'flgPool'  - (Log)  Treat every row of winCalc as one condition
%                               and return a single column. {false}
%           'basepath' - (Char) Save location. {pwd}
%           'flgSave'  - (Log)  Save <basename>.burstStats.mat. {false}
%
%   OUTPUTS:
%       stats    - (Struct) Fields are [nUnits x nCol]:
%                    .bN       Number of bursts
%                    .br       Burst event rate [Hz]
%                    .fr       Total firing rate [Hz]
%                    .frBurst  Burst spike firing rate [Hz]
%                    .frSingle Single spike firing rate [Hz]
%                    .dur      Mean burst duration [s]
%                    .freq     Mean intra-burst frequency [Hz]
%                    .ibi      Mean inter-burst interval [s]
%                    .bSize    Mean spikes per burst
%                    .pBurst   Fraction of spikes in bursts
%
%   NOTES:
%       - A burst counts toward a condition only if it lies wholly inside a
%         single window of it, so a burst straddling a bout edge is dropped
%         rather than split.
%       - IBI is the gap preceding a burst. It belongs to the condition only
%         when the previous burst sits in the same window; otherwise the gap
%         spans time the condition excludes. The first burst of a window
%         therefore never contributes an IBI.
%       - Rates carry a 1/wDur floor, the rate of a single event over the
%         window, so a silent unit is bounded away from zero on a log scale.
%
%   DEPENDENCIES:
%       intervals, backup_file.
%
%   HISTORY:
%       260719    flgPool added; ibi restricted to same-window pairs.
%
%   See also: BURST_DETECT, BURST_DYNAMICS, SPK_BYCOND

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'burst', @isstruct);
addRequired(p, 'spktimes', @iscell);
addParameter(p, 'winCalc', [], @isnumeric);
addParameter(p, 'flgPool', false, @islogical);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgSave', false, @islogical);

parse(p, burst, spktimes, varargin{:});
winCalc  = p.Results.winCalc;
flgPool  = p.Results.flgPool;
basepath = p.Results.basepath;
flgSave  = p.Results.flgSave;


%% ========================================================================
%  INITIALIZE
%  ========================================================================

nUnits = length(burst.times);

if isempty(winCalc)
    maxTime = max(cellfun(@(x) max([0; x(:)]), spktimes));
    winCalc = [0, maxTime];
end

% pooled bouts arrive unordered and may touch; consolidating them keeps the
% edge vector below strictly monotonic, which is what discretize needs
if flgPool
    winCalc = winCalc(winCalc(:, 2) > winCalc(:, 1), :);
    winCalc = intervals(winCalc);
    winCalc = winCalc.consolidate();
    winCalc = winCalc.ints;
end

nWin = size(winCalc, 1);
if flgPool
    winGrp = {1 : nWin};
else
    winGrp = num2cell(1 : nWin);
end
nCol = numel(winGrp);

% Initialize output matrices [nUnits x nCol]
stats.bN       = zeros(nUnits, nCol);
stats.br       = zeros(nUnits, nCol);
stats.fr       = zeros(nUnits, nCol);
stats.frBurst  = zeros(nUnits, nCol);
stats.frSingle = zeros(nUnits, nCol);
stats.pBurst   = zeros(nUnits, nCol);
stats.bSize    = nan(nUnits, nCol);
stats.dur      = nan(nUnits, nCol);
stats.freq     = nan(nUnits, nCol);
stats.ibi      = nan(nUnits, nCol);

% Info
stats.info.input   = p.Results;
stats.info.winCalc = winCalc;
stats.info.flgPool = flgPool;


%% ========================================================================
%  COMPUTE LOOP
%  ========================================================================

for iUnit = 1 : nUnits

    % Access burst properties
    times = burst.times{iUnit};
    nBspk = burst.size{iUnit};
    dur   = burst.dur{iUnit};
    freq  = burst.freq{iUnit};
    ibi   = burst.ibi{iUnit};
    st    = spktimes{iUnit};

    if isempty(times)
        bStart = [];
        bEnd   = [];
    else
        bStart = times(:, 1);
        bEnd   = times(:, 2);
    end

    for iCol = 1 : nCol

        win   = winCalc(winGrp{iCol}, :);
        wDur  = sum(win(:, 2) - win(:, 1));
        edges = reshape(win', [], 1);

        % membership by edge search: an odd bin index means inside a window
        % of this condition, an even one the gap between two
        iB1 = discretize(bStart, edges);
        iB2 = discretize(bEnd, edges);
        bIdx = ~isnan(iB1) & iB1 == iB2 & mod(iB1, 2) == 1;

        % Floor of detection (1 event per window)
        c = 1 / wDur;

        % Count & event rate
        nb = sum(bIdx);
        stats.bN(iUnit, iCol) = nb;
        stats.br(iUnit, iCol) = (nb / wDur) + c;

        % Structural means
        if nb > 0
            stats.bSize(iUnit, iCol) = mean(nBspk(bIdx));
            stats.dur(iUnit, iCol)   = mean(dur(bIdx));
            stats.freq(iUnit, iCol)  = mean(freq(bIdx));

            % the preceding burst must share the window, else the gap
            % crosses time this condition excludes
            iPrev = [NaN; iB1(1 : end - 1)];
            ibiIdx = bIdx & iPrev == iB1;
            stats.ibi(iUnit, iCol) = mean(ibi(ibiIdx), 'omitnan');
        end

        % Firing rates & partitioning
        iS = discretize(st, edges);
        nst = sum(~isnan(iS) & mod(iS, 2) == 1);

        frRawTot = nst / wDur;
        stats.fr(iUnit, iCol) = frRawTot + c;

        if nst > 0
            % a burst is wholly inside the window, so all of its spikes are
            bSpks = sum(nBspk(bIdx));
            frBspk = bSpks / wDur;

            stats.frBurst(iUnit, iCol) = frBspk + c;
            stats.pBurst(iUnit, iCol)  = bSpks / nst;
        else
            frBspk = 0;
            stats.frBurst(iUnit, iCol) = 0 + c;
            stats.pBurst(iUnit, iCol)  = NaN;
        end

        % Single spike rate
        stats.frSingle(iUnit, iCol) = (frRawTot - frBspk) + c;
    end
end


%% ========================================================================
%  SAVE
%  ========================================================================

if flgSave
    [~, basename] = fileparts(basepath);
    fname = fullfile(basepath, [basename, '.burstStats.mat']);
    backup_file(fname);
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
