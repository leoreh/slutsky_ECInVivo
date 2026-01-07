function brst = brst_detect(spktimes, varargin)
% BRST_DETECT Detects bursts using the Max Interval method.
%
%   brst = BRST_DETECT(SPKTIMES, ...) implements the Max Interval burst
%   detection algorithm, defined by fixed thresholds for inter-spike
%   intervals (ISIs) and burst properties.
%
%   Default parameters were extracted by visualizing the isi histogram from
%   WT MEA recordings (December 25, see brst_isiValley.m)
%
%   INPUTS:
%       spktimes    - (cell) Spike times per unit (e.g., {unit1, unit2}).
%                     Times should be in seconds.
%       varargin    - (param/value) Optional parameters:
%                     'isiStart'  : (num) Max ISI to start burst {0.013}
%                     'isiEnd'    : (num) Max ISI within burst {0.02}
%                     'minIbi'    : (num) Min Inter-Burst Interval {0.05}
%                     'minDur'    : (num) Min burst duration {0.015}
%                     'minSpks'   : (num) Min spikes in burst {3}
%                     'flgPlot'   : (log) Plot raster with bursts {false}
%                     'basepath'  : (char) Recording path {pwd}
%                     'flgSave'   : (log) Save to file {true}
%                     'flgForce'  : (log) Analyze even if exists {false}
%
%   OUTPUTS:
%       brst        - (struct) Burst event data.
%                     .times    : (nUnits x 1 cell) Burst start/end times
%                     .nBspk    : (nUnits x 1 cell) Count of spikes per burst
%                     .dur      : (nUnits x 1 cell) Duration of bursts (s)
%                     .freq     : (nUnits x 1 cell) Intra-burst frequency (Hz)
%                     .ibi      : (nUnits x 1 cell) Pre-burst interval (s)
%                     .spktimes : (nUnits x 1 cell) Spike times within bursts
%                     .params   : (struct) Parameters used
%
%   See also: BRST_DYNAMICS, BRST_STATS

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'spktimes', @iscell);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'isiStart', 0.013, @isnumeric);
addParameter(p, 'isiEnd', 0.02, @isnumeric);
addParameter(p, 'minIbi', 0.05, @isnumeric);
addParameter(p, 'minDur', 0.015, @isnumeric);
addParameter(p, 'minSpks', 3, @isnumeric);
addParameter(p, 'flgPlot', false, @islogical);
addParameter(p, 'flgSave', true, @islogical);
addParameter(p, 'flgForce', false, @islogical);

parse(p, spktimes, varargin{:});
basepath  = p.Results.basepath;
isiStart  = p.Results.isiStart;
isiEnd    = p.Results.isiEnd;
minIbi    = p.Results.minIbi;
minDur    = p.Results.minDur;
minSpks   = p.Results.minSpks;
flgPlot   = p.Results.flgPlot;
flgSave   = p.Results.flgSave;
flgForce  = p.Results.flgForce;
params    = p.Results;


%% ========================================================================
%  INITIALIZE
%  ========================================================================

% Check existence
[~, basename] = fileparts(basepath);
saveFile = fullfile(basepath, [basename, '.brst.mat']);

if exist(saveFile, 'file') && ~flgForce && ~flgPlot
    load(saveFile, 'brst');
    return;
end

nUnits = length(spktimes);

% Initialize Output
brst.times    = cell(nUnits, 1);
brst.nBspk    = cell(nUnits, 1);
brst.dur      = cell(nUnits, 1);
brst.freq     = cell(nUnits, 1);
brst.ibi      = cell(nUnits, 1);
brst.spktimes = cell(nUnits, 1);
brst.params   = params;


%% ========================================================================
%  COMPUTE LOOP
%  ========================================================================

for iUnit = 1 : nUnits

    st = spktimes{iUnit};
    if isempty(st) || length(st) < minSpks
        continue;
    end

    % Ensure column vector
    st = st(:);
    isi = diff(st);
    nBspk = length(st);

    % ---------------------------------------------------------------------
    % Burst Detection
    % ---------------------------------------------------------------------

    bspkIdx = [];
    inBurst = false;
    iStart = -1;

    iSpk = 1;
    while iSpk < nBspk
        if ~inBurst
            % Start criteria
            if isi(iSpk) <= isiStart
                inBurst = true;
                iStart = iSpk;
                iSpk = iSpk + 1;
            else
                iSpk = iSpk + 1;
            end
        else
            % Continuation criteria
            if isi(iSpk) <= isiEnd
                iSpk = iSpk + 1;
            else
                % End of burst
                bspkIdx = [bspkIdx; iStart, iSpk]; 
                inBurst = false;
                iSpk = iSpk + 1;
            end
        end
    end

    % Close last burst if active
    if inBurst
        bspkIdx = [bspkIdx; iStart, nBspk]; 
    end

    if isempty(bspkIdx)
        continue;
    end

    % ---------------------------------------------------------------------
    % Merge Bursts
    % ---------------------------------------------------------------------

    mrgIdx = [];
    prevStart = bspkIdx(1, 1);
    prevEnd   = bspkIdx(1, 2);
    for ib = 2:size(bspkIdx, 1)        
        currStart = bspkIdx(ib, 1);
        currEnd   = bspkIdx(ib, 2);

        % IBI: Time between end of previous and start of next
        tEnd   = st(prevEnd);
        tStart = st(currStart);
        currIbi = tStart - tEnd;

        if currIbi < minIbi
            % Merge
            prevEnd = currEnd;
        else
            % Commit & Start New
            mrgIdx = [mrgIdx; prevStart, prevEnd]; 
            prevStart = currStart;
            prevEnd   = currEnd;
        end
    end
    mrgIdx = [mrgIdx; prevStart, prevEnd];

    % ---------------------------------------------------------------------
    % Filter (Duration & nSpikes)
    % ---------------------------------------------------------------------

    finIdx = [];
    for ib = 1:size(mrgIdx, 1)
        currStart = mrgIdx(ib, 1);
        currEnd = mrgIdx(ib, 2);

        n = currEnd - currStart + 1;
        d = st(currEnd) - st(currStart);

        if n >= minSpks && d >= minDur
            finIdx = [finIdx; currStart, currEnd]; 
        end
    end

    if isempty(finIdx)
        continue;
    end

    % ---------------------------------------------------------------------
    % Collect Statistics
    % ---------------------------------------------------------------------

    nb = size(finIdx, 1);
    bTimes = zeros(nb, 2);
    bNspks = zeros(nb, 1);
    bDur   = zeros(nb, 1);
    bFreq  = zeros(nb, 1);
    bIbi   = nan(nb, 1);
    bSt    = [];

    for ib = 1:nb
        currStart = finIdx(ib, 1);
        currEnd = finIdx(ib, 2);

        bTimes(ib, :) = [st(currStart), st(currEnd)];
        bNspks(ib)    = currEnd - currStart + 1;
        bDur(ib)      = st(currEnd) - st(currStart);
        bFreq(ib)     = bNspks(ib) / bDur(ib);

        bChunk = st(currStart:currEnd);
        bSt = [bSt; bChunk]; 
    end

    % IBI Calculation
    if nb > 1
        bIbi(2:end) = bTimes(2:end, 1) - bTimes(1:end-1, 2);
    end

    brst.times{iUnit}    = bTimes;
    brst.nBspk{iUnit}    = bNspks;
    brst.dur{iUnit}      = bDur;
    brst.freq{iUnit}     = bFreq;
    brst.ibi{iUnit}      = bIbi;
    brst.spktimes{iUnit} = bSt;

end


%% ========================================================================
%  PLOT & SAVE
%  ========================================================================

if flgPlot

    [hFig, hAx] = plot_axSize('szOnly', false, 'flgFullscreen', true, ...
        'flgPos', true);
    winPlot = [10, 20];
    lnH = 0.8;
    lnW = 0.8;

    % Plot all spikes
    clr = [0.6 0.6 0.6];
    [hAx, ~] = plot_raster(spktimes, 'PlotType', 'vertline', ...
        'lineHeight', lnH, ...
        'lineWidth', lnW, ...
        'hAx', hAx, ...
        'clr', clr, ...
        'xLim', winPlot);
    hold(hAx, 'on');

    % Plot burst spikes in color (Red)
    clr = [1 0 0];
    plot_raster(brst.spktimes, 'PlotType', 'vertline', ...
        'lineHeight', 0.8, ...
        'lineWidth', lnW, ...
        'hAx', hAx, ...
        'clr', clr, ...
        'xLim', winPlot);

    title(hAx, 'Max Interval Burst Detection');
    xlabel('Time (s)')
    ylabel('Unit No.')

end

if flgSave
    save(saveFile, 'brst');
end

end     % EOF


%% ========================================================================
%  NOTE: REFERENCES
%  ========================================================================
%  - Cotterill, E., et al. (2016). A comparison of computational methods
%    for detecting bursts in neuronal spike trains and their application to
%    human stem cell-derived neuronal networks. J. Neurophysiol.
%  - NeuroExplorer Manual, section 2.24.
%  ========================================================================

%% ========================================================================
%  NOTE: ALGORITHM
%  ========================================================================
%  The Max Interval algorithm identifies bursts based on Inter-Spike
%  Intervals (ISIs) using a four-step process. This implementation follows
%  the standard described by Cotterill et al. (2016) and the NeuroExplorer
%  manual.
%
%  CORE DETECTION: Find "cores" of potential bursts. A core must begin
%  with an ISI strictly less than 'maxISI_start'.
%
%  EXTENSION: Extend the core to subsequent spikes as long as each ISIs
%  remains less than or equal to 'maxISI_end'. The burst ends when an ISI
%  exceeds this threshold.
%
%  MERGING: If the interval between the end of one burst and the start of
%  the next (the Inter-Burst Interval, IBI) is less than 'minIBI', the two
%  bursts are merged into a single event.
%
%  FILTERING: Discard any candidate bursts that do not meet the minimum
%  duration ('minDur') or minimum spike count ('minSpks') criteria.
%  ========================================================================

%% ========================================================================
%  NOTE: MINSPKS
%  ========================================================================
% The choice of `minSpks = 3` is a structural requirement that ensures
% detected events possess a true internal frequency. A two-spike event
% consists of only a single interval, which provides no information
% about whether a neuron has entered a sustained bursting state or
% simply fired a random doublet.
%
% * Reduces Poisson noise
% * Requires internal structure
% * Defines sustained states
%
% By requiring at least three spikes, we ensure that the event contains
% at least two intervals. This drastically increases the statistical
% probability that the event represents a coordinated physiological
% burst rather than an accidental coincidence of spikes within a
% random (Poisson-distributed) spike train.
%  ========================================================================

%% ========================================================================
%  NOTE: BURST START VS. BURST END
%  ========================================================================
% The Max Interval (MI) method uses two distinct temporal parameters to
% model the life cycle of a burst. The `maxISI_start` acts as a high-
% frequency trigger that identifies the initiation of a burst. Once
% initiated, the `maxISI_end` allows the burst to continue even if the
% firing rate slightly tapers off.
%
% This dual-threshold logic accounts for the biological reality that
% bursts often begin with an intense high-frequency onset and
% gradually decelerate before terminating. Cotterill et al. (2016)
% maintained a fixed offset where `maxISI_end` was 0.130s longer than
% `maxISI_start` to preserve this "tapering" window while testing
% robustness across different temporal scales.
%
% This fixed relationship ensures that the detector remains sensitive
% to the primary frequency of the network without being overly
% sensitive to the exact point of termination. It prevents a single
% slightly longer interval from prematurely splitting a single burst
% into two separate events.
%  ========================================================================

%% ========================================================================
%  NOTE: MINIMUM DURATION
%  ========================================================================
% The minDur parameter acts as a qualitative filter that is conceptually
% distinct from frequency-based triggers like maxISI_start. While minSpks
% and maxISI define the mathematical entry criteria for a burst, minDur
% ensures the event represents a sustained physiological state rather than
% a high-frequency "glitch."
% * Frequency vs. Stability: A burst might meet the minSpks count at an
%   extremely high frequency but last only a few milliseconds. The minDur
%   requirement forces the detector to ignore these unstable transients.
% * Filtering Poisson Noise: It helps differentiate true coordinated bursts
%   from accidental coincidences where multiple spikes happen to occur in
%   very close proximity within a high-rate tonic stream.
% * Scaling with Biology: In hippocampal preparations, true bursts have a
%   characteristic temporal envelope; minDur allows the user to exclude
%   events that are too brief to be biologically relevant to the network's
%   signaling.
%  ========================================================================

%% ========================================================================
%  NOTE: MAXISI_END
%  ========================================================================
% The `maxISI_end` parameter defines the termination criteria for a burst
% and is extracted by focusing on the right-hand decay of the first peak
% (the burst peak) in the log-ISI histogram. In biological networks like
% hippocampal cultures, bursts often exhibit a high-frequency onset
% followed by a gradual deceleration in firing rate before the event
% concludes. To identify this value visually, locate the "shoulder" of the
% first peak where the steep negative slope begins to flatten toward the
% horizontal floor of the valley.
%
% This specific shoulder point represents the threshold where a neuron
% has slowed significantly from its peak intra-burst velocity but has not
% yet transitioned into a tonic firing state. Setting the
% parameter at this transition ensures the detector captures the natural
% tapering "tail" of the burst without prematurely splitting a single
% event into multiple smaller segments.
%  ========================================================================

%% ========================================================================
%  NOTE: MINIBI
%  ========================================================================
% The `minIBI` (Minimum Inter-Burst Interval) parameter determines the
% independence of firing events and is extracted by focusing on the
% left-hand rising edge of the second peak (the tonic peak) in the log-ISI
% histogram. While the valley floor represents the statistical minimum
% between states, the `minIBI` should be placed further to the right where
% the histogram bins begin to rise significantly toward the tonic peak.
%
% This "rising edge" visually identifies the shortest intervals that are
% statistically likely to belong to the background background activity
% rather than the internal structure of a burst. selecting this value
% ensureS that two bursts occurring in rapid succession are correctly
% identified as separate physiological events instead of being erroneously
% merged into a single long, reverberating burst.
%  ========================================================================

%% ========================================================================
%  NOTE: BRST_MEA
%  ========================================================================
%  This function can be configured to replicate the behavior of standard
%  simple threshold burst detectors (eg, brst_mea.m) by forcing the
%  two ISI thresholds to be identical and disabling the secondary filters.
%
%  To achieve this equivalence:
%    - Set 'maxISI_start' and 'maxISI_end' to the same value (e.g., 0.02).
%    - Set 'minIBI' to 0 (disables burst merging).
%    - Set 'minDur' to 0 (disables duration filtering).
%
%  Example:
%    brst = brst_maxInt(spktimes, 'maxISI_start', 0.02, 'maxISI_end', 0.02, ...
%           'minIBI', 0, 'minDur', 0, 'minSpks', 2);
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
