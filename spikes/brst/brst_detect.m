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
%                     'verbose'   : (log) Print progress to command window {false}
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
addParameter(p, 'verbose', false, @islogical);

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
verbose   = p.Results.verbose;
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

% Initialize parallel pool if needed
if isempty(gcp('nocreate'))
    parpool('local', 6);
end
dq = parallel.pool.DataQueue;
if verbose
    afterEach(dq, @(msg) fprintf('%s\n', msg));
end

% Temporary cell arrays for parfor slicing
bTimes    = cell(nUnits, 1);
bNspks    = cell(nUnits, 1);
bDur      = cell(nUnits, 1);
bFreq     = cell(nUnits, 1);
bIbi      = cell(nUnits, 1);
bSpktimes = cell(nUnits, 1);

parfor iUnit = 1 : nUnits

    st = spktimes{iUnit};

    % Ensure column vector 
    st = st(:);
    nSpk = length(st);
    isi = diff(st);
    
    % Identify intervals smaller than isiEnd
    isEnd = isi <= isiEnd;

    % Skip if not enough spikes
    if isempty(st) || length(st) < minSpks || ~any(isEnd)
        continue;
    end

    % ---------------------------------------------------------------------
    % Burst Detection
    % ---------------------------------------------------------------------

    % Find connected components (potential bursts)
    % Pad with 0 to detect edges
    edges = diff([0; isEnd; 0]);

    % Start indices (in isi vector)
    bStartIsi = find(edges == 1);
    % End indices (in isi vector)
    bEndIsi   = find(edges == -1) - 1;

    nPot = length(bStartIsi);
    validStarts = false(nPot, 1);

    % Refine start times: Must start with isi <= isiStart
    % Trim leading ISIs in the group that are > isiStart
    for iPot = 1:nPot
        idxStart = bStartIsi(iPot);
        idxEnd   = bEndIsi(iPot);

        % Find first ISI in this group that meets Start criteria
        % Indices relative to the group
        chunk   = isi(idxStart:idxEnd);
        validIs = find(chunk <= isiStart, 1, 'first');

        if ~isempty(validIs)
            % Update start index to the first valid start ISI
            bStartIsi(iPot) = idxStart + validIs - 1;
            validStarts(iPot) = true;
        end
    end

    bStartIsi = bStartIsi(validStarts);
    bEndIsi   = bEndIsi(validStarts);

    if isempty(bStartIsi)
        continue;
    end

    % Convert ISI indices to Spike indices
    % If ISI k is the start, then Spike k is the first spike of the burst.
    % If ISI m is the end, then Spike m+1 is the last spike of the burst.
    spkStart = bStartIsi;
    spkEnd   = bEndIsi + 1;

    % ---------------------------------------------------------------------
    % Merge Bursts
    % ---------------------------------------------------------------------


    % Calculate IBIs: Time diff between Start of Next and End of Curr
    tEnd   = st(spkEnd(1:end-1));
    tStart = st(spkStart(2:end));
    currIbi = tStart - tEnd;

    shouldMerge = currIbi < minIbi;

    if any(shouldMerge)
        % Identify groups of bursts to merge

        % Keep a start if it is the first one (idx 1) OR the previous merge was false.
        keepStart = [true; ~shouldMerge];

        % Keep an end if it is the last one (idx end) OR the next merge is false.
        keepEnd   = [~shouldMerge; true];

        spkStart = spkStart(keepStart);
        spkEnd   = spkEnd(keepEnd);
    end

    % ---------------------------------------------------------------------
    % Filter (Duration & nSpikes)
    % ---------------------------------------------------------------------

    dur = st(spkEnd) - st(spkStart);
    n   = spkEnd - spkStart + 1;

    isVld = (n >= minSpks) & (dur >= minDur);

    spkStart = spkStart(isVld);
    spkEnd   = spkEnd(isVld);
    dur      = dur(isVld);
    n        = n(isVld);

    if isempty(spkStart)
        continue;
    end

    % ---------------------------------------------------------------------
    % Collect Statistics
    % ---------------------------------------------------------------------

    nb = length(spkStart);

    % Frequency
    freq = n ./ dur;

    % Pre-burst IBI
    % By definition, the first burst has NaN IBI if we look at previous burst
    ibiVal = nan(nb, 1);
    if nb > 1
        tS = st(spkStart(2:end));
        tE = st(spkEnd(1:end-1));
        ibiVal(2:end) = tS - tE;
    end

    % Times
    timesVal = [st(spkStart), st(spkEnd)];

    % Collect all spikes in bursts (Vertical concatenation of chunks)
    bSpks = cell(nb, 1);
    for ib = 1:nb
        bSpks{ib} = st(spkStart(ib):spkEnd(ib));
    end
    bSpks = vertcat(bSpks{:});

    % Store in temps
    bTimes{iUnit}    = timesVal;
    bNspks{iUnit}    = n;
    bDur{iUnit}      = dur;
    bFreq{iUnit}     = freq;
    bIbi{iUnit}      = ibiVal;
    bSpktimes{iUnit} = bSpks;

    % Periodic update
    if verbose && mod(iUnit, 10) == 0
        send(dq, sprintf('Processed Unit %d / %d', iUnit, nUnits));
    end
end

% Assign back to structure
brst.times    = bTimes;
brst.nBspk    = bNspks;
brst.dur      = bDur;
brst.freq     = bFreq;
brst.ibi      = bIbi;
brst.spktimes = bSpktimes;


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

