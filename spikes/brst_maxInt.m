function brst = brst_maxInt(spktimes, varargin)
% BRST_MAXINT Detects bursts using the Max Interval method.
%
%   brst = BRST_MAXINT(SPKTIMES, ...) implements the Max Interval burst
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
%                     'maxISI_start' : (num) Max ISI to start burst {0.013s}
%                     'maxISI_end'   : (num) Max ISI within burst {0.022s}
%                     'minIBI'       : (num) Min Inter-Burst Interval {0.05s}
%                     'minDur'       : (num) Min burst duration {0.015s}
%                     'minSpks'      : (num) Min spikes in burst {3}
%                     'flgPlot'      : (log) Plot raster with bursts {false}
%                     'basepath'     : (char) Recording path {pwd}
%                     'flgSave'      : (log) Save to file {true}
%                     'flgForce'     : (log) Analyze even if exists {false}
%
%   OUTPUTS:
%       brst        - (struct) Burst statistics and event data.
%                     Fields: detect, nspks, brstDur, freq, bspks, ibi, all.
%                     (All vector fields are nUnits x 1 column vectors)
%

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'maxISI_start', 0.013, @isnumeric);
addParameter(p, 'maxISI_end', 0.02, @isnumeric);
addParameter(p, 'minIBI', 0.05, @isnumeric);
addParameter(p, 'minDur', 0.015, @isnumeric);
addParameter(p, 'minSpks', 3, @isnumeric);
addParameter(p, 'flgPlot', false, @islogical);
addParameter(p, 'flgSave', true, @islogical);
addParameter(p, 'flgForce', false, @islogical);

parse(p, varargin{:});
basepath     = p.Results.basepath;
maxISI_start = p.Results.maxISI_start;
maxISI_end   = p.Results.maxISI_end;
minIBI       = p.Results.minIBI;
minDur       = p.Results.minDur;
minSpks      = p.Results.minSpks;
flgPlot      = p.Results.flgPlot;
flgSave      = p.Results.flgSave;
flgForce     = p.Results.flgForce;


%% ========================================================================
%  PREPARATIONS
%  ========================================================================

% Load if exists
[~, basename] = fileparts(basepath);
brstfile = fullfile(basepath, [basename, '.brst.mat']);

if exist(brstfile, 'file') && ~flgForce && ~flgPlot
    load(brstfile, 'brst');
    return;
end

nunits = length(spktimes);

% Initialize output structure
% CHANGED: Ensure column vectors (nUnits x 1)
brst.info.runtime   = datetime("now");
brst.info.algorithm = 'MaxInterval';
brst.info.input     = p.Results;

brst.detect         = zeros(nunits, 1);
brst.nspks          = nan(nunits, 1);
brst.brstDur        = nan(nunits, 1);
brst.freq           = nan(nunits, 1);
brst.bspks          = zeros(nunits, 1);
brst.ibi            = nan(nunits, 1);
brst.all            = cell(nunits, 1);

% For plotting
spktimes_bursts = cell(nunits, 1);


%% ========================================================================
%  RUN DETECTION
%  ========================================================================

for iunit = 1 : nunits

    st = spktimes{iunit};
    if isempty(st) || length(st) < minSpks
        continue;
    end

    % Ensure column vector
    st = st(:);
    uISI = diff(st);
    nSpks = length(st);

    % ---------------------------------------------------------------------
    % Phase 1: Core Burst Detection
    % ---------------------------------------------------------------------

    raw_bursts_idx = [];
    inBurst = false;
    currentBurstStartIdx = -1;

    k = 1;
    while k < nSpks
        if ~inBurst
            % Start criteria
            if uISI(k) <= maxISI_start
                inBurst = true;
                currentBurstStartIdx = k;
                k = k + 1;
            else
                k = k + 1;
            end
        else
            % Continuation criteria
            if uISI(k) <= maxISI_end
                k = k + 1;
            else
                % End of burst
                currentBurstEndIdx = k;
                raw_bursts_idx = [raw_bursts_idx; currentBurstStartIdx, currentBurstEndIdx]; %#ok<AGROW>
                inBurst = false;
                k = k + 1;
            end
        end
    end

    % Close last burst if active
    if inBurst
        raw_bursts_idx = [raw_bursts_idx; currentBurstStartIdx, nSpks]; %#ok<AGROW>
    end

    if isempty(raw_bursts_idx)
        continue;
    end

    % ---------------------------------------------------------------------
    % Phase 2: Merge Bursts
    % ---------------------------------------------------------------------

    merged_bursts_idx = [];
    if size(raw_bursts_idx, 1) > 0
        curr_start = raw_bursts_idx(1, 1);
        curr_end   = raw_bursts_idx(1, 2);

        for b = 2:size(raw_bursts_idx, 1)
            next_start = raw_bursts_idx(b, 1);
            next_end   = raw_bursts_idx(b, 2);

            % IBI: Time between end of previous and start of next
            t_end   = st(curr_end);
            t_start = st(next_start);
            ibi     = t_start - t_end;

            if ibi < minIBI
                % Merge
                curr_end = next_end;
            else
                % Commit & Start New
                merged_bursts_idx = [merged_bursts_idx; curr_start, curr_end]; %#ok<AGROW>
                curr_start = next_start;
                curr_end   = next_end;
            end
        end
        merged_bursts_idx = [merged_bursts_idx; curr_start, curr_end]; %#ok<AGROW>
    end

    % ---------------------------------------------------------------------
    % Phase 3: Filter (Duration & nSpikes)
    % ---------------------------------------------------------------------

    final_bursts_idx = [];
    for b = 1:size(merged_bursts_idx, 1)
        s_idx = merged_bursts_idx(b, 1);
        e_idx = merged_bursts_idx(b, 2);

        n_spks = e_idx - s_idx + 1;
        dur    = st(e_idx) - st(s_idx);

        if n_spks >= minSpks && dur >= minDur
            final_bursts_idx = [final_bursts_idx; s_idx, e_idx]; %#ok<AGROW>
        end
    end

    if isempty(final_bursts_idx)
        continue;
    end

    % ---------------------------------------------------------------------
    % Collect Statistics
    % ---------------------------------------------------------------------

    nb = size(final_bursts_idx, 1);

    b_struct.times = zeros(nb, 2);
    b_struct.nspks = zeros(nb, 1);
    b_struct.dur   = zeros(nb, 1);
    b_struct.freq  = zeros(nb, 1);
    b_struct.ibi   = nan(nb, 1);

    % For IBI calc
    fl_times = zeros(nb, 2);

    for b = 1:nb
        s_idx = final_bursts_idx(b, 1);
        e_idx = final_bursts_idx(b, 2);

        b_struct.times(b, :) = [st(s_idx), st(e_idx)];
        b_struct.nspks(b)    = e_idx - s_idx + 1;
        b_struct.dur(b)      = st(e_idx) - st(s_idx);
        b_struct.freq(b)     = b_struct.nspks(b) / b_struct.dur(b);

        fl_times(b, :) = [st(s_idx), st(e_idx)];

        if flgPlot
            spktimes_bursts{iunit} = [spktimes_bursts{iunit}; st(s_idx:e_idx)]; %#ok<AGROW>
        end
    end

    if nb > 1
        b_struct.ibi(2:end) = fl_times(2:end, 1) - fl_times(1:end-1, 2);
    end

    total_spikes_in_bursts = sum(b_struct.nspks);

    brst.detect(iunit)  = nb;
    brst.nspks(iunit)   = mean(b_struct.nspks);
    brst.brstDur(iunit) = mean(b_struct.dur);
    brst.freq(iunit)    = mean(b_struct.freq);
    brst.bspks(iunit)   = total_spikes_in_bursts / length(st);
    brst.ibi(iunit)     = mean(b_struct.ibi, 'omitnan');
    brst.all{iunit}     = b_struct;
end


%% ========================================================================
%  PLOTTING
%  ========================================================================

if flgPlot

    [hFig, hAx] = plot_axSize('szOnly', false, 'flgFullscreen', true, ...
        'flgPos', true);
    winPlot = [10, 20];
    lnHeight = 0.8;
    lnWidth = 0.8;

    % Plot all spikes
    clr = [0.6 0.6 0.6];
    [hAx, hPlt] = plot_raster(spktimes, 'PlotType', 'vertline', ...
        'lineHeight', lnHeight, ...
        'lineWidth', lnWidth, ...
        'hAx', hAx, ...
        'clr', clr, ...
        'xLim', winPlot);
    hold(hAx, 'on');

    % Plot burst spikes in color (Red)
    clr = [1 0 0];
    [hAx, hPlt] = plot_raster(spktimes_bursts, 'PlotType', 'vertline', ...
        'lineHeight', 0.8, ...
        'lineWidth', lnWidth, ...
        'hAx', hAx, ...
        'clr', clr, ...
        'xLim', winPlot);

    title(hAx, 'Max Interval Burst Detection');
    xlabel('Time (s)')
    ylabel('Unit No.')

end


%% ========================================================================
%  SAVE
%  ========================================================================

if flgSave
    save(brstfile, 'brst');
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
%  NOTE: MINIMUM SPIKE COUNT (minSpks)
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
%  NOTE: PARAMETER DIVERSITY
%  ========================================================================
% The default Max Interval (MI) threshold of 170ms (~6 Hz) used in
% Cotterill et al. (2016) was optimized for mouse RGCs and human iPSC
% networks, which exhibit significantly lower burst prevalence and higher
% variability than rodent hippocampal cultures.
%
% In hippocampal preparations, the "natural valley" separating bursts
% from tonic activity typically occurs at a much higher frequency
% (e.g., ~76 Hz or 13ms).
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