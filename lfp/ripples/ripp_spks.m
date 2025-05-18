function ripp = ripp_spks(ripp, varargin)
% RIPP_SPKS Analyzes single-unit (SU) and multi-unit (MU) spiking activity 
% relative to hippocampal ripple events.
%
% SUMMARY:
% This function takes ripple event data and neuronal spiking data to 
% characterize how neuronal firing rates change during ripples compared to
% baseline periods (control events). It performs the following main steps:
%   1. Loads necessary data (ripple structure, session info, spike times).
%   2. Generates control events: periods outside ripples, matched for 
%      duration and optionally limited to a specific vigilance state.
%   3. Calculates SU firing rates during each ripple and control event.
%   4. Performs statistical comparison (Wilcoxon sign-rank test) between 
%      ripple and control firing rates for each SU, optionally limited to 
%      a specified vigilance state. Calculates spike gain (z-scored rate change).
%   5. Generates Peri-Event Time Histograms (PETHs) for MU and SU activity 
%      aligned to ripple peaks and control event centers (outputting spike counts).
%   6. Stores results in the input 'ripp' structure and optionally saves it.
%   7. Optionally plots summary graphics.
%
% INPUT:
%   ripp                Structure containing ripple detection results from
%                       getRipples.m. Must include fields: 
%                         times: [N x 2] ripple start/end times.
%                         peakTime: [N x 1] ripple peak times.
%                         dur: [N x 1] ripple durations.
%                         info.fs: LFP sampling frequency.
%                         info.recWin: [1 x 2] recording window used.
%                       If 'limState' is used, requires:
%                         states.idx: Logical matrix/cell identifying ripples within states.
%   basepath            (Optional) Path to recording session directory {pwd}.
%   flgGraphics         (Optional) Logical flag to plot results {true}.
%   flgSaveVar          (Optional) Logical flag to save/update the ripp
%                       structure in a .mat file {true}.
%   limState            (Optional) Numeric index (e.g., 1-5) specifying a
%                       vigilance state (from sleep_states.mat) to restrict
%                       the *statistical analysis* (gain, p-value) to.
%                       {Default: [], analyze all ripple events for stats}.
%
% OUTPUT:
%   ripp                Updated structure containing ripple analysis results.
%                       See spks_initialize for detailed field descriptions.
%
% DEPENDENCIES:
%   basepaths2vars (custom)
%   SubtractIntervals (FMAToolbox)
%   Restrict (FMAToolbox)
%   plot_rippleSpks (custom) 
%   Requires FMAToolbox in path for SubtractIntervals, Restrict.
%
% HISTORY:
% 06 Aug 24 LH - Major refactor: Changed control event selection, state-limiting application, naming conventions.
% 05 Aug 24 LH - Refactored from rippleSpks.m, integrated bz_getRipSpikes logic,
%                added state-matching, spike gain calculation, and stats.
% 
%%% consider adding phase locking
% https://www.sciencedirect.com/science/article/pii/S2352289521000357#sec2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING & VALIDATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section defines optional input arguments and parses them using 
% inputParser. It sets default values and validates input types.
p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'flgGraphics', true, @islogical);
addOptional(p, 'flgSaveVar', true, @islogical);
addOptional(p, 'limState', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));

parse(p, varargin{:});
basepath        = p.Results.basepath;
flgGraphics     = p.Results.flgGraphics;
flgSaveVar      = p.Results.flgSaveVar;
limState        = p.Results.limState;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARATIONS & DATA LOADING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section sets up file paths, checks for essential fields in the input 
% 'ripp' structure, validates state information if requested, loads 
% necessary session variables (session info, spike data, sleep states) using 
% getSessionVars, and defines key parameters for the analysis.

% files
cd(basepath);
[~, basename] = fileparts(basepath);
rippfile = fullfile(basepath, [basename, '.ripp.mat']);

% Check if essential ripple info exists
if ~isfield(ripp, 'times') || ~isfield(ripp, 'peakTime') || ~isfield(ripp, 'info')
     error('Input ''ripp'' structure is missing required fields.');
end
% Check for state info if limState is requested
if ~isempty(limState) && size(ripp.states.idx, 2) < limState 
     warning('ripp.states.idx{%d} not found.', limState);
     limState = []; % Reset limState if state data isn't available/valid
end

% load session vars
vars = ["session"; "spikes"; "spktimes"; "sleep_states"]; 
v = basepaths2vars('basepaths', {basepath}, 'vars', vars);

% Check if required variables loaded
if ~isfield(v, 'session') || (~isfield(v, 'spikes') && ~isfield(v, 'spktimes'))
    error('Failed to load required session variables (session, spikes/spktimes).');
end

% Check if single units exist
if isfield(v, 'spikes') && ~isempty(v.spikes) && isfield(v.spikes, 'times')
    nunits = length(v.spikes.times);
else
    nunits = 0;
end

% params
fsSpk = v.session.extracellular.sr;
fsLfp = ripp.info.fs;
rippTimes = ripp.times;             % Ripple start/end times 
rippPeakTime = ripp.peakTime;       % Ripple peak times 
rippDur = ripp.dur / 1000;          % Ripple durations (s)
nRipples = size(rippTimes, 1);

% Peri-event map parameters. Use ripple mapDur if already defined,
% otherwise default
if isfield(ripp.maps, 'mapDur')
    mapDur = ripp.maps.mapDur;
else
    mapDur = [-0.075 0.075];        
end
nBinsMap = floor(fsLfp * diff(mapDur) / 2) * 2 + 1;     % Ensure odd number of bins

% Analysis window
recWin = ripp.info.recWin;
if isfield(v.session.general, 'duration') && recWin(2) > floor(v.session.general.duration)
    recWin(2) = floor(v.session.general.duration);
end

% Initialize output structure
spks = spks_initialize();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE CONTROL (NON-RIPPLE) EVENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section generates control time intervals used as a baseline for 
% comparison against ripple events. Control events are matched to ripple 
% durations and are selected from periods outside of ripples (plus a buffer).
% If 'limState' is specified, control intervals are restricted to that state.
% The current implementation sequentially places control events in the first
% available valid interval, ensuring temporal separation (padding) and 
% optionally proximity (maxDist) to the corresponding ripple.

% Define intervals to exclude: ripple events +/- map window buffer
% Use mapDur to define exclusion zone around ripples
buffer = max(rippDur); % Use the longest ripple as a buffer
intervalsExclude = rippTimes + [-buffer, buffer];

% Define intervals to include: only within bouts of a specific state. This
% means all control events will be limited to one state (e.g., NREM)
if limState
    intervalsInclude = v.ss.bouts.times{limState};
else
    intervalsInclude = recWin;
end

% Calculate non-event intervals within the entire recording window
intervalsCtrl = SubtractIntervals(intervalsInclude, intervalsExclude);

% Remove intervals shorter than the longest ripple duration 
durIdx = diff(intervalsCtrl, 1, 2) < max(rippDur);
intervalsCtrl(durIdx, :) = [];

% Check if any valid control intervals exist
if isempty(intervalsCtrl)
    error('No valid non-event intervals found.');
end

% Initialize control structures and parameters
ctrlTimes = nan(nRipples, 2);
padding = 0.1;

% Lazy fix to balance between proximity of ripples and control events, and
% availability of intervals for control. 
maxDist = 4 * 60 * 60;     

for iEvent = 1:nRipples

    % Place control event at the start of the chosen interval
    ctrlStart = intervalsCtrl(1, 1);
    ctrlEnd = ctrlStart + rippDur(iEvent);
    ctrlTimes(iEvent, :) = [ctrlStart, ctrlEnd];
    
    % subtract event from available intervals 
    rmEvent = [ctrlTimes(iEvent, 1) - padding, ctrlTimes(iEvent, 2) + padding];
    intervalsCtrl = SubtractIntervals(intervalsCtrl, rmEvent);
    
    % and make sure next event starts after previous one
    intervalsCtrl(intervalsCtrl(:, 2) < ctrlEnd, :) = [];

    % make sure updated intervals are close to ripple positions 
    rmIdx = (rippPeakTime(iEvent) - intervalsCtrl(:, 1)) > maxDist;
    intervalsCtrl(rmIdx, :) = [];
end

% Calculate control peak positions (e.g., midpoint) and ensure duration
% fits with ripples
ctrlPeakTime = ctrlTimes(:, 1) + diff(ctrlTimes, 1, 2) / 2;
ctrlDur = diff(ctrlTimes, 1, 2); % Recalculate duration based on actual placed events

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTEGRATED SPIKE-IN-EVENT CALCULATION & STATISTICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section calculates the firing rate of each single unit during every 
% ripple event and every control event using times2rate. It then performs 
% statistical comparisons (Wilcoxon signed-rank test) between the distributions 
% of ripple rates and control rates for each unit. If 'limState' is specified, 
% only events within that state are included in the statistical tests and 
% spike gain calculation. 

% Initialize SU results fields
spks.su.pVal = nan(nunits, 1);
spks.su.sigMod = false(nunits, 1);
spks.su.frGain = nan(nunits, 1);
spks.su.frModulation = nan(nunits, 1);
spks.su.frPrct = nan(nunits, 1);

% Calculate SU firing rates during all ripple and control events for all units at once
if nunits > 0
    spks.su.rippRates = times2rate(v.spikes.times, 'winCalc', rippTimes, 'binsize', Inf);
    spks.su.ctrlRates = times2rate(v.spikes.times, 'winCalc', ctrlTimes, 'binsize', Inf);
    
    % Calculate peak rates using a narrow window around ripple and control peaks
    peakWin = [-0.003 0.003]; % +/- 3ms window around peaks
    rippPeakTimes = [rippPeakTime + peakWin(1), rippPeakTime + peakWin(2)];
    spks.su.peakRates = times2rate(v.spikes.times, 'winCalc', rippPeakTimes, 'binsize', Inf);
end

% Determine state index for filtering statistics
if isempty(limState)
    idxState = true(nRipples, 1);
else
    idxState = ripp.states.idx(:, limState);
end

% Statistical Comparison & Gain Calculation (State-Filtered)
for iunit = 1:nunits
    
    % Filter rates based on state index
    rippUnit = spks.su.rippRates(iunit, idxState);
    ctrlUnit = spks.su.ctrlRates(iunit, idxState);

    % Perform statistical test to see if distribution of FRs in ripples
    % is different than in control events. The use of Wilcoxon
    % (dependent samples) is becuase each ripp and ctrl event are
    % matched by vigalance state and event duration.
    % Perform test only if there are non-NaN values in both vectors
    if any(~isnan(rippUnit)) && any(~isnan(ctrlUnit))
        [p, h] = signrank(rippUnit, ctrlUnit);
        spks.su.pVal(iunit) = p;
        spks.su.sigMod(iunit) = h;
    else
        spks.su.pVal(iunit) = NaN;
        spks.su.sigMod(iunit) = false;
    end

    % calculate firing rate increase parameters
    meanCtrlRate = mean(ctrlUnit, 'omitnan');
    sdCtrlRate = std(ctrlUnit, 0, 'omitnan');
    meanRippRate = mean(rippUnit, 'omitnan');
    diffRate = meanRippRate - meanCtrlRate;
    sumRate = meanRippRate + meanCtrlRate;

    spks.su.frGain(iunit) = diffRate / sdCtrlRate;         % X-score
    spks.su.frModulation(iunit) = diffRate / sumRate;      % Modulation index
    spks.su.frPrct(iunit) = diffRate / meanCtrlRate * 100; % Percent change

end % Unit loop 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE PERI-EVENT SPIKE COUNT MAPS (MU and SU)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section generates Peri-Event Time Histograms (PETHs) by aligning 
% multi-unit (MU) and single-unit (SU) spike times to the peak of each 
% ripple event and the center of each control event. The helper function 
% 'sync_spksMap' is used to compute spike counts within time bins 
% relative to the alignment points across all events.

% Single-Unit Maps
spks.su.rippMap = nan(nunits, nRipples, nBinsMap);
spks.su.ctrlMap = nan(nunits, nRipples, nBinsMap);

for iunit = 1:nunits
    unitSpkTimes = v.spikes.times{iunit};
    unitSpkTimes = Restrict(unitSpkTimes, recWin); % Ensure SU spikes are within analysis window

    spks.su.rippMap(iunit, :, :) = sync_spksMap(unitSpkTimes, rippPeakTime, mapDur, nBinsMap);
    spks.su.ctrlMap(iunit, :, :) = sync_spksMap(unitSpkTimes, ctrlPeakTime, mapDur, nBinsMap);
end

% Multi-Unit Maps
spks.mu = []; % Indicate no MU data
if isfield(v, 'spktimes') && ~isempty(v.spktimes)
    muSpkTimes = sort(vertcat(v.spktimes{:})) / fsSpk;      % Combine all tetrode spikes
    muSpkTimes = Restrict(muSpkTimes, recWin);              % Ensure MU spikes are within analysis window

    spks.mu.rippMap = sync_spksMap(muSpkTimes, rippPeakTime, mapDur, nBinsMap);
    spks.mu.ctrlMap = sync_spksMap(muSpkTimes, ctrlPeakTime, mapDur, nBinsMap);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINALIZE OUTPUT, SAVE & GRAPHICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section stores the generated control event times and peak positions 
% in the output structure. If requested ('saveVar' is true), the updated 
% 'ripp' structure is saved to a .mat file. If requested ('graphics' 
% is true), summary plots are generated using 'plot_rippleSpks'.

% store and organize output
spks.info.ctrlPeakTime = ctrlPeakTime; 
spks.info.ctrlTimes = ctrlTimes;   
spks.info.limState = limState;
spks.info.mapDur = mapDur;
spks.info.nBinsMap = nBinsMap;
ripp.spks = spks;

% save
if flgSaveVar
    save(rippfile, 'ripp', '-v7.3');
end

% graphics
if flgGraphics 
    ripp_plotSpks(ripp, 'basepath', basepath, 'flgSaveFig', true);
end

end % End of main function ripp_spks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% sync_spksMap (Non-overlapping windows)
% -------------------------------------------------------------------------
function syncMap = sync_spksMap(spikeTimes, eventTimes, mapDur, nBinsMap)
% Creates a synchronized map (PETH) of spike counts around event times
% using a single histcounts call, assuming non-overlapping event windows.
%
% INPUTS:
%   spikeTimes   Vector of spike timestamps.
%   eventTimes   Vector of event timestamps to align to.
%   mapDur       2-element vector specifying window duration [start, end]
%                relative to eventTime (e.g., [-0.075, 0.075]).
%   nBinsMap     Number of bins for the output map.
%
% OUTPUT:
%   syncMap      Matrix [nEvents x nBinsMap] of spike counts per bin.

    nEvents = length(eventTimes);

    % Initialize output map with zeros
    syncMap = zeros(nEvents, nBinsMap);

    if isempty(spikeTimes) || isempty(eventTimes)
        return; % Return zero map if no spikes or events
    end

    % Ensure eventTimes is a column vector
    eventTimes = eventTimes(:);

    % Create relative bin edges
    relativeEdges = linspace(mapDur(1), mapDur(2), nBinsMap + 1);

    % Create absolute bin edges for all events.
    absoluteEdges = eventTimes + relativeEdges; % nEvents x (nBinsMap+1) matrix

    % Pre-filter spikes to the overall time range covering all event windows for efficiency
    spikeEvents = Restrict(spikeTimes, [absoluteEdges(:, 1), absoluteEdges(:, end);]);
    
    % If no spikes remain after restriction, return the zero map
    if isempty(spikeEvents)
        return;
    end
    
    % counts spikes for each event
    for iEvent = 1 : nEvents
        syncMap(iEvent, :) = histcounts(spikeEvents, absoluteEdges(iEvent, :));
    end

end % End of sync_spksMap

% -------------------------------------------------------------------------
% spks_initialize
% -------------------------------------------------------------------------
function spks = spks_initialize()
% Creates the initial spks structure with all fields set to default/empty.

spks = struct();

% Info substruct (Metadata)
spks.info = struct();
spks.info.limState = [];       % Limiting state index used for stats
spks.info.mapDur = [];         % Window duration used for PETH maps [start, end] in sec
spks.info.nBinsMap = [];       % Number of bins used for PETH maps
spks.info.ctrlPeakTime = [];   % Timestamps used as control event centers
spks.info.ctrlTimes = [];      % Start/end times of generated control events

% Multi-unit (MU) substruct
spks.mu = struct();
spks.mu.rippMap = [];     % [N_ripples x nBinsMap] PETH map of MU spike counts around ripple peaks
spks.mu.ctrlMap = [];     % [N_controls x nBinsMap] PETH map of MU spike counts around control event centers

% Single-unit (SU) substruct
spks.su = struct();
spks.su.peakRates = [];   % [nUnits x N_ripples] firing rates during +/-3ms around ripple peaks
spks.su.rippRates = [];   % [nUnits x N_ripples] firing rates during ripple events
spks.su.ctrlRates = [];   % [nUnits x N_controls] firing rates during control events
spks.su.pVal = [];        % [nUnits x 1] P-value from Wilcoxon sign-rank test
spks.su.sigMod = [];      % [nUnits x 1] Logical flag for significant modulation
spks.su.frGain = [];      % [nUnits x 1] Mean z-scored firing rate during ripples
spks.su.frModulation = [];% [nUnits x 1] Modulation index ((Ripp-Ctrl)/(Ripp+Ctrl))
spks.su.frPrct = [];      % [nUnits x 1] Percent change relative to control rate
spks.su.rippMap = [];     % [nUnits x N_ripples x nBinsMap] PETH map of SU spike counts
spks.su.ctrlMap = [];     % [nUnits x N_controls x nBinsMap] PETH map of SU spike counts

end % End of spks_initialize

% EOF