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
%                         ripp.times: [N x 2] ripple start/end times.
%                         ripp.peakTime: [N x 1] ripple peak times.
%                         ripp.dur: [N x 1] ripple durations.
%                         ripp.info.fs: LFP sampling frequency.
%                         ripp.info.recWin: [1 x 2] recording window used.
%                       If 'limState' is used, requires:
%                         ripp.states.idx: Logical matrix/cell identifying ripples within states.
%   basepath            (Optional) Path to recording session directory {pwd}.
%   flgGraphics         (Optional) Logical flag to plot results {true}.
%   flgSaveVar          (Optional) Logical flag to save/update the ripp
%                       structure in a .ripp.mat file {true}.
%   limState            (Optional) Numeric index (e.g., 1-5) specifying a
%                       vigilance state (from sleep_states.mat) to restrict
%                       the *statistical analysis* (gain, p-value) to.
%                       {Default: [], analyze all ripple events for stats}.
%
% OUTPUT:
%   ripp                Updated structure containing ripple analysis results,
%                       with the following fields added or updated under
%                       ripp.spks:
%                       ripp.spks.info.limState: Limiting state index used for stats.
%                       ripp.spks.info.mapDur: Window duration used for PETH maps [start, end] in sec.
%                       ripp.spks.info.nBinsMap: Number of bins used for PETH maps.
%                       ripp.spks.rippIdxState: Logical index of ripples included in state-specific statistical analysis (if limState used).
%                       ripp.spks.ctrlPeakTime: Timestamps used as control event centers.
%                       ripp.spks.ctrlTimes: Start/end times of generated control events.
%                       ripp.spks.mu.rippMap: [N_ripples x nBinsMap] PETH map of MU spike counts around ripple peaks.
%                       ripp.spks.mu.ctrlMap: [N_controls x nBinsMap] PETH map of MU spike counts around control event centers.
%                       ripp.spks.su.rippRates: {nUnits x 1} cell array, each cell [N_ripples x 1] firing rates during ripple events.
%                       ripp.spks.su.ctrlRates: {nUnits x 1} cell array, each cell [N_controls x 1] firing rates during control events.
%                       ripp.spks.su.pVal: [nUnits x 1] P-value from Wilcoxon sign-rank test comparing ripple vs control rates per unit (*state-filtered*).
%                       ripp.spks.su.sigMod: [nUnits x 1] Logical flag indicating significant modulation per unit (pVal < 0.05, *state-filtered*).
%                       ripp.spks.su.spikeGain: [nUnits x 1] Mean z-scored firing rate during ripples relative to control rate distribution per unit (*state-filtered*).
%                       ripp.spks.su.frGain: [nUnits x 1] Mean z-scored firing rate during ripples relative to control rate distribution per unit (*state-filtered*).
%                       ripp.spks.su.frModulation: [nUnits x 1] Modulation index ((Ripp-Ctrl)/(Ripp+Ctrl)) per unit (*state-filtered*).
%                       ripp.spks.su.frPrct: [nUnits x 1] Percent change relative to control rate ((Ripp-Ctrl)/Ctrl*100) per unit (*state-filtered*).
%                       ripp.spks.su.rippMap: [nUnits x N_ripples x nBinsMap] PETH map of SU spike counts around ripple peaks.
%                       ripp.spks.su.ctrlMap: [nUnits x N_controls x nBinsMap] PETH map of SU spike counts around control event centers.
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
% (Add current date/initials if modifying significantly)

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
     warning('ripp.states.idx{%d} not found or empty. Run rippleStates.m first or check state index. Proceeding without state limiting.', limState);
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
rippTimes = ripp.times;           % Ripple start/end times [N x 2]
rippPeakTime = ripp.peakTime;       % Ripple peak times [N x 1]
rippDur = ripp.dur;                 % Ripple durations [N x 1]
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
ripp.spks = [];
ripp.spks.info.limState = limState;
ripp.spks.info.mapDur = mapDur;
ripp.spks.info.nBinsMap = nBinsMap;
ripp.spks.rippIdxState = []; 

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
    error('No valid non-event intervals found for control event calculation within the specified state/window.');
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
% ripple event and every control event using the 'calculate_event_rates' 
% helper. It then performs statistical comparisons (Wilcoxon signed-rank test) 
% between the distributions of ripple rates and control rates for each unit. 
% If 'limState' is specified, only events within that state are included in 
% the statistical tests and spike gain calculation. Spike gain is calculated 
% as the z-scored difference between mean ripple rate and mean control rate.

if nunits > 0
    
    % Initialize SU results fields
    ripp.spks.su.rippRates = cell(nunits, 1);
    ripp.spks.su.ctrlRates = cell(nunits, 1);
    ripp.spks.su.pVal = nan(nunits, 1);
    ripp.spks.su.sigMod = false(nunits, 1);
    ripp.spks.su.spikeGain = nan(nunits, 1);
    ripp.spks.su.frGain = nan(nunits, 1);
    ripp.spks.su.frModulation = nan(nunits, 1);
    ripp.spks.su.frPrct = nan(nunits, 1);

    % Calculate SU firing rates during all ripple and control events
    for iunit = 1:nunits
        unitSpkTimes = v.spikes.times{iunit};
        unitSpkTimes = Restrict(unitSpkTimes, recWin); % Ensure SU spikes are within analysis window

        % Calculate rates using helper function
        ripp.spks.su.rippRates{iunit} = calculate_event_rates(unitSpkTimes, rippTimes);
        ripp.spks.su.ctrlRates{iunit} = calculate_event_rates(unitSpkTimes, ctrlTimes);

    end % End unit loop for rate calculation

    % Determine state index for filtering statistics
    if isempty(limState)
        idxState = true(nRipples, 1);
    else
        idxState = ripp.states.idx(:, limState); 
    end

    % Statistical Comparison & Gain Calculation (State-Filtered)
    for iunit = 1:nunits
        % Filter rates based on state index
        rippRatesUnit = ripp.spks.su.rippRates{iunit}(idxState);
        ctrlRatesUnit = ripp.spks.su.ctrlRates{iunit}(idxState);
        
        % Perform statistical test to see if distribution of FRs in ripples
        % is different than in control events. The use of Wilcoxon
        % (dependent samples) is becuase each ripp and ctrl event are
        % matched by vigalance state and event duration. 
        % Perform test only if there are non-NaN values in both vectors
        if any(~isnan(rippRatesUnit)) && any(~isnan(ctrlRatesUnit))
            [p, h] = signrank(rippRatesUnit, ctrlRatesUnit);
            ripp.spks.su.pVal(iunit) = p;
            ripp.spks.su.sigMod(iunit) = h;
        else
            ripp.spks.su.pVal(iunit) = NaN;
            ripp.spks.su.sigMod(iunit) = false;
        end

        % calculate firing rate increase parameters
        meanCtrlRate = mean(ctrlRatesUnit, 'omitnan');
        sdCtrlRate = std(ctrlRatesUnit, 0, 'omitnan');
        meanRippRate = mean(rippRatesUnit, 'omitnan');
        diffRate = meanRippRate - meanCtrlRate;
        sumRate = meanRippRate + meanCtrlRate;

        ripp.spks.su.frGain(iunit) = diffRate / sdCtrlRate;         % X-score
        ripp.spks.su.frModulation(iunit) = diffRate / sumRate;      % Modulation index
        ripp.spks.su.frPrct(iunit) = diffRate / meanCtrlRate * 100; % Percent change

    end % End unit loop for stats
end % End if nunits > 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE PERI-EVENT SPIKE COUNT MAPS (MU and SU)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section generates Peri-Event Time Histograms (PETHs) by aligning 
% multi-unit (MU) and single-unit (SU) spike times to the peak of each 
% ripple event and the center of each control event. The helper function 
% 'create_sync_map' is used to compute spike counts within time bins 
% relative to the alignment points across all events.

% Multi-Unit Maps
if isfield(v, 'spktimes') && ~isempty(v.spktimes)
    muSpkTimes = sort(vertcat(v.spktimes{:})) / fsSpk;      % Combine all tetrode spikes
    muSpkTimes = Restrict(muSpkTimes, recWin);              % Ensure MU spikes are within analysis window

    ripp.spks.mu.rippMap = create_sync_map(muSpkTimes, rippPeakTime, mapDur, nBinsMap);
    ripp.spks.mu.ctrlMap = create_sync_map(muSpkTimes, ctrlPeakTime, mapDur, nBinsMap);
else
    ripp.spks.mu = []; % Indicate no MU data
end

% Single-Unit Maps
if nunits > 0

    % Initialize
    ripp.spks.su.rippMap = nan(nunits, nRipples, nBinsMap);
    ripp.spks.su.ctrlMap = nan(nunits, nRipples, nBinsMap);

    for iunit = 1:nunits
        unitSpkTimes = v.spikes.times{iunit};
        unitSpkTimes = Restrict(unitSpkTimes, recWin); % Ensure SU spikes are within analysis window

        ripp.spks.su.rippMap(iunit, :, :) = create_sync_map(unitSpkTimes, rippPeakTime, mapDur, nBinsMap);
        ripp.spks.su.ctrlMap(iunit, :, :) = create_sync_map(unitSpkTimes, ctrlPeakTime, mapDur, nBinsMap);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINALIZE OUTPUT, SAVE & GRAPHICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section stores the generated control event times and peak positions 
% in the output structure. If requested ('saveVar' is true), the updated 
% 'ripp' structure is saved to a .ripp.mat file. If requested ('graphics' 
% is true), summary plots are generated using 'plot_rippleSpks'.

% store and organize output
ripp.spks.ctrlPeakTime = ctrlPeakTime; 
ripp.spks.ctrlTimes = ctrlTimes;   

% save
if flgSaveVar
    save(rippfile, 'ripp', '-v7.3');
end

% graphics
if flgGraphics 
    plot_rippleSpks(ripp, 'basepath', basepath, 'saveFig', flgSaveVar); 
end


end % End of main function ripp_spks

% -------------------------------------------------------------------------
% Helper Function: calculate_event_rates
% -------------------------------------------------------------------------
function eventRates = calculate_event_rates(spikeTimes, eventTimes)
% Calculates firing rates within specified event intervals using histcounts.
%
% INPUTS:
%   spikeTimes      Vector of spike timestamps.
%   eventIntervals  [N x 2] matrix of event start and end times.
%
% OUTPUT:
%   eventRates      [N x 1] vector of firing rates (spikes/sec) per event.

    if isempty(spikeTimes) || isempty(eventTimes)
        eventRates = nan(size(eventTimes, 1), 1);
        return;
    end
    
    % Create bin edges for histcounts: [start1, end1, start2, end2, ...]
    edges = reshape([eventTimes(:, 1), eventTimes(:, 2)]', [], 1);

    % Use histcounts
    counts = histcounts(spikeTimes, edges);

    % Counts within intervals are in the odd bins (1, 3, 5, ...)
    intervalCounts = counts(1:2:end)'; % Ensure column vector

    % Calculate rates
    eventDur = diff(eventTimes, [], 2);
    eventRates = intervalCounts ./ eventDur;

end % End of helper function calculate_event_rates

% -------------------------------------------------------------------------
% Helper Function: create_sync_map (Non-overlapping windows)
% -------------------------------------------------------------------------
function syncMap = create_sync_map(spikeTimes, eventTimes, mapDur, nBinsMap)
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

end % End of helper function create_sync_map

% EOF