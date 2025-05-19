function [spkLfp] = spklfp_calc(varargin)
% SPKLFP_CALC Analyzes spike-LFP coupling characteristics.
%
% SUMMARY:
% This function investigates the relationship between neuronal spiking activity
% and local field potential (LFP) oscillations, primarily within a specific
% frequency band. It calculates several key metrics to characterize this coupling:
%   1. Phase Coupling: Assesses if spikes preferentially occur at a particular
%      phase of the LFP oscillation. Metrics include mean phase, mean resultant
%      length (MRL), Von Mises concentration (kappa), and Rayleigh test p-value.
%      LFP segments for this analysis can be restricted to high-power bouts.
%   2. Power-Phase Rate Maps: Generates a 2D histogram illustrating spike rate
%      as a function of LFP phase and LFP power (amplitude).
%   3. Rate-Magnitude Correlation: Computes the correlation between binned
%      spike rates and the instantaneous LFP magnitude (envelope).
%   4. Population Synchrony: Examines the synchrony of neuronal populations
%      (e.g., regular-spiking vs. fast-spiking units) and its relationship
%      to LFP phase and magnitude.
%   5. Spike-Triggered LFP (Optional): Computes the average LFP waveform
%      aligned to spike times, typically using the filtered LFP.
%
% METHODOLOGY:
% The LFP signal can be provided directly or loaded from disk. If loaded, it
% is filtered for the specified frequency band. The Hilbert transform is then
% applied to the (filtered) LFP to extract its instantaneous phase and
% amplitude (power). Spikes are assigned the LFP phase and power
% corresponding to their time of occurrence. Analysis windows can be further
% refined based on LFP power.
%
% THEORETICAL & PRACTICAL CONSIDERATIONS:
% - Builds upon concepts from bz_GenSpikeLFPCoupling, bz_PhaseModulation,
%   and bz_PowerPhaserateMap.
% - Alternative spike-field coupling tools include Chronux's `coherencypt`
%   (for multitaper cross-spectrograms) and FieldTrip's STA methods, often
%   used for task-based data with repeated trials.
% - Filtering: Employs Butterworth filtering if LFP is loaded. Other methods
%   (e.g., wavelet convolution) exist for different analytical goals.
% - Multi-channel LFP: If multiple LFP channels are specified for loading,
%   they are averaged. Advanced handling (e.g., source localization,
%   channel selection based on signal quality) is not implemented.
% - Spike Bleed-through: LFP recorded on the same channel as spikes may be
%   contaminated by spike waveform energy, potentially impacting coupling
%   measures. This is a common concern in spike-field analysis.
%
% INPUT (Optional Key-Value Pairs):
%   basepath    (char) Full path to the recording folder. {pwd}
%   sig         (numeric vector) Pre-filtered LFP signal for the band of
%               interest. If empty, LFP will be loaded and filtered.
%   spktimes    (cell array) Spike times for each unit (in seconds).
%               Assumed relative to 'sig' start if 'sig' is provided,
%               otherwise absolute times from the recording.
%   fs          (numeric) Sampling frequency of LFP data [Hz]. If 'sig' is
%               empty and 'fs' is not provided or invalid, it's loaded
%               from 'session.mat'. {1250}
%   fRange      (1x2 numeric) Frequency range [low high] for LFP filtering
%               and analysis [Hz]. {[1.5 100]}
%   lfpTimes    (numeric matrix) N x 2 matrix of [start end] times (seconds)
%               defining specific windows for analysis. If 'sig' is empty,
%               these times (if finite) guide LFP loading. Internally, these
%               define analysis intervals, overriding 'powThr' if specified.
%               Assumed relative to 'sig' if 'sig' is provided, otherwise
%               absolute. {[] -> full duration}
%   ch          (numeric vector) LFP channel(s) to load if 'sig' is empty.
%               If multiple, they are averaged. {[] -> defaults to ripple
%               channel from 'session.mat' or channel 1}
%   bit2uv      (numeric) Conversion factor from bits to microvolts for LFP
%               loading. {0.195}
%   flgGraphics (logical) Flag to generate and save summary plots. {true}
%   flgSave     (logical) Flag to save the output 'spkLfp' structure. {true}
%   flgStl      (logical) Flag to compute and include spike-triggered LFP. {false}
%   powThr      (numeric) STD above mean power in band to define high-power
%               LFP bouts for phase coupling analysis. If empty or 'lfpTimes'
%               is provided, no power thresholding is applied. {[]}
%
% OUTPUT:
%   spkLfp      Structure containing spike-LFP coupling analysis results:
%     .info     Metadata (runtime, fRange, fs, lfpTimes, powThr, etc.).
%     .phase    Phase coupling (dist, kappa, theta, mrl, p, spkPhase).
%     .rateMap  Power-phase rate maps (occupancy, counts, rate, bins).
%     .rateMag  Rate-magnitude correlation (r, p).
%     .pop      Population synchrony results (phase coupling, synchrony-magnitude correlation).
%     .stl      Spike-triggered LFP results (if flgStl is true).
%
% DEPENDENCIES:
%   CircularDistribution (FMAT), NormToInt (buzcode),
%   bz_SpktToSpkmat (buzcode), SubtractIntervals (FMAT), binary2bouts (custom),
%   filterLFP (custom), binary_load (custom/buzcode), fastrms (custom/buzcode),
%   InIntervals (buzcode), Restrict (buzcode).
%
% TO DO LIST:
%   - Consider jitter/surrogate methods for statistical significance testing.
%   - Option for sorting units in plots based on coupling strength or other metrics.
%
% HISTORY:
%   Feb 2022 LH - Original version.
%   Aug 2024 (AI Assisted) - Added LFP loading, enhanced documentation,
%                            variable name changes (camelCase), silent mode,
%                            output structure refinement.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING & INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section defines and parses input arguments using MATLAB's inputParser.
% Default values are assigned, and basic parameters for internal calculations
% and file paths are established.
p = inputParser;
addParameter(p, 'basepath', pwd, @ischar)
addParameter(p, 'sig', [], @isnumeric)
addParameter(p, 'spktimes', {}, @iscell)
addParameter(p, 'fs', 1250, @isnumeric)
addParameter(p, 'fRange', [0 0]);
addParameter(p, 'lfpTimes', [], @isnumeric)
addParameter(p, 'ch', [], @isnumeric)
addParameter(p, 'bit2uv', 0.195, @isnumeric);
addParameter(p, 'flgGraphics', true, @islogical)
addParameter(p, 'flgSave', true, @islogical)
addParameter(p, 'flgStl', false, @islogical)
addParameter(p, 'powThr', [], @isnumeric);

parse(p, varargin{:});
basepath        = p.Results.basepath;
sig             = p.Results.sig;
spkTimes        = p.Results.spktimes;
fs              = p.Results.fs;
fRange          = p.Results.fRange;
lfpTimes        = p.Results.lfpTimes;
ch              = p.Results.ch;
bit2uv          = p.Results.bit2uv;
flgGraphics     = p.Results.flgGraphics;
flgSave         = p.Results.flgSave;
flgStl          = p.Results.flgStl;
powThr          = p.Results.powThr;

% Parameters for internal calculations
nBinsPhase = 180;               % Bins for phase distribution
nBinsMap = 20;                  % Bins for power/phase in rate map

% File paths
[~, basename] = fileparts(basepath);
cd(basepath);
sessionFile = fullfile(basepath, [basename, '.session.mat']);
lfpFile = fullfile(basepath, [basename, '.lfp']);
unitsFile = fullfile(basepath, [basename, '.units.mat']);
spksfile = fullfile(basepath, [basename, '.spikes.cellinfo.mat']);
spkLfpFile = fullfile(basepath, sprintf('%s.spklfp_%.0f-%.0fHz.mat',...
    basename, fRange(1), fRange(2)));

% LFP file parameters
sessionData = load(sessionFile, 'session');
if isempty(fs)
    fs = sessionData.session.extracellular.srLfp;
end
nChans = sessionData.session.extracellular.nChannels;
if isempty(ch)
    if isfield(sessionData.session.channelTags, 'Ripple') && ~isempty(sessionData.session.channelTags.Ripple)
        ch = sessionData.session.channelTags.Ripple;
    else
        ch = 1;
    end
end

% load spike times
if isempty(spkTimes)
    load(spksfile, 'spikes')
    spkTimes = spikes.times;
end
nUnits = length(spikes.times);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA LOADING & PREPARATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section manages the LFP signal ('sig').
% If 'sig' is not provided as input, LFP data is loaded from binary files
% for the specified 'ch' (channel(s)) and 'lfpTimes' (analysis windows).
% The loaded LFP is then filtered within the 'fRange'.
% If 'sig' is provided, it's assumed to be already filtered for 'fRange'.
% Note: Spike times are adjusted relative to 'lfpTimes' and 'sigOffset' later,
% in the 'PREPARE SPIKES FOR ANALYSIS' section, after 'sig' is finalized.

sigOffset = 0; % Offset (in seconds) applied to spkTimes and lfpTimes if LFP is loaded from a specific start time.
loadDur = Inf; % Duration of LFP to load (in seconds). Infinite by default, can be restricted by lfpTimes.

if isempty(sig)
    % If 'sig' (LFP signal) is not provided, load it from 'lfpFile'.
    % 'lfpTimes', if provided, determines the start offset and duration for loading.
    if ~isempty(lfpTimes) 
        sigOffset = min(lfpTimes(:));
        loadDur = max(lfpTimes(:)) - sigOffset;
        lfpTimes = lfpTimes - sigOffset;
    end
    
    % Load LFP data, average multi-channel, and bandpass filter.
    sig = binary_load(lfpFile, 'duration', loadDur, 'fs', fs, 'nCh', nChans,...
        'start', sigOffset, 'ch', ch, 'downsample', 1, 'bit2uv', bit2uv);
    sig = mean(sig, 2);
    sig = filterLFP(sig, 'fs', fs, 'type', 'butter', 'dataOnly', true,...
        'order', 3, 'passband', fRange, 'graphics', false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIGNAL CHARACTERIZATION & ANALYSIS WINDOW DEFINITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section characterizes the LFP signal and defines the precise analysis
% windows ('lfpTimes'). It calculates timestamps for the LFP signal, performs
% the Hilbert transform to obtain instantaneous amplitude and phase, and
% normalizes the power. Optionally, 'lfpTimes' can be refined by
% thresholding LFP power, restricting analysis to high-power oscillatory bouts.

tStamps = (0:length(sig)-1)' / fs; % Timestamps relative to the start of 'sig'.

% Hilbert transform for analytical signal.
sigH = hilbert(sig);

% Instantaneous power (Z-scored log-amplitude) and phase [0, 2*pi].
sigPow = NormToInt(log10(abs(sigH)), 'Z', [0 Inf], fs);
sigPhase = mod(angle(sigH), 2 * pi);

% Optionally restrict 'lfpTimes' to high-power LFP bouts if 'powThr' is set
% and 'lfpTimes' were not provided as input.
if isempty(lfpTimes) && ~isempty(powThr) 
    minCycles = 2;                  % Minimum duration of LFP bout in cycles of slowest freq
    minSamples = ceil((fs / fRange(2)) * minCycles);

    sigRms = fastrms(sig, ceil(fs ./ fRange(1)), 1); % Smoothed RMS.
    minRms = mean(sigRms) + std(sigRms) * powThr; % Power threshold.

    % Identify and exclude low-power periods.
    lowPow = binary2bouts('vec', sigRms < minRms, 'minDur', minSamples,...
        'maxDur', Inf, 'interDur', 0);
    lfpTimes = SubtractIntervals(lfpTimes, lowPow / fs);
end

% Validate 'lfpTimes': clip to signal boundaries, remove empty/invalid intervals.
lfpTimes(lfpTimes < 0) = 0;
lfpTimes(lfpTimes > tStamps(end)) = tStamps(end);
lfpTimes = lfpTimes(lfpTimes(:,2) > lfpTimes(:,1), :);

% clear memory
clear sigRms lowPow sigH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARE SPIKES FOR ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section prepares spike data for coupling analysis.
% Spike times ('spkTimes') are adjusted: spikes outside the analysis windows
% defined by 'lfpTimes' (relative to the start of 'sig') are removed,
% and the remaining spike times are made relative to 'sigOffset' if LFP was loaded.
% Then, LFP instantaneous power and phase are interpolated at each adjusted spike time.

% Adjust spike times: select spikes within 'lfpTimes' and make them
% relative to 'sigOffset'. This ensures spike times are aligned with the
% potentially loaded and offset 'sig' LFP data.
spkTimes = cellfun(@(x) x - sigOffset,...
    spkTimes, 'uni', false);
spkTimes = cellfun(@(x) x(InIntervals(x, lfpTimes)),...
    spkTimes, 'uni', false);

% Interpolate LFP power at each adjusted spike time.
% 'nearest' interpolation is used.
spkPow = cellfun(@(x) interp1(tStamps, sigPow, x, 'nearest'),...
    spkTimes, 'uni', false);

% Interpolate LFP phase at each adjusted spike time.
spkPhase = cellfun(@(x) interp1(tStamps, sigPhase, x, 'nearest'),...
    spkTimes, 'uni', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE MODULATION PER CELL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section quantifies the phase preference of each cell's spiking activity
% relative to the LFP oscillations. It calculates key circular statistics
% for spikes occurring within the defined 'lfpTimes' 

phase = struct();
phase.dist      = zeros(nBinsPhase, nUnits);
phase.kappa     = zeros(1, nUnits);     % Von Mises concentration parameter
phase.theta     = zeros(1, nUnits);     % Mean phase
phase.mrl       = zeros(1, nUnits);     % Mean resultant length
phase.pVal      = zeros(1, nUnits);     % Rayleigh test p-value for non-uniformity
phase.bins      = [];                   % Phase bins for distribution

for iUnit = 1 : nUnits
    
    % Require a minimum number of spikes for robust analysis.
    % spkTimes{iUnit} now contains only spikes within lfpTimes, relative to sigOffset.
    if length(spkTimes{iUnit}) < 10
        continue
    end

    % Calculate circular statistics using pre-interpolated phases for the unit's spikes.
    % spkPhase{iUnit} corresponds to spikes already filtered by lfpTimes.
    [phase.dist(:, iUnit), phase.bins, tmpStats] = ...
        CircularDistribution(spkPhase{iUnit}, 'nBins', nBinsPhase);
    phase.kappa(iUnit) = tmpStats.k;
    phase.theta(iUnit) = mod(tmpStats.m, 2*pi); 
    phase.mrl(iUnit) = tmpStats.r;
    phase.pVal(iUnit) = tmpStats.p;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RATE MAP AS A FUNCTION OF POWER AND PHASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section constructs a 2D rate map (or PETH) showing how neuronal firing
% rate varies as a function of LFP power and LFP phase. This analysis uses
% all spikes (already filtered by lfpTimes) and the LFP segment defined by
% 'lfpTimes' for occupancy calculations.

% Initialize structure for rate map results.
rateMap = struct();

% Define LFP power and phase bins for the rate map. Power bins can be
% percentiles of LFP power during 'lfpTimes' (alt = 1) or percentiles based
% on power during spikes (alt = 2)
powAlt = 2;
switch powAlt
    case 1
        lfpIdx = Restrict(1 : length(sigPow), round(lfpTimes * fs));
        lfpPow = sigPow(lfpIdx);
        powEdges = linspace(prctile(lfpPow, 5), prctile(lfpPow, 97), nBinsMap + 1);
    case 2
        spkPowVec = vertcat(spkPow{:});
        powEdges = linspace(prctile(spkPowVec, 5), prctile(spkPowVec, 97), nBinsMap + 1);
end
phaseEdges = linspace(0, 2 * pi, nBinsMap + 1);
rateMap.phaseBins = phaseEdges(1 : end - 1) + 0.5 .* diff(phaseEdges(1 : 2));
rateMap.powBins = powEdges(1 : end - 1) + 0.5 .* diff(powEdges(1 : 2));
powEdges(1) = -Inf; powEdges(end) = Inf; % Capture all power values.

% Calculate occupancy: duration LFP spent in each power-phase bin within 'lfpTimes'.
rateMap.occupancy = ...
    histcounts2(sigPow(lfpIdx), sigPhase(lfpIdx), powEdges, phaseEdges) / fs;

% Count spikes in each power-phase bin for each unit.
% 'spkPow' and 'spkPhase' contain power/phase for spikes already within lfpTimes.
rateMap.counts = cellfun(@(x, y) histcounts2(x, y, powEdges, phaseEdges),...
    spkPow, spkPhase, 'uni', false);

% normalize counts to rate by dividing with occupancy
rateMap.rate = cellfun(@(x) x ./ rateMap.occupancy,...
    rateMap.counts, 'uni', false);

% convert to 3d mat (power x phase x cell)
rateMap.counts = cat(3, rateMap.counts{:});
rateMap.rate = cat(3, rateMap.rate{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RATE - MAGNITUDE CORRELATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section assesses the correlation between neuronal firing rates (binned
% spike counts) and the instantaneous LFP magnitude (absolute power).
% Spike counts are binned, and Spearman correlation with LFP power is calculated.

% Define parameters for binning spike counts.
movWin = 0.02;                  % Moving window [s].
stepSize = 0.005;               % Step size for spike counts [s]

% Create a matrix of binned spike counts.
% spkTimes are already filtered by lfpTimes and relative to sigOffset.
spkcnts = bz_SpktToSpkmat(spkTimes, 'binsize', movWin, 'dt', stepSize);
spkcnts.data = double(spkcnts.data); % Ensure double type.

% Interpolate LFP power/phase to timestamps of binned spike counts.
spkcnts.lfpPow = interp1(tStamps, sigPow, spkcnts.timestamps, 'nearest');
spkcnts.lfpPhase = interp1(tStamps, sigPhase, spkcnts.timestamps, 'nearest');

% Spearman correlation: binned spike counts vs. LFP power.
[rateMag.r, rateMag.pVal] = corr(spkcnts.data, abs(spkcnts.lfpPow),...
    'type', 'spearman', 'rows', 'complete');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POPULATION SYNCHRONY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section analyzes collective activity (synchrony) of neuronal
% subpopulations (e.g., putative pyramidal vs. interneurons) and
% relates this synchrony to LFP phase and magnitude.

pop = struct();
pop.names = ["pPYR"; "pPV"]; 
nPop = length(pop.names);

% Load unit classification and assign to populations.
load(unitsFile, 'units');
pop.unitIdx = zeros(1, nUnits); 
pop.unitIdx(units.clean(1, :)) = 1; % e.g., pPV
pop.unitIdx(units.clean(2, :)) = 2; % e.g., pPYR

% Initialize fields for population phase coupling.
pop.phase.dist      = nan(nBinsPhase, nPop);
pop.phase.kappa     = nan(1, nPop);
pop.phase.theta     = nan(1, nPop);
pop.phase.mrl       = nan(1, nPop);
pop.phase.pVal      = nan(1, nPop);
pop.phase.bins      = phase.bins; % Use same phase bins as single unit.

% Analyze each defined population.
for iPop = 1 : nPop
      
    % Calculate population synchrony (fraction of active cells per bin).
    popIdx = pop.unitIdx == iPop;
    popSyn = sum(spkcnts.data(:, popIdx) > 0, 2) ./...
        sum(popIdx);
  
    % Standardize synchrony (e.g., to its mean).
    popSyn = popSyn ./ mean(popSyn, 'omitnan');
    
    % Correlate population synchrony with LFP magnitude.
    [pop.synMag.r(:, iPop), pop.synMag.pVal(:, iPop)] =...
        corr(popSyn, abs(spkcnts.lfpPow), 'type', 'spearman',...
        'rows', 'complete');

    % Phase coupling of population synchrony.
    popPhase = spkcnts.lfpPhase(popSyn > 0); % LFP phases at synchrony events.

    [pop.phase.dist(:, iPop), pop.phase.bins, tmpStats] =...
        CircularDistribution(popPhase, 'nBins', nBinsPhase);
    pop.phase.kappa(iPop) = tmpStats.k;
    pop.phase.theta(iPop) = mod(tmpStats.m, 2*pi);
    pop.phase.mrl(iPop) = tmpStats.r;
    pop.phase.pVal(iPop) = tmpStats.p;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPIKE TRIGGERED LFP (OPTIONAL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section computes the Spike-Triggered LFP (STL) if 'flgStl' is true.
% It uses the *filtered* LFP ('sig'). For a conventional wideband STA,
% the original broadband LFP should be used.

stl = []; 
if flgStl
    % Call 'spklfp_stl' with current filtered 'sig' and adjusted 'spkTimes'.
    % 'lfpTimes' restricts the analysis window for STL.
    stl = spklfp_stl('basepath', basepath, 'fs', fs, ...
        'sig', sig, 'spktimes', spkTimes, ...                   % Pass current sig and relative spkTimes
        'mapWin', [-0.1 0.1], 'lfpTimes', lfpTimes, ...         % Use same lfpTimes
        'graphics', false, 'flgSave', false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORGANIZE OUTPUT STRUCTURE AND SAVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section consolidates all analysis results into 'spkLfp' and saves it
% if 'flgSave' is true. The structure includes metadata and coupling metrics.

spkLfp = struct();
spkLfp.info.runtime         = datetime("now"); % Analysis timestamp.
spkLfp.info.fRange          = fRange;          % Frequency range [Hz].
spkLfp.info.lfpFs           = fs;              % LFP sampling rate [Hz].
spkLfp.info.lfpTimes        = lfpTimes;        % Analysis windows [s].
spkLfp.info.powThr          = powThr;          % Power threshold [STD].
spkLfp.info.nUnits          = nUnits;          % Number of units.
spkLfp.info.sigOffset       = sigOffset;       % LFP signal offset if loaded [s].
spkLfp.info.basepath        = basepath;        % Data basepath.

spkLfp.phase                = phase;           % Single unit phase coupling.
spkLfp.phase.spks           = spkPhase;        % LFP phase at each spike.
spkLfp.rateMap              = rateMap;         % Power-phase rate maps.
spkLfp.rateMag              = rateMag;         % Rate-magnitude correlation.
spkLfp.pop                  = pop;             % Population synchrony.
spkLfp.stl                  = stl;             % Spike-triggered LFP.

if flgSave
    % Save results to .mat file (includes frequency band in filename).
    save(spkLfpFile, 'spkLfp', '-v7.3');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphics (OPTIONAL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section generates and saves summary plots if 'flgGraphics' is true,
% using the 'spklfp_plot' function.

if flgGraphics 
    spklfp_plot(spkLfp); % Generate plots from results structure.
end

end     % EOF


