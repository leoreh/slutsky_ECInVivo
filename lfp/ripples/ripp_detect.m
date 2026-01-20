function ripp = ripp_detect(varargin)
% RIPP_DETECT Detects hippocampal ripple events from LFP data using adaptive thresholding.
%
% SUMMARY:
% This function identifies ripple events in LFP signals. Key steps include:
%   1. Loading LFP data (if not provided) or using provided LFP.
%   2. Filtering the LFP in the specified ripple band.
%   3. Calculating a detection signal (smoothed amplitude, squared signal, or TEO).
%   4. Applying an adaptive threshold based on local signal statistics.
%   5. Refining events using multiple criteria (peak power, continuous power duration, artifact rejection using optional provided EMG).
%   6. Calculating ripple properties (peak time, power, duration, frequency).
%   7. Generating peri-event maps and calculating statistics (ACG, rate, correlations).
%   8. Optionally suggesting refined detection parameters based on results.
%
% INPUT:
%   basepath        (Optional) Path to recording session directory {pwd}.
%   sig             (Optional) Numeric vector of LFP data. If empty, loads from basename.lfp.
%   emg             (Optional) Numeric vector of EMG data (assumed same fs and length as sig) for artifact rejection.
%   fs              (Optional) LFP sampling frequency [Hz]. If empty, loads from session info.
%   rippCh          (Optional) Channel(s) to load/average for ripple detection {uses session info}.
%   recWin          (Optional) [start end] time window to analyze [s] {[0 Inf]}.
%   flgPlot         (Optional) Logical flag to plot results {true}.
%   flgSave         (Optional) Logical flag to save results to .ripp.mat file {true}.
%   thr             (Optional) Numeric 5-element vector for thresholds {[1.5, 2.5, 2, 200, 100]}:
%                     thr(1): Initial detection (std above local mean).
%                     thr(2): Peak power (std above local mean).
%                     thr(3): Continuous power (std above local mean).
%                     thr(4): Artifact power (max allowed std).
%                     thr(5): EMG percentile for exclusion (requires 'emg' input).
%   passband        (Optional) [low high] frequency band for LFP filtering [Hz] {[100 300]}.
%   limDur          (Optional) Numeric 4-element vector for duration limits [ms] {[20, 300, 20, 10]}:
%                     limDur(1): Min ripple duration.
%                     limDur(2): Max ripple duration.
%                     limDur(3): Min inter-ripple interval.
%                     limDur(4): Min continuous power duration (above thr(3)).
%   detectAlt       (Optional) Detection method {3}: 1=Amp, 2=SqSig, 3=TEO.
%   bit2uv          (Optional) bit2uv conversion factor for LFP data {0.195}.
%
% OUTPUT:
%   ripp            Structure (see ripp_initialize).
%
% DEPENDENCIES:
%   CE_sessionTemplate (custom), bz_LoadBinary (buzcode), filterLFP (custom),
%   binary2bouts (custom), Sync (buzcode), SyncMap (buzcode), CCG (FMA),
%   times2rate (custom), plot_ripples (custom), ripp_initialize (local)
%
% HISTORY:
%   Aug 2024 LH - Major refactor and simplification.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING & INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defines input parameters, parses them using inputParser, and initializes
% the output structure.

% Define input parameters and their defaults/validation functions
p = inputParser;
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'sig', [], @isnumeric);      % LFP signal input
addParameter(p, 'emg', [], @isnumeric);      % EMG signal input (must match sig)
addParameter(p, 'recWin', [0 Inf], @(x) isnumeric(x) && numel(x)==2);
addParameter(p, 'rippCh', [], @isnumeric);    % Ripple channel(s) if loading LFP
addParameter(p, 'fs', [], @isnumeric);        % Sampling rate
addParameter(p, 'flgPlot', true, @islogical);
addParameter(p, 'flgSave', true, @islogical);
addParameter(p, 'thr', [1.5, 2.5, 2, 200, 100], @(x) isnumeric(x) && numel(x)==5);
addParameter(p, 'passband', [100 300], @(x) isnumeric(x) && numel(x)==2);
addParameter(p, 'limDur', [20, 300, 20, 10], @(x) isnumeric(x) && numel(x)==4);
addParameter(p, 'detectAlt', 3, @(x) isnumeric(x) && isscalar(x) && ismember(x, [1, 2, 3]));
addParameter(p, 'bit2uv', 0.195, @(x) isnumeric(x) && isscalar(x));

% Parse input arguments
parse(p, varargin{:});
basepath    = p.Results.basepath;
sig         = p.Results.sig;
emg         = p.Results.emg;
recWin      = p.Results.recWin;
rippCh      = p.Results.rippCh;
fs          = p.Results.fs;
flgPlot     = p.Results.flgPlot;
flgSave     = p.Results.flgSave;
thr         = p.Results.thr;
passband    = p.Results.passband;
limDur      = p.Results.limDur;
detectAlt   = p.Results.detectAlt;
bit2uv      = p.Results.bit2uv;

% Initialize output structure
ripp = ripp_initialize();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SESSION & PARAMETER SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sets up file paths, loads session information, determines sampling rates,
% converts time durations to samples, and defines calculation constants.

% Set basepath and define file paths
cd(basepath);
[~, basename] = fileparts(basepath);
sessionfile = fullfile(basepath, [basename, '.session.mat']);
rippfile    = fullfile(basepath, [basename, '.ripp.mat']);
clufile     = fullfile(basepath, [basename, '.ripp.clu.1']);
resfile     = fullfile(basepath, [basename, '.ripp.res.1']);

% Load session info, generating if necessary
if ~exist(sessionfile, 'file')
    warning('Session file not found: %s. Running CE_sessionTemplate.', sessionfile);
    session = CE_sessionTemplate(pwd, 'which', 'all', 'showGUI', false, ...
        'trigger', 'N', 'updateMode', 'update', 'makeData', false, 'saveMat', true);
else
    load(sessionfile, 'session');
end

% Assign default ripple channel if not provided and update session file
if isempty(sig)
    if isempty(rippCh)
        error('Ripple channel(s) must be provided.');
    else
        session.channelTags.Ripple = rippCh;
        save(sessionfile, 'session');
    end
end

% Determine LFP sampling frequency (fs) if not provided
if isempty(fs)
    fs = session.extracellular.srLfp;
end
fsSpks = session.extracellular.sr; % Spike sampling frequency
nchans = session.extracellular.nChannels;

% Convert duration limits from ms to samples
limDur_Samples = round(limDur / 1000 * fs);

% Define window length for local statistics (adaptive thresholding)
movLen = round(10 * fs); % 10 second moving window

% Define binsize for rate calculation
ripp.rate.binsize = 60; % [s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA LOADING & PREPARATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loads LFP signal from binary file if not provided directly, validates
% the time window, and checks consistency of optional EMG signal.

% Load LFP signal if 'sig' is empty
loadDur = Inf;
if isempty(sig)
    fprintf('Loading LFP data...');
    if isfinite(recWin(2))
        loadDur = recWin(2) - recWin(1);
    end

    fname = fullfile(basepath, [basename, '.lfp']);
    sig = binary_load(fname, 'duration', loadDur,...
        'fs', fs, 'nCh', nchans, 'start', recWin(1),...
        'ch', rippCh, 'downsample', 1, 'bit2uv', bit2uv);

    % Average if multiple ripple channels provided
    if size(sig, 2) > 1
        sig = mean(sig, 2);
    end
end

% Adjust recording window if loaded duration is shorter than requested
actualDur = size(sig,1)/fs;
if actualDur < loadDur - 1/fs
    recWin(2) = recWin(1) + actualDur;
    warning('Requested recording window truncated to actual LFP duration: [%.2f %.2f] s', recWin(1), recWin(2));
end

% Create timestamps vector
timestamps = recWin(1) + (0:length(sig)-1)' / fs;

% Validate EMG length if provided
if ~isempty(emg) && length(emg) ~= length(sig)
    error('EMG length (%d) does not match LFP length (%d).', length(emg), length(sig));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIGNAL PROCESSING FOR DETECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filters the LFP signal, calculates instantaneous amplitude, phase, and
% frequency, computes the chosen detection signal (Amp, SqSig, or TEO),
% and calculates the adaptive threshold signal based on local statistics.

% Filter LFP in the specified ripple passband
sigFilt = filterLFP(sig, 'fs', fs, 'type', 'butter', 'dataOnly', true,...
    'order', 5, 'passband', passband, 'graphics', false);

% Calculate Hilbert transform components (amplitude and phase)
sigH = hilbert(sigFilt);
sigAmp = abs(sigH);
sigPhase = angle(sigH);
sigUnwrapped = unwrap(sigPhase);

% Calculate instantaneous frequency
dt = diff(timestamps);
t_diff = timestamps(1 : end - 1) + dt / 2; % Midpoints for differentiation
d0 = diff(medfilt1(sigUnwrapped, 12)) ./ dt; % Differentiate phase after median filtering
d1 = interp1(t_diff, d0, timestamps(2 : end - 1, 1)); % Interpolate back to original timepoints
sigFreq = [d0(1); d1; d0(end)] / (2 * pi); % Combine and convert rad/s to Hz

% Select and calculate base detection signal based on 'detectAlt'
switch detectAlt
    case 1 % Smoothed Amplitude
        baseSignal = sigAmp;
        winSmooth = round(5 / (1000 / fs)); % 5ms smoothing window

    case 2 % Smoothed Squared Signal
        baseSignal = sigFilt .^ 2;
        winSmooth = round(10 / (1000 / fs)); % 10ms smoothing window

    case 3 % Smoothed Teager Energy Operator (TEO)
        sigFiltPad = [sigFilt(1); sigFilt; sigFilt(end)]; % Pad for TEO calculation
        baseSignal = sigFiltPad(2:end-1).^2 - sigFiltPad(1:end-2) .* sigFiltPad(3:end);
        baseSignal(baseSignal < 0) = 0; % Rectify TEO signal
        winSmooth = round(5 / (1000 / fs)); % 5ms smoothing window

    otherwise
        error('Unknown detection method (detectAlt).');
end
baseSignal = smoothdata(baseSignal, 'gaussian', winSmooth); % Smooth the chosen signal

% Calculate local mean and standard deviation
localMean = movmean(baseSignal, movLen);
localStd = movstd(baseSignal, movLen);
stdFloor = 1e-6 * mean(localStd, 'omitnan'); % Prevent division by zero or near-zero std
localStd(localStd < stdFloor) = stdFloor;

% Compute the final detection signal (z-scored relative to local stats)
sigDetect = (baseSignal - localMean) ./ localStd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RIPPLE EVENT DETECTION & REFINEMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs initial event detection based on the primary threshold, then
% refines events by excluding based on optional EMG activity, peak power
% limits, and minimum continuous duration above a secondary threshold.

% Initial detection of candidate events
eventSamples = binary2bouts('vec', sigDetect > thr(1), 'minDur', limDur_Samples(1),...
    'maxDur', limDur_Samples(2), 'interDur', limDur_Samples(3));
nEvents = size(eventSamples, 1);

% Initialize discard index and intermediate variables
idxDiscard = false(nEvents, 1); % Master discard index
emgRms = nan(nEvents, 1);
peakPower = nan(nEvents, 1);
contPowDur = nan(nEvents, 1);
contPowAvg = nan(nEvents, 1);
avgFreq = nan(nEvents, 1);

% EMG
if ~isempty(emg)

    % calculate EMG RMS per event
    for iEvent = 1:nEvents
        emgRms(iEvent) = rms(emg(eventSamples(iEvent, 1) : eventSamples(iEvent, 2)));
    end

    % Apply EMG exclusion
    idxDiscard = idxDiscard | (emgRms > prctile(emgRms, thr(5)));
end

% Remove events with high average frequency, potentially indicative of EMG
% artifacts. The frequency threshold is dynamically determined based on
% events coinciding with high EMG activity.
if ~isempty(emg)

    % Calculate mean instantaneous frequency per event
    for iEvent = 1:nEvents
        avgFreq(iEvent) = mean(sigFreq(eventSamples(iEvent, 1) : eventSamples(iEvent, 2)), 'omitnan');
    end

    % calculate the median mean frequency during events of high emg
    emgIdx = emgRms > prctile(emgRms, 99);
    thrFreq = mean(avgFreq(emgIdx));
    idxDiscard = idxDiscard | (avgFreq > thrFreq);
end

% Apply Peak Power and Continuous Duration Thresholds
for iEvent = 1:nEvents

    % Skip previously discarded events
    if idxDiscard(iEvent)
        continue
    end

    % Extract detection signal during the event
    eventPow = sigDetect(eventSamples(iEvent, 1) : eventSamples(iEvent, 2));

    % Check peak power threshold (minimum and maximum)
    peakPower(iEvent) = max(eventPow);
    if peakPower(iEvent) < thr(2) || peakPower(iEvent) > thr(4)
        idxDiscard(iEvent) = true;
        continue;
    end

    % Check continuous duration above secondary threshold (thr(3))
    aboveThr = eventPow > thr(3);
    if ~any(aboveThr)
        contPowDur(iEvent) = 0; % No continuous power above threshold
        contPowAvg(iEvent) = NaN;
    else
        % Find longest continuous run above threshold
        switches = diff([0; aboveThr(:); 0]);
        startRun = find(switches == 1);
        endRun = find(switches == -1) - 1;
        runLengths = endRun - startRun + 1;
        [contPowDur(iEvent), maxIdx] = max(runLengths);
        % Calculate average power during that longest run
        contPowAvg(iEvent) = mean(eventPow(startRun(maxIdx):endRun(maxIdx)));
    end

    % Discard if longest continuous run is too short
    if contPowDur(iEvent) < limDur_Samples(4)
        idxDiscard(iEvent) = true;
    end
end

% Remove discarded events and associated data
eventSamples(idxDiscard, :) = [];
avgFreq(idxDiscard, :) = [];
peakPower(idxDiscard) = [];
emgRms(idxDiscard) = [];
contPowDur(idxDiscard) = [];
contPowAvg(idxDiscard) = [];

% Check if any events remain
nEvents = size(eventSamples, 1);
if nEvents == 0
    warning('No valid ripples detected in %s.', basename);
    return;
else
    fprintf('Final number of ripple events for %s: %d\n', basename, nEvents);
end

% Clear intermediate detection variables to free memory
clear sigDetect aboveThr sigH sigUnwrapped baseSignal localMean localStd sigFiltPad eventPow switches startRun endRun runLengths maxIdx stdFloor


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE FINAL RIPPLE PROPERTIES & STATISTICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates properties (peak time, filtered peak, duration, frequency, etc.)
% for the refined set of ripple events, computes peri-event maps, ACG,
% firing rate, and correlations between ripple features.

% Find peak sample index within each ripple (minimum of filtered LFP)
peakSample = zeros(nEvents, 1);
for iEvent = 1:nEvents
    eventFilt = sigFilt(eventSamples(iEvent, 1) : eventSamples(iEvent, 2));
    [~, peakIdx] = min(eventFilt); % Peak defined as minimum of filtered LFP trough
    peakSample(iEvent) = eventSamples(iEvent, 1) + peakIdx - 1;
end

% Convert event boundaries and peak samples to times
peakTime = timestamps(peakSample);
times = timestamps(eventSamples);
dur = times(:, 2) - times(:, 1); % Ripple duration in seconds

% Generate Peri-Event Maps using Sync/SyncMap
nbinsMap = floor(fs * diff(ripp.maps.durWin) / 2) * 2 + 1; % Ensure odd number of bins
centerBin = ceil(nbinsMap / 2);
mapArgs = {'durations', ripp.maps.durWin, 'nbins', nbinsMap, 'smooth', 0};

[r, ir] = Sync([timestamps sig], peakTime, 'durations', ripp.maps.durWin);
ripp.maps.raw = SyncMap(r, ir, mapArgs{:}); % Raw LFP map

[r, ir] = Sync([timestamps sigFilt], peakTime, 'durations', ripp.maps.durWin);
ripp.maps.filt = SyncMap(r, ir, mapArgs{:}); % Filtered LFP map

[r, ir] = Sync([timestamps sigAmp], peakTime, 'durations', ripp.maps.durWin);
ripp.maps.amp = SyncMap(r, ir, mapArgs{:}); % Amplitude map

[r, ir] = Sync([timestamps sigPhase], peakTime, 'durations', ripp.maps.durWin);
ripp.maps.phase = SyncMap(r, ir, mapArgs{:}); % Phase map

[r, ir] = Sync([timestamps sigFreq], peakTime, 'durations', ripp.maps.durWin);
ripp.maps.freq = SyncMap(r, ir, mapArgs{:}); % Frequency map

% Extract ripple properties from maps or calculated previously
peakAmp = ripp.maps.amp(:, centerBin);                 % Amplitude at peak time
peakFreq = ripp.maps.freq(:, centerBin);               % Frequency at peak time
peakFilt = ripp.maps.filt(:, centerBin);               % Filtered LFP value at peak time
maxFreq = max(ripp.maps.freq, [], 2, 'omitnan');       % Max frequency during ripple map window

% Calculate range of filtered LFP around peak (+/- 5 samples)
idxWin = centerBin + (-5 : 5);
ripp.peakRng = range(ripp.maps.filt(:, idxWin), 2);

% Calculate peakEnergy using Mean Square (RMS^2) of filtered LFP around peak
ripp.peakEnergy = mean(ripp.maps.filt(:, idxWin).^2, 2);

% Calculate Autocorrelogram (ACG)
[ripp.acg.data, ripp.acg.t] = CCG(peakTime, ones(nEvents, 1),...
    'binSize', 0.01, 'duration', 1);                   % 1s ACG, 10ms bins

% Calculate Instantaneous Ripple Rate
[ripp.rate.rate, ripp.rate.binedges, ripp.rate.timestamps] = ...
    times2rate(peakTime, 'binsize', ripp.rate.binsize, 'winCalc', recWin, 'c2r', true);

% Calculate Correlations between ripple properties
ripp.corr.AmpFreq = corr(peakAmp, peakFreq, 'rows', 'complete');
ripp.corr.DurFreq = corr(dur, peakFreq, 'rows', 'complete');
ripp.corr.DurAmp = corr(dur, peakAmp, 'rows', 'complete');

% Clear large signal arrays to free memory
clear sig sigFilt sigAmp sigPhase sigFreq timestamps r ir eventFilt peakIdx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER TUNING SUGGESTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyzes the detected ripple properties (peak power, continuous power avg,
% continuous power length) to suggest potentially refined threshold and
% duration parameters based on percentiles of the detected events.

% Suggest power thresholds
ripp.info.thrData(1) = prctile(contPowAvg, 1);
ripp.info.thrData(2) = prctile(peakPower, 10);
ripp.info.thrData(3) = prctile(contPowAvg, 10);
ripp.info.thrData(4) = prctile(peakPower, 99) * 2;
ripp.info.thrData(5) = thr(5); % Keep original EMG threshold suggestion
ripp.info.thrData = round(ripp.info.thrData, 1);

% Suggest duration limits
ripp.info.limDurData = limDur; % Start with original limits
ripp.info.limDurData(4) = round(prctile(contPowDur, 5) / fs * 1000);

% Suggest passband limits
ripp.info.passbandData = [prctile(avgFreq, 2), prctile(peakFreq, 99.9)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINALIZE AND SAVE OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Populates the final output structure 'ripp' with all calculated data
% and metadata, saves the structure to a .ripp.mat file, optionally creates
% .res and .clu files for NeuroScope, and generates summary plots if requested.

% Populate Output Structure
ripp                = workspace2struct(ripp);
ripp.info           = workspace2struct(ripp.info);

% update specific fields
ripp.contPowDur     = contPowDur / fs * 1000;
ripp.info.runtime   = datetime("now");
ripp.contPowDur     = contPowDur / fs * 1000;
ripp.dur            = ripp.dur * 1000;         % covert to ms

% Save Results
if flgSave
    save(rippfile, 'ripp', '-v7.3');

    % Create .res and .clu files for NeuroScope visualization
    % Note: This currently assumes recWin(1) = 0 for correct sample
    % alignment. In addition, for tdt (non-integer fsSpks), there is a
    % drift in the timing of events (only during visualization).
    if recWin(1) ~= 0
        warning('Neuroscope .res/.clu file generation might be inaccurate because recWin does not start at 0.');
    end

    res = round([ripp.times(:, 1); ripp.times(:, 2); ripp.peakTime] * fsSpks);
    [res, sort_idx] = sort(res);
    fid = fopen(resfile, 'w');
    fprintf(fid, '%d\n', res);
    rc = fclose(fid);
    if rc == -1
        error('failed to close res file')
    end

    clu = [ones(nEvents, 1); ones(nEvents, 1) * 2; ones(nEvents, 1) * 3]; % 1=start, 2=end, 3=peak
    clu = clu(sort_idx);
    fid = fopen(clufile, 'w');
    fprintf(fid, '%d\n', 3);
    fprintf(fid, '%d\n', clu);
    rc = fclose(fid);
    if rc == -1
        error('failed to close clu file')
    end
end

% Plot Graphics
if flgPlot
    ripp_plot(ripp, 'basepath', basepath, 'flgSaveFig', true);
end

end % END FUNCTION ripp_detect

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION: INITIALIZE RIPPLE STRUCT
% Creates the initial empty ripple structure with predefined fields.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ripp = ripp_initialize()
% Creates the initial ripp structure with all fields set to default/empty.

ripp = struct();
% Event Timing & Core Properties
ripp.times              = []; % Nx2 matrix, [start end] times in seconds
ripp.dur                = []; % Nx1 vector, duration of ripples in (ms)
ripp.peakTime           = []; % Nx1 vector, time of ripple peak (min of filtered LFP) in seconds
ripp.peakFilt           = []; % Nx1 vector, value of filtered LFP at peak time
ripp.peakPower          = []; % Nx1 vector, value of detection signal (z-score) at its peak within the event
ripp.peakAmp            = []; % Nx1 vector, LFP envelope amplitude at peakTime
ripp.peakEnergy         = []; % Nx1 vector, energy (mean squared amplitude) of filtered LFP around peak
ripp.peakRng            = []; % Nx1 vector, range of filtered LFP around peak
ripp.peakFreq           = []; % Nx1 vector, instantaneous frequency at peakTime [Hz]
ripp.maxFreq            = []; % Nx1 vector, maximum instantaneous frequency within map window [Hz]
ripp.emgRms             = []; % Nx1 vector, RMS of EMG during ripple (NaN if EMG not provided/used)
ripp.contPowDur         = []; % Nx1 vector, length of longest continuous period above thr(3) [ms]
ripp.contPowAvg         = []; % Nx1 vector, average detection signal during longest continuous period above thr(3)

% Info substruct (Metadata)
ripp.info               = struct();
ripp.info.basename      = '';
ripp.info.rippCh        = [];
ripp.info.limDur        = []; % [minRippDur maxRippDur minInterRipp minContPowerDur] in ms
ripp.info.limDurData    = []; % Suggested duration limits based on detected data in ms
ripp.info.recWin        = []; % [start end] analysis window in seconds
ripp.info.runtime       = []; % datetime object
ripp.info.thr           = []; % [detect peak cont artifact emgPercentile] thresholds used
ripp.info.thrData       = []; % Suggested thresholds based on detected data
ripp.info.fs            = []; % LFP sampling rate in Hz
ripp.info.fsSpks        = []; % Spike sampling rate in Hz
ripp.info.passband      = []; % [low high] filter passband in Hz
ripp.info.passbandData  = []; % [low high] filter passband in Hz
ripp.info.detectAlt     = []; % Detection method: 1=Amp, 2=SqSig, 3=TEO
ripp.info.movLen        = []; % Moving window size for adaptive threshold in seconds
ripp.info.bit2uv        = []; % bit2uv conversion factor used

% Maps substruct (Peri-Event Averages)
ripp.maps               = struct();
ripp.maps.durWin        = [-75 75] / 1000; % Default map window [-0.075 +0.075] seconds
ripp.maps.raw           = []; % Map of raw LFP signal
ripp.maps.filt          = []; % Map of filtered LFP signal
ripp.maps.amp           = []; % Map of LFP amplitude envelope
ripp.maps.phase         = []; % Map of LFP phase
ripp.maps.freq          = []; % Map of instantaneous LFP frequency

% Rate substruct
ripp.rate               = struct();
ripp.rate.rate          = []; % Instantaneous ripple rate vector [Hz]
ripp.rate.binedges      = []; % Bin edges for rate calculation [s]
ripp.rate.timestamps    = []; % Timestamps for rate calculation [s]
ripp.rate.binsize       = []; % binsize for rate calculation [s]

% ACG substruct (Autocorrelogram)
ripp.acg                = struct();
ripp.acg.data           = []; % ACG counts per bin
ripp.acg.t              = []; % Time lags for ACG [s]

% Correlations substruct
ripp.corr               = struct();
ripp.corr.AmpFreq       = NaN; % Correlation between peak Amplitude and peak Frequency
ripp.corr.DurFreq       = NaN; % Correlation between Duration and peak Frequency
ripp.corr.DurAmp        = NaN; % Correlation between Duration and peak Amplitude

end     % EOF


%% ========================================================================
%  NOTE: PEAK AMPLITUDE VS. PEAK POWER (Z-SCORE)
%  ========================================================================
%  In ripple analysis, distinguishing between "Raw Magnitude" and
%  "Statistical Magnitude" is critical
%
%  1. Peak Amplitude (peakAmp) & Energy (peakEnergy):
%     These metrics represent the *physical* strength of the oscillation.
%     - 'peakAmp' is the envelope amplitude in microvolts (uV).
%     - 'peakEnergy' is the squared amplitude or RMS power (uV^2).
%     These values reflect the absolute magnitude of the synchronous
%     population firing. If a manipulation (like MCU-KO) increases neuronal
%     excitability or recruitment, these raw metrics will increase directly.
%
%  2. Peak Power (peakPower):
%     In the context of this detection algorithm, 'peakPower' is a Z-SCORE.
%     It does not measure absolute energy. Instead, it measures the
%     *Signal-to-Noise Ratio* (SNR) relative to the local background:
%
%         peakPower = (RawAmp - BaselineMean) / BaselineSD
%
%     This metric quantifies how much the ripple "pops out" from the
%     background activity. It is a measure of detection reliability, not
%     physiological strength.
%
%  THE "SNR PARADOX":
%     It is possible for ripples to become physically much larger (higher
%     peakAmp) while their Z-score (peakPower) remains unchanged. This
%     occurs if the background "noise" (BaselineSD) increases proportionally
%     with the signal.
%
%     Biologically, this may indicate a global increase in network gain.
%     The network is "noisier" or more active even during non-ripple
%     periods (higher baseline high-frequency activity), so the larger
%     ripples do not stand out any more clearly against this louder
%     background than smaller ripples do against a quiet background.
%
%  DECISION RULE:
%     - Use 'peakAmp' or 'peakEnergy' to test hypotheses about **Event
%       Magnitude** (e.g., "Does the manipulation make ripples stronger?").
%     - Use 'peakPower' (Z-score) to test hypotheses about **Event
%       Saliency** or detection quality (e.g., "Are ripples harder to detect
%       in this group?").
%  ========================================================================


