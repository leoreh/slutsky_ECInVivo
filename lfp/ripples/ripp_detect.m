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
%   flgGraphics     (Optional) Logical flag to plot results {true}.
%   flgSaveVar      (Optional) Logical flag to save results to .ripp.mat file {true}.
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
%
% OUTPUT:
%   ripp            Structure (see initialize_ripple_struct).
%
% DEPENDENCIES:
%   CE_sessionTemplate (custom), bz_LoadBinary (buzcode), filterLFP (custom),
%   binary2bouts (custom), Sync (buzcode), SyncMap (buzcode), CCG (FMA),
%   times2rate (custom), plot_ripples (custom), initialize_ripple_struct (local)
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
addParameter(p, 'flgGraphics', true, @islogical);
addParameter(p, 'flgSaveVar', true, @islogical);
addParameter(p, 'thr', [1.5, 2.5, 2, 200, 100], @(x) isnumeric(x) && numel(x)==5);
addParameter(p, 'passband', [100 300], @(x) isnumeric(x) && numel(x)==2);
addParameter(p, 'limDur', [20, 300, 20, 10], @(x) isnumeric(x) && numel(x)==4);
addParameter(p, 'detectAlt', 3, @(x) isnumeric(x) && isscalar(x) && ismember(x, [1, 2, 3]));

% Parse input arguments
parse(p, varargin{:});
basepath    = p.Results.basepath;
sig         = p.Results.sig;
emg         = p.Results.emg;
recWin      = p.Results.recWin;
rippCh      = p.Results.rippCh;
fs          = p.Results.fs;
flgGraphics = p.Results.flgGraphics;
flgSaveVar  = p.Results.flgSaveVar;
thr         = p.Results.thr;
passband    = p.Results.passband;
limDur      = p.Results.limDur;
detectAlt   = p.Results.detectAlt;

% Initialize output structure
ripp = initialize_ripple_struct();

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
binsizeRate = 60; % [s]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA LOADING & PREPARATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loads LFP signal from binary file if not provided directly, validates
% the time window, and checks consistency of optional EMG signal.

% Load LFP signal if 'sig' parameter is empty
loadDur = Inf;
if isempty(sig)
    fprintf('Loading LFP data...');
    if isfinite(recWin(2))
        loadDur = recWin(2) - recWin(1);
    end

    sig = double(bz_LoadBinary(fullfile(basepath, [basename, '.lfp']), ...
        'duration', loadDur, 'frequency', fs, 'nchannels', nchans, ...
        'start', recWin(1), 'channels', rippCh, 'downsample', 1));

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
     error('Provided EMG length (%d) does not match LFP length (%d).', length(emg), length(sig));
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
h = hilbert(sigFilt);
sigAmp = abs(h);
sigPhase = angle(h);
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
rippEmg = nan(nEvents, 1);
peakPower = nan(nEvents, 1);
contPowerLength = nan(nEvents, 1);
contPowerAvg = nan(nEvents, 1);
avgFreq = nan(nEvents, 1);

% EMG 
if ~isempty(emg) 
    
    % calculate EMG RMS per event
    for iEvent = 1:nEvents
        rippEmg(iEvent) = rms(emg(eventSamples(iEvent, 1) : eventSamples(iEvent, 2)));
    end
    
    % Apply EMG exclusion
    idxDiscard = idxDiscard | (rippEmg > prctile(rippEmg, thr(5)));
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
    emgIdx = rippEmg > prctile(rippEmg, 99);
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
        contPowerLength(iEvent) = 0; % No continuous power above threshold
        contPowerAvg(iEvent) = NaN;
    else
        % Find longest continuous run above threshold
        switches = diff([0; aboveThr(:); 0]);
        startRun = find(switches == 1);
        endRun = find(switches == -1) - 1;
        runLengths = endRun - startRun + 1;
        [contPowerLength(iEvent), maxIdx] = max(runLengths);
        % Calculate average power during that longest run
        contPowerAvg(iEvent) = mean(eventPow(startRun(maxIdx):endRun(maxIdx)));
    end

    % Discard if longest continuous run is too short
    if contPowerLength(iEvent) < limDur_Samples(4)
        idxDiscard(iEvent) = true;
    end
end

% Remove discarded events and associated data
eventSamples(idxDiscard, :) = [];
avgFreq(idxDiscard, :) = [];
peakPower(idxDiscard) = [];
rippEmg(idxDiscard) = [];
contPowerLength(idxDiscard) = [];
contPowerAvg(idxDiscard) = [];

% Check if any events remain
nEvents = size(eventSamples, 1);
if nEvents == 0
    warning('No valid ripples detected after filtering for %s. Returning empty struct.', basename);
    % Keep ripp struct initialized but empty
    return;
else
    fprintf(' Final number of ripple events for %s: %d\n', basename, nEvents);
end

% Clear intermediate detection variables to free memory
clear sigDetect aboveThr h sigUnwrapped baseSignal localMean localStd sigFiltPad eventPow switches startRun endRun runLengths maxIdx stdFloor


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
mapDurWin = ripp.maps.durWin; % Time window for maps [s]
nbinsMap = floor(fs * diff(mapDurWin) / 2) * 2 + 1; % Ensure odd number of bins
centerBin = ceil(nbinsMap / 2);
mapArgs = {'durations', mapDurWin, 'nbins', nbinsMap, 'smooth', 0};

[r, ir] = Sync([timestamps sig], peakTime, 'durations', mapDurWin);
mapSig = SyncMap(r, ir, mapArgs{:}); % Raw LFP map

[r, ir] = Sync([timestamps sigFilt], peakTime, 'durations', mapDurWin);
mapFilt = SyncMap(r, ir, mapArgs{:}); % Filtered LFP map

[r, ir] = Sync([timestamps sigAmp], peakTime, 'durations', mapDurWin);
mapAmp = SyncMap(r, ir, mapArgs{:}); % Amplitude map

[r, ir] = Sync([timestamps sigPhase], peakTime, 'durations', mapDurWin);
mapPhase = SyncMap(r, ir, mapArgs{:}); % Phase map

[r, ir] = Sync([timestamps sigFreq], peakTime, 'durations', mapDurWin);
mapFreq = SyncMap(r, ir, mapArgs{:}); % Frequency map

% Extract ripple properties from maps or calculated previously
peakAmp = mapAmp(:, centerBin);                 % Amplitude at peak time
peakFreq = mapFreq(:, centerBin);               % Frequency at peak time
peakFilt = mapFilt(:, centerBin);               % Filtered LFP value at peak time
maxFreq = max(mapFreq, [], 2, 'omitnan');       % Max frequency during ripple map window

% Calculate Autocorrelogram (ACG)
[acgData, acgT] = CCG(peakTime, ones(nEvents, 1), 'binSize', 0.01, 'duration', 1); % 1s ACG, 10ms bins

% Calculate Instantaneous Ripple Rate
[rateRate, rateBins, rateTimestamps] = ...
    times2rate(peakTime, 'binsize', binsizeRate, 'winCalc', recWin, 'c2r', true);

% Calculate Correlations between ripple properties
rippCorr.AmpFreq = corr(peakAmp, peakFreq, 'rows', 'complete');
rippCorr.DurFreq = corr(dur, peakFreq, 'rows', 'complete');
rippCorr.DurAmp = corr(dur, peakAmp, 'rows', 'complete');

% Clear large signal arrays to free memory
clear sig sigFilt sigAmp sigPhase sigFreq timestamps r ir eventFilt peakIdx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER TUNING SUGGESTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyzes the detected ripple properties (peak power, continuous power avg,
% continuous power length) to suggest potentially refined threshold and
% duration parameters based on percentiles of the detected events.

% Suggest power thresholds 
ripp.info.thrData(1) = prctile(contPowerAvg, 1);   
ripp.info.thrData(2) = prctile(peakPower, 10);  
ripp.info.thrData(3) = prctile(contPowerAvg, 10); 
ripp.info.thrData(4) = prctile(peakPower, 99) * 2; 
ripp.info.thrData(5) = thr(5); % Keep original EMG threshold suggestion
ripp.info.thrData = round(ripp.info.thrData, 1);

% Suggest duration limits 
ripp.info.limDurData = limDur; % Start with original limits
ripp.info.limDurData(4) = round(prctile(contPowerLength, 5) / fs * 1000); 

% Suggest passband limits 
ripp.info.passbandData = [prctile(avgFreq, 2), prctile(peakFreq, 99.9)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION 8: FINALIZE OUTPUT, SAVE RESULTS, & PLOT GRAPHICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Populates the final output structure 'ripp' with all calculated data
% and metadata, saves the structure to a .ripp.mat file, optionally creates
% .res and .clu files for NeuroScope, and generates summary plots if requested.

% Populate Output Structure ('ripp')
ripp.times              = times;           % [startTime endTime]
ripp.peakTime           = peakTime;         % Time of ripple peak (min of filtered LFP) in seconds
ripp.peakFilt           = peakFilt;         % Filtered LFP value at peak time
ripp.peakPower          = peakPower;        % Peak detection signal value (z-score) at its peak within the event
ripp.dur                = dur;              % Duration in seconds
ripp.rippEmg            = rippEmg;          % RMS of EMG during ripple (NaN if EMG not provided/used)
ripp.contPowerLength    = contPowerLength / fs * 1000; % Length of longest continuous run above thr(3) [ms]
ripp.contPowerAvg       = contPowerAvg;     % Average detection signal during longest run above thr(3)
ripp.peakFreq           = peakFreq;         % Instantaneous frequency at peak time [Hz]
ripp.peakAmp            = peakAmp;          % LFP envelope amplitude at peak time
ripp.maxFreq            = maxFreq;          % Max instantaneous frequency within map window [Hz]
ripp.avgFreq            = avgFreq;
ripp.maps.raw           = mapSig;           % Peri-event map: Raw LFP
ripp.maps.ripp          = mapFilt;          % Peri-event map: Filtered LFP
ripp.maps.amp           = mapAmp;           % Peri-event map: Amplitude
ripp.maps.phase         = mapPhase;         % Peri-event map: Phase
ripp.maps.freq          = mapFreq;          % Peri-event map: Frequency
ripp.rate.rate          = rateRate;         % Instantaneous ripple rate [Hz]
ripp.rate.binedges      = rateBins;         % Bin edges for rate calculation [s]
ripp.rate.timestamps    = rateTimestamps;   % Timestamps for rate calculation [s]
ripp.acg.data           = acgData;          % Autocorrelogram data
ripp.acg.t              = acgT;             % Time lags for ACG [s]
ripp.corr               = rippCorr;         % Correlations between ripple features
ripp.info.basename      = basename;         % Basename of the session
ripp.info.rippCh        = rippCh;           % Channel(s) used for detection
ripp.info.limDur        = limDur;           % Original duration limits used [ms]
ripp.info.recWin        = recWin;           % Recording window analyzed [s]
ripp.info.runtime       = datetime("now");  % Timestamp of detection run
ripp.info.thr           = thr;              % Original thresholds used
ripp.info.binsizeRate   = binsizeRate;      % Binsize used for rate calc [s]
ripp.info.fs            = fs;               % LFP sampling frequency [Hz]
ripp.info.fsSpk         = fsSpks;           % Spike sampling frequency [Hz]
ripp.info.passband      = passband;         % Filter passband used [Hz]
ripp.info.detectAlt     = detectAlt;        % Detection method used (1=Amp, 2=SqSig, 3=TEO)
ripp.info.movLen        = movLen / fs;      % Moving window size for adaptive threshold [s]

% Save Results 
if flgSaveVar
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
if flgGraphics
    ripp_plot(ripp, 'basepath', basepath, 'flgSaveFig', true);
end

end % END FUNCTION ripp_detect

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION: INITIALIZE RIPPLE STRUCT
% Creates the initial empty ripple structure with predefined fields.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ripp = initialize_ripple_struct()
% Creates the initial ripp structure with all fields set to default/empty.

    ripp = struct();
    % Event Timing & Core Properties
    ripp.times             = []; % Nx2 matrix, [start end] times in seconds
    ripp.peakTime           = []; % Nx1 vector, time of ripple peak (min of filtered LFP) in seconds
    ripp.peakFilt           = []; % Nx1 vector, value of filtered LFP at peak time
    ripp.peakPower          = []; % Nx1 vector, value of detection signal (z-score) at its peak within the event
    ripp.dur                = []; % Nx1 vector, duration of ripples in seconds
    ripp.rippEmg            = []; % Nx1 vector, RMS of EMG during ripple (NaN if EMG not provided/used)
    ripp.contPowerLength    = []; % Nx1 vector, length of longest continuous period above thr(3) [ms]
    ripp.contPowerAvg       = []; % Nx1 vector, average detection signal during longest continuous period above thr(3)
    ripp.peakFreq           = []; % Nx1 vector, instantaneous frequency at peakTime [Hz]
    ripp.peakAmp            = []; % Nx1 vector, LFP envelope amplitude at peakTime
    ripp.maxFreq            = []; % Nx1 vector, maximum instantaneous frequency within map window [Hz]

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
    ripp.info.binsizeRate   = []; % binsize for rate calculation in seconds
    ripp.info.fs            = []; % LFP sampling rate in Hz
    ripp.info.fsSpk         = []; % Spike sampling rate in Hz
    ripp.info.passband      = []; % [low high] filter passband in Hz
    ripp.info.passbandData  = []; % [low high] filter passband in Hz
    ripp.info.detectAlt     = []; % Detection method: 1=Amp, 2=SqSig, 3=TEO
    ripp.info.movLen        = []; % Moving window size for adaptive threshold in seconds

    % Maps substruct (Peri-Event Averages)
    ripp.maps               = struct();
    ripp.maps.durWin        = [-75 75] / 1000; % Default map window [-0.075 +0.075] seconds
    ripp.maps.raw           = []; % Map of raw LFP signal
    ripp.maps.ripp          = []; % Map of filtered LFP signal
    ripp.maps.amp           = []; % Map of LFP amplitude envelope
    ripp.maps.phase         = []; % Map of LFP phase
    ripp.maps.freq          = []; % Map of instantaneous LFP frequency

    % Rate substruct
    ripp.rate               = struct();
    ripp.rate.rate          = []; % Instantaneous ripple rate vector [Hz]
    ripp.rate.binedges      = []; % Bin edges for rate calculation [s]
    ripp.rate.timestamps    = []; % Timestamps for rate calculation [s]

    % ACG substruct (Autocorrelogram)
    ripp.acg                = struct();
    ripp.acg.data           = []; % ACG counts per bin
    ripp.acg.t              = []; % Time lags for ACG [s]

    % Correlations substruct
    ripp.corr               = struct();
    ripp.corr.AmpFreq       = NaN; % Correlation between peak Amplitude and peak Frequency
    ripp.corr.DurFreq       = NaN; % Correlation between Duration and peak Frequency
    ripp.corr.DurAmp        = NaN; % Correlation between Duration and peak Amplitude

end
% EOF