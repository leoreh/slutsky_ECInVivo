function ripp = ripp_detect(rippSig, fs, varargin)
% RIPP_DETECT Detects hippocampal ripple events from LFP data.
%
% SUMMARY:
% Uses provided signal 'sig' (assumed filtered) for event detection.
% Calculates detection signal, applies adaptive thresholds, refines events
% based on power/duration/EMG, and calculates event properties.
%
% INPUT:
%   rippSig         Structure with fields: lfp, filt, amp, freq, z, emg.
%   fs              Sampling frequency [Hz].
%   varargin        'basepath', 'thr', 'limDur', 'detectAlt', 'flgPlot', 'flgSave'.
%
% OUTPUT:
%   ripp            Structure (times, peakTime, info, maps, etc).
%
% DEPENDENCIES: binary2bouts, Sync, SyncMap, CCG, times2rate.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'rippSig', @isstruct);
addRequired(p, 'fs', @isnumeric);
addParameter(p, 'emg', [], @isnumeric);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'thr', [1.5, 2.5, 2, 200, 50], @isnumeric);
addParameter(p, 'limDur', [15, 300, 20, 10], @isnumeric);
addParameter(p, 'detectAlt', 3, @isnumeric);
addParameter(p, 'flgPlot', true, @islogical);
addParameter(p, 'flgSave', true, @islogical);

parse(p, rippSig, fs, varargin{:});
rippSig     = p.Results.rippSig;
fs          = p.Results.fs;
emg         = p.Results.emg;
basepath    = p.Results.basepath;
thr         = p.Results.thr;
limDur      = p.Results.limDur;
flgPlot     = p.Results.flgPlot;
flgSave     = p.Results.flgSave;



%% ========================================================================
%  SETUP
%  ========================================================================

% Initialize
ripp = struct();
ripp.times = []; ripp.dur = []; ripp.peakTime = [];
ripp.info = struct();
ripp.maps = struct('durWin', [-0.075 0.075]);
ripp.rate = struct();
ripp.acg = struct();

% File
cd(basepath);
[~, basename] = fileparts(basepath);
rippfile = fullfile(basepath, [basename, '.ripp.mat']);

timestamps = (0:length(rippSig.lfp)-1)' / fs;
ripp.rate.binsize = 60;

% Duration limits to samples
limDur_Samples = round(limDur / 1000 * fs);


%% ========================================================================
%  DETECTION
%  ========================================================================
eventSamples = binary2bouts('vec', rippSig.z > thr(1), 'minDur', limDur_Samples(1),...
    'maxDur', limDur_Samples(2), 'interDur', limDur_Samples(3));
nEvents = size(eventSamples, 1);

idxDiscard = false(nEvents, 1);
emgRms = nan(nEvents, 1);
peakPower = nan(nEvents, 1);
contPowDur = nan(nEvents, 1);
contPowAvg = nan(nEvents, 1);
avgFreq = nan(nEvents, 1);

% EMG Exclusion
% if ~isempty(emg)
%     for iEvent = 1:nEvents
%         idx = eventSamples(iEvent,1):eventSamples(iEvent,2);
%         emgRms(iEvent) = rms(emg(idx));
%     end
%     idxDiscard = idxDiscard | (emgRms > prctile(emgRms, thr(5)));
% end

% Frequency Exclusion
% if ~isempty(emg)
%     for iEvent = 1:nEvents
%         idx = eventSamples(iEvent,1):eventSamples(iEvent,2);
%         avgFreq(iEvent) = mean(rippSig.freq(idx), 'omitnan');
%     end
% 
%     highEmgIdx = emgRms > prctile(emgRms, 99);
%     if any(highEmgIdx)
%         thrFreq = mean(avgFreq(highEmgIdx));
%         idxDiscard = idxDiscard | (avgFreq > thrFreq);
%     end
% end

% Power and Continuity thresholds
for iEvent = 1:nEvents
    if idxDiscard(iEvent), continue; end

    idx = eventSamples(iEvent,1):eventSamples(iEvent,2);
    evtPow = rippSig.z(idx);

    % Peak Power
    peakPower(iEvent) = max(evtPow);
    if peakPower(iEvent) < thr(2) || peakPower(iEvent) > thr(4)
        idxDiscard(iEvent) = true;
        continue;
    end

    % Continuous Power
    above = evtPow > thr(3);
    if ~any(above)
        contPowDur(iEvent) = 0;
        contPowAvg(iEvent) = NaN;
    else
        d = [0; above(:); 0];
        on = find(diff(d)==1);
        off = find(diff(d)==-1);
        lens = off - on;
        [maxLen, maxIdx] = max(lens);

        contPowDur(iEvent) = maxLen;
        contPowAvg(iEvent) = mean(evtPow(on(maxIdx):off(maxIdx)-1));
    end

    if contPowDur(iEvent) < limDur_Samples(4)
        idxDiscard(iEvent) = true;
    end
end

% Cleanup discarded
eventSamples(idxDiscard, :) = [];
avgFreq(idxDiscard) = [];
peakPower(idxDiscard) = [];
contPowAvg(idxDiscard) = [];
contPowDur(idxDiscard) = [];
emgRms(idxDiscard) = [];

nEvents = size(eventSamples, 1);
if nEvents == 0
    warning('No valid ripples detected.');
    return;
else
    fprintf('Detected %d events.\n', nEvents);
end


%% ========================================================================
%  PROPERTIES & STATISTICS
%  ========================================================================

% Peak Finding
peakSample = zeros(nEvents, 1);
for iEvent = 1:nEvents
    idx = eventSamples(iEvent,1):eventSamples(iEvent,2);
    [~, peakRel] = min(rippSig.filt(idx));
    peakSample(iEvent) = eventSamples(iEvent,1) + peakRel - 1;
end

peakTime = timestamps(peakSample);
times = timestamps(eventSamples);
dur = times(:,2) - times(:,1);

% Maps
nbinsMap = floor(fs * diff(ripp.maps.durWin) / 2) * 2 + 1;
centerBin = ceil(nbinsMap / 2);
mapArgs = {'durations', ripp.maps.durWin, 'nbins', nbinsMap, 'smooth', 0};

[r, ir] = Sync([timestamps rippSig.filt], peakTime, 'durations', ripp.maps.durWin);
ripp.maps.raw = SyncMap(r, ir, mapArgs{:});
ripp.maps.filt = ripp.maps.raw;

[r, ir] = Sync([timestamps rippSig.amp], peakTime, 'durations', ripp.maps.durWin);
ripp.maps.amp = SyncMap(r, ir, mapArgs{:});

[r, ir] = Sync([timestamps rippSig.sigPhase], peakTime, 'durations', ripp.maps.durWin);
ripp.maps.phase = SyncMap(r, ir, mapArgs{:});

[r, ir] = Sync([timestamps rippSig.freq], peakTime, 'durations', ripp.maps.durWin);
ripp.maps.freq = SyncMap(r, ir, mapArgs{:});

% Store Properties
ripp.times = times;
ripp.peakTime = peakTime;
ripp.dur = dur;
ripp.peakAmp = ripp.maps.amp(:, centerBin);
ripp.peakFreq = ripp.maps.freq(:, centerBin);
ripp.peakFilt = ripp.maps.filt(:, centerBin);

idxWin = centerBin + (-5:5);
ripp.peakEnergy = mean(ripp.maps.filt(:, idxWin).^2, 2);

% ACG
[ripp.acg.data, ripp.acg.t] = CCG(peakTime, ones(nEvents, 1), 'binSize', 0.01, 'duration', 1);

% Rate
[ripp.rate.rate, ripp.rate.binedges, ripp.rate.timestamps] = ...
    times2rate(peakTime, 'binsize', ripp.rate.binsize, 'winCalc', [timestamps(1) timestamps(end)], 'c2r', true);


%% ========================================================================
%  PARAMETER SUGGESTIONS
%  ========================================================================
% Analyzes detected properties to suggest data-driven thresholds

ripp.info.thrData(1) = prctile(contPowAvg, 1);
ripp.info.thrData(2) = prctile(peakPower, 10);
ripp.info.thrData(3) = prctile(contPowAvg, 10);
ripp.info.thrData(4) = prctile(peakPower, 99) * 2;
ripp.info.thrData(5) = thr(5);
ripp.info.thrData = round(ripp.info.thrData, 1);

%% ========================================================================
%  FINALIZE
%  ========================================================================
ripp.contPowDur = contPowDur * 1000 / fs; % Convert to ms if storing
ripp.dur = ripp.dur * 1000; % ms

if flgSave
    save(rippfile, 'ripp', '-v7.3');
end

if flgPlot
    ripp_plot(ripp, 'basepath', basepath, 'flgSaveFig', true);
end

end


