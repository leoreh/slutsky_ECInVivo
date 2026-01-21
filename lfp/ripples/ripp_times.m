function ripp = ripp_times(rippSig, fs, varargin)
% RIPP_TIMES Preliminary detection of hippocampal ripple events.
%
% SUMMARY:
% Detects potential ripple events based on thresholding the Z-scored signal.
% Applies duration and power constraints. Returns times and peak times.
%
% INPUT:
%   rippSig         Structure with fields: lfp, filt, amp, freq, z
%   fs              Sampling frequency [Hz].
%   varargin        'thr', 'limDur'
%
% OUTPUT:
%   ripp            Structure with fields:
%                   .times      (nEvents x 2) start and end times [s]
%                   .peakTime   (nEvents x 1) peak time [s]
%
% DEPENDENCIES: binary2bouts

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'rippSig', @isstruct);
addRequired(p, 'fs', @isnumeric);
addParameter(p, 'thr', [1.5, 2.5, 2, 200, 50], @isnumeric);
addParameter(p, 'limDur', [15, 300, 20, 10], @isnumeric);

parse(p, rippSig, fs, varargin{:});
rippSig     = p.Results.rippSig;
fs          = p.Results.fs;
thr         = p.Results.thr;
limDur      = p.Results.limDur;

%% ========================================================================
%  SETUP
%  ========================================================================
ripp = struct();
ripp.times = [];
ripp.peakTime = [];

timestamps = (0:length(rippSig.lfp)-1)' / fs;

% Duration limits in samples
limDur_Samples = round(limDur / 1000 * fs);

%% ========================================================================
%  DETECTION
%  ========================================================================
% Initial bout detection
eventSamples = binary2bouts('vec', rippSig.z > thr(1), 'minDur', limDur_Samples(1),...
    'maxDur', limDur_Samples(2), 'interDur', limDur_Samples(3));
nEvents = size(eventSamples, 1);

idxDiscard = false(nEvents, 1);
emgRms = nan(nEvents, 1);
peakPower = nan(nEvents, 1);
contPowDur = nan(nEvents, 1);
contPowAvg = nan(nEvents, 1);

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

% Cleanup discarded events
eventSamples(idxDiscard, :) = [];

% (Unused stats are cleared implicitly by not returning them, but we keep the logic consistent)
nEvents = size(eventSamples, 1);

if nEvents == 0
    % Return empty structure
    return;
end

%% ========================================================================
%  FINALIZE TIMES
%  ========================================================================
% Peak Finding
peakSample = zeros(nEvents, 1);
for iEvent = 1:nEvents
    idx = eventSamples(iEvent,1):eventSamples(iEvent,2);
    [~, peakRel] = min(rippSig.filt(idx));
    peakSample(iEvent) = eventSamples(iEvent,1) + peakRel - 1;
end

ripp.peakTime = timestamps(peakSample);
ripp.times = timestamps(eventSamples);

end
