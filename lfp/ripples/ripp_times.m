function ripp = ripp_times(rippSig, fs, varargin)
% RIPP_TIMES Preliminary detection of hippocampal ripple events.
%
%   ripp = RIPP_TIMES(rippSig, fs, varargin)
%
%   SUMMARY:
%       Detects potential ripple events by thresholding the Z-scored signal.
%       Applies duration, power, and continuity constraints.
%
%   INPUTS:
%       rippSig     - (Struct) Signal structure from ripp_sigPrep.
%                              Must contain .z (Z-scored) and .filt fields.
%       fs          - (Num)    Sampling frequency [Hz].
%       varargin    - Parameter/Value pairs:
%           'thr'    - (Vec) Thresholds [Start, Peak, Cont, Max, Min_Cont].
%                            1: Start Trigger   (Z > 1.5)
%                            2: Peak Threshold  (Z > 2.5)
%                            3: Continuity Thr  (Z > 2.0) - must stay above this for...
%                            4: Max Threshold   (Z < 200) - artifact rejection.
%                            5: Min Cont Duration (ms) - time detection must obey Thr3.
%           'limDur' - (Vec) Duration limits [Min, Max, Inter, MinCont] (ms).
%                            1: Min Event Duration [15]
%                            2: Max Event Duration [300]
%                            3: Min Inter-Event Interval [20] (merge if closer)
%                            4: Min Continuity Duration [10]
%
%   OUTPUTS:
%       ripp        - (Struct) Detection results:
%           .times      - (N x 2) Start and End times [s].
%           .peakTime   - (N x 1) Peak times [s].
%           .info       - (Struct) Parameters used (.thr, .limDur, .fs).
%
%   DEPENDENCIES:
%       binary2bouts
%
%   HISTORY:
%       Updated: 23 Jan 2026

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

% Initialize
ripp = struct();
ripp.times = [];
ripp.peakTime = [];
ripp.info.fs = fs;
ripp.info.thr = thr;
ripp.info.limDur = limDur;

% Signal timestamps
tstamps = (0:length(rippSig.lfp)-1)' / fs;

% Duration limits in samples
limDur_Samples = round(limDur / 1000 * fs);

%% ========================================================================
%  DETECTION
%  ========================================================================
% 1. Initial Bout Detection
% Find all intervals where signal exceeds Start Threshold (Thr 1)
eventSamples = binary2bouts('vec', rippSig.z > thr(1), 'minDur', limDur_Samples(1),...
    'maxDur', limDur_Samples(2), 'interDur', limDur_Samples(3));
nEvents = size(eventSamples, 1);

idxDiscard = false(nEvents, 1);
peakPower = nan(nEvents, 1);
contPowDur = nan(nEvents, 1);
contPowAvg = nan(nEvents, 1);

% 2. Power and Continuity Validation
% Iterate through candidate events to enforce Peak and Continuity thresholds
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

% Check for Empty Results

% (Unused stats are cleared implicitly by not returning them, but we keep the logic consistent)
nEvents = size(eventSamples, 1);

if nEvents == 0
    % Return empty structure
    return;
end

%% ========================================================================
%  FINALIZE TIMES
%  ========================================================================

% 3. Find Exact Peak Times
% Locate the minimum (trough) of the filtered LFP within each event window
peakSample = zeros(nEvents, 1);
for iEvent = 1:nEvents
    idx = eventSamples(iEvent,1):eventSamples(iEvent,2);
    [~, peakRel] = min(rippSig.filt(idx));
    peakSample(iEvent) = eventSamples(iEvent,1) + peakRel - 1;
end

ripp.peakTime = tstamps(peakSample);
ripp.times = tstamps(eventSamples);

end
