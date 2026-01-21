function qa = ripp_qa(ripp, rippSig, varargin)
% RIPP_QA Quality assurance for ripple detection.
%
% SUMMARY:
% Filters detected ripples based on EMG activity, Frequency content, and
% optional Spiking Gain.
%
% INPUT:
%   ripp            Structure with .times (nEvents x 2).
%   rippSig         Structure with .emg and .freq (and .z for consistency).
%   varargin        'rippSpks' (struct), 'thrEmg', 'thrGain'.
%
% OUTPUT:
%   qa              Structure with fields:
%                   .goodIdx    (logical nEvents x 1) Combined accepted events.
%                   .badIdx     (struct of logicals) Individual exclusion flags.
%                   .metrics    (struct) Calculated QA metrics per event.
%
% DEPENDENCIES: (None standard)

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'ripp', @isstruct);
addRequired(p, 'rippSig', @isstruct);
addParameter(p, 'spkGain', [], @isstruct);
addParameter(p, 'emg', [], @isnumeric); % Percentile threshold
addParameter(p, 'thrEmg', 90, @isnumeric); % Percentile threshold
addParameter(p, 'thrGain', 0, @isnumeric); % Z-score threshold
addParameter(p, 'fs', 1250, @isnumeric);    % Needed for sample conversion if not using ripp.times implicitly

parse(p, ripp, rippSig, varargin{:});
ripp        = p.Results.ripp;
rippSig     = p.Results.rippSig;
spkGain     = p.Results.spkGain;
emg         = p.Results.emg;
thrEmg      = p.Results.thrEmg;
thrGain     = p.Results.thrGain;
fs          = p.Results.fs;

%% ========================================================================
%  SETUP
%  ========================================================================
nEvents = size(ripp.times, 1);

qa = struct();
qa.goodIdx = true(nEvents, 1);
qa.badIdx = struct();
qa.badIdx.emg = false(nEvents, 1);
qa.badIdx.freq = false(nEvents, 1);
qa.badIdx.gain = false(nEvents, 1);

qa.metrics = struct();
qa.metrics.emgRms = nan(nEvents, 1);
qa.metrics.avgFreq = nan(nEvents, 1);
qa.metrics.rippGain = nan(nEvents, 1);

% Timestamps for signal indexing
timestamps = (0:length(rippSig.lfp)-1)' / fs;
% Convert times to samples for faster indexing (ensure bounds)
evtStart = round(ripp.times(:,1) * fs) + 1;
evtEnd = round(ripp.times(:,2) * fs) + 1;
evtStart = max(1, evtStart);
evtEnd = min(length(timestamps), evtEnd);

%% ========================================================================
%  EMG EXCLUSION
%  ========================================================================
% Logic: Calculate RMS of EMG during ripple event.
% Exclude if > percentile threshold (relative to detected events distribution).

if ~isempty(emg)
    for iEvent = 1:nEvents
        idx = evtStart(iEvent):evtEnd(iEvent);
        qa.metrics.emgRms(iEvent) = rms(emg(idx));
    end

    threshVal = prctile(qa.metrics.emgRms, thrEmg);
    qa.badIdx.emg = qa.metrics.emgRms > threshVal;

    fprintf('QA: %d events flagged for High EMG (> prctile %d).\n', ...
        sum(qa.badIdx.emg), thrEmg);
else
    warning('RIPP_QA: No EMG signal provided. Skipping EMG exclusion.');
end

%% ========================================================================
%  FREQUENCY EXCLUSION
%  ========================================================================
% Logic: Identify "High EMG" events (top 1% of EMG RMS).
% Calculate mean frequency of these high-noise events.
% Exclude ANY event with frequency higher than this noise-floor frequency.

if isfield(rippSig, 'freq') && ~isempty(rippSig.freq)
    for iEvent = 1:nEvents
        idx = evtStart(iEvent):evtEnd(iEvent);
        qa.metrics.avgFreq(iEvent) = mean(rippSig.freq(idx), 'omitnan');
    end

    % Only proceed if we have valid EMG metrics to define "High EMG"
    if ~all(isnan(qa.metrics.emgRms))
        highEmgThresh = prctile(qa.metrics.emgRms, 99);
        highEmgIdx = qa.metrics.emgRms > highEmgThresh;

        if any(highEmgIdx)
            thrFreq = mean(qa.metrics.avgFreq(highEmgIdx));
            qa.badIdx.freq = qa.metrics.avgFreq > thrFreq;

            fprintf('QA: %d events flagged for High Freq (> %.1f Hz).\n', ...
                sum(qa.badIdx.freq), thrFreq);
        end
    end
else
    warning('RIPP_QA: No Frequency signal provided. Skipping Frequency exclusion.');
end

%% ========================================================================
%  GAIN EXCLUSION
%  ========================================================================
% Logic: Check if spiking gain (modulation) is sufficient.

if ~isempty(spkGain)

    qa.badIdx.gain = spkGain < thrGain;
    qa.badIdx.gain = qa.badIdx.gain(:);

    fprintf('QA: %d events flagged for Low Gain (< %.2f).\n', ...
        sum(qa.badIdx.gain), thrGain);
end

%% ========================================================================
%  FINALIZE
%  ========================================================================
% Combine exclusions
qa.goodIdx = ~(qa.badIdx.emg | qa.badIdx.freq | qa.badIdx.gain);

end
