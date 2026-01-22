function goodIdx = ripp_qa(ripp, rippSig, varargin)
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
%   goodIdx         (logical nEvents x 1) Combined accepted events.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'ripp', @isstruct);
addRequired(p, 'rippSig', @isstruct);
addParameter(p, 'spkGain', [], @isnumeric);
addParameter(p, 'emg', [], @isnumeric);        
addParameter(p, 'thrEmg', 2, @isnumeric);       
addParameter(p, 'thrGain', 0, @isnumeric);      
addParameter(p, 'fs', 1250, @isnumeric);        % Needed for sample conversion if not using ripp.times implicitly

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

badIdx = struct();
badIdx.emg = false(nEvents, 1);
badIdx.freq = false(nEvents, 1);
badIdx.gain = false(nEvents, 1);

% Convert times to samples for faster indexing (ensure bounds)
evtStart = round(ripp.times(:,1) * fs) + 1;
evtEnd = round(ripp.times(:,2) * fs) + 1;
evtStart = max(1, evtStart);
evtEnd = min(length(rippSig.lfp), evtEnd);

%% ========================================================================
%  EMG EXCLUSION
%  ========================================================================
% Logic: Calculate RMS of EMG during ripple event.
% Exclude if > percentile threshold (relative to detected events distribution).

emgRms = nan(nEvents, 1);
if ~isempty(emg)
    for iEvent = 1:nEvents
        idx = evtStart(iEvent):evtEnd(iEvent);
        emgRms(iEvent) = rms(emg(idx));
    end

    % Z-Score Method
    mu = mean(emgRms, 'omitnan');
    sigma = std(emgRms, 'omitnan');
    threshVal = mu + thrEmg * sigma;

    badIdx.emg = emgRms > threshVal;

    fprintf('QA: %d events flagged for High EMG (> %.2f std, Z-Score).\n', ...
        sum(badIdx.emg), thrEmg);
else
    warning('RIPP_QA: No EMG signal provided. Skipping EMG exclusion.');
end

%% ========================================================================
%  FREQUENCY EXCLUSION
%  ========================================================================
% Logic: Identify "High EMG" events (top 1% of EMG RMS).
% Calculate mean frequency of these high-noise events.
% Exclude ANY event with frequency higher than this noise-floor frequency.

% for iEvent = 1:nEvents
%     idx = evtStart(iEvent):evtEnd(iEvent);
%     qa.metrics.avgFreq(iEvent) = mean(rippSig.freq(idx), 'omitnan');
% end
%
% % Only proceed if we have valid EMG metrics to define "High EMG"
% if ~all(isnan(qa.metrics.emgRms))
%     highEmgThresh = prctile(qa.metrics.emgRms, 99);
%     highEmgIdx = qa.metrics.emgRms > highEmgThresh;
%
%     if any(highEmgIdx)
%         thrFreq = mean(qa.metrics.avgFreq(highEmgIdx));
%         badIdx.freq = qa.metrics.avgFreq > thrFreq;
%
%         fprintf('QA: %d events flagged for High Freq (> %.1f Hz).\n', ...
%             sum(badIdx.freq), thrFreq);
%     end
% end


%% ========================================================================
%  GAIN EXCLUSION
%  ========================================================================
% Logic: Check if spiking gain (modulation) is sufficient.

if ~isempty(spkGain)

    badIdx.gain = spkGain < thrGain;
    badIdx.gain = badIdx.gain(:);

    fprintf('QA: %d events flagged for Low Gain (< %.2f).\n', ...
        sum(badIdx.gain), thrGain);
end

%% ========================================================================
%  FINALIZE
%  ========================================================================
% Combine exclusions
goodIdx = ~(badIdx.emg | badIdx.freq | badIdx.gain);

end
