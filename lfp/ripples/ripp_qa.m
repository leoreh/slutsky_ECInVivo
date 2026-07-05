function [idxQA, met] = ripp_qa(rippTimes, nremTimes, muTimes, emg, rippState, fs)
% RIPP_QA Quality-assurance filter for candidate ripple events.
%
%   [idxQA, met] = RIPP_QA(rippTimes, nremTimes, muTimes, emg, rippState, fs)
%
%   SUMMARY:
%       Scores each candidate ripple on three criteria and returns the
%       intersection as the pass mask. Each criterion is skipped (all events
%       pass) when its supporting data is absent, so unlabelled / spikeless /
%       EMG-less sessions are handled gracefully:
%           - State: keep events in QWAKE / LSLEEP / NREM.
%           - Spike gain: MUA rate in-event z-scored against control > 0.
%           - EMG: in-event EMG z-scored against the NREM baseline < 2.
%
%   INPUTS:
%       rippTimes - (Mat)  [N x 2] Event start/end times (s).
%       nremTimes - (Mat)  [M x 2] NREM bout times for the EMG baseline (s).
%       muTimes   - (Cell) {sorted_vec} Multi-unit spike times (s).
%       emg       - (Vec)  EMG signal.
%       rippState - (Cat)  [N x 1] State labels.
%       fs        - (Num)  Sampling rate (Hz).
%
%   OUTPUTS:
%       idxQA     - (Struct) Logical masks: .state / .gain / .emg / .good.
%       met       - (Struct) Per-event metrics: .emg (z), .gain (z).
%
%   DEPENDENCIES:
%       evt_ctrlTimes, times2rate.
%
%   HISTORY:
%       Updated: 05 Jul 2026 (promoted from a ripp_wrapper subfunction to a file).

% EMG Analysis
% -------------------------------------------------------------------------
% Skipped (all events pass) when no EMG or no NREM baseline is available.
nRipp = size(rippTimes, 1);
rippEmgZ = nan(nRipp, 1);
idxEmg = true(nRipp, 1);

if ~isempty(emg) && ~isempty(nremTimes)
    emgAbs = abs(emg);
    rippEmg = nan(nRipp, 1);
    for iRipp = 1:nRipp
        bStart = max(1, round(rippTimes(iRipp,1) * fs));
        bEnd = min(length(emg), round(rippTimes(iRipp,2) * fs));
        rippEmg(iRipp) = mean(emgAbs(bStart : bEnd));
    end

    % NREM Baseline
    nremMask = false(size(emg));
    for iBout = 1:size(nremTimes, 1)
        bStart = max(1, round(nremTimes(iBout,1) * fs));
        bEnd = min(length(emg), round(nremTimes(iBout,2) * fs));
        nremMask(bStart : bEnd) = true;
    end
    nremAmp = emgAbs(nremMask);

    % Z-Score against the NREM baseline
    if ~isempty(nremAmp)
        rippEmgZ = (rippEmg - mean(nremAmp, 'omitnan')) / std(nremAmp, 'omitnan');
        idxEmg = rippEmgZ < 2;
    end
end

% Spike Gain
% -------------------------------------------------------------------------
% Skipped (all events pass) when no multi-unit spikes are available.
spkGain = nan(nRipp, 1);
idxGain = true(nRipp, 1);

muEmpty = isempty(muTimes) || all(cellfun(@isempty, muTimes));
if ~muEmpty
    % Generate control intervals locally
    ctrlTimes = evt_ctrlTimes(rippTimes);

    % Calculate Rates (muTimes is expected to be {vector})
    rippRates = times2rate(muTimes, 'winCalc', rippTimes, 'binsize', Inf);
    ctrlRates = times2rate(muTimes, 'winCalc', ctrlTimes, 'binsize', Inf);

    % Gain = (Ripple Rate - Mean Control) / Std Control
    muCtrl = mean(ctrlRates, 'all', 'omitnan');
    sdCtrl = std(ctrlRates, [], 'all', 'omitnan');
    if sdCtrl == 0, sdCtrl = 1; end

    spkGain = (rippRates - muCtrl) ./ sdCtrl;
    % Ensure gain is column vector for consistency
    if size(spkGain, 2) > size(spkGain, 1), spkGain = spkGain'; end

    idxGain = spkGain > 0;
end

% States
% -------------------------------------------------------------------------
% Keep events in valid vigilance states. When no states are assigned at all
% (all undefined), the criterion is skipped so unlabelled sessions pass.
validStates = {'QWAKE', 'LSLEEP', 'NREM'};
if iscategorical(rippState) && all(isundefined(rippState))
    idxState = true(nRipp, 1);
else
    idxState = ismember(rippState, validStates);
    idxState = idxState(:);
end

% Combine
% -------------------------------------------------------------------------
idxGood = idxState & idxGain & idxEmg;

% Output
idxQA.state = idxState;
idxQA.gain = idxGain;
idxQA.emg = idxEmg;
idxQA.good = idxGood;

met.emg = rippEmgZ;
met.gain = spkGain;

end     % EOF
