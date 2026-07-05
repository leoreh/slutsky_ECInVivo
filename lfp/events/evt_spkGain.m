function gain = evt_spkGain(muTimes, evtTimes)
% EVT_SPKGAIN Per-event multi-unit spike gain (z-score vs matched controls).
%
%   gain = EVT_SPKGAIN(muTimes, evtTimes)
%
%   SUMMARY:
%       Scores how strongly pooled multi-unit activity (MUA) increases during
%       each event: the in-event MUA rate z-scored against the rate over matched
%       control intervals (evt_ctrlTimes). A positive gain means the event
%       carries above-baseline spiking, the signature of a genuine population
%       event. Returns NaN per event when no MUA is available (criterion skipped).
%
%   INPUTS:
%       muTimes  - (Cell) {1 x 1} Pooled MUA spike times [s].
%       evtTimes - (Mat)  [N x 2] Event start/end times [s].
%
%   OUTPUTS:
%       gain     - (Vec)  [N x 1] MUA gain z-score per event; NaN when no MUA.
%
%   DEPENDENCIES:
%       evt_ctrlTimes, times2rate.
%
%   HISTORY:
%       Created: 05 Jul 2026 (extracted from ripp_qa's spike-gain block).

N = size(evtTimes, 1);
gain = nan(N, 1);

if isempty(muTimes) || all(cellfun(@isempty, muTimes))
    return;                         % no MUA -> skip (all NaN)
end

% Matched control intervals (whole recording; mirrors the detection-time QA)
ctrlTimes = evt_ctrlTimes(evtTimes);

evtRates  = times2rate(muTimes, 'winCalc', evtTimes,  'binsize', Inf);
ctrlRates = times2rate(muTimes, 'winCalc', ctrlTimes, 'binsize', Inf);

muCtrl = mean(ctrlRates, 'all', 'omitnan');
sdCtrl = std(ctrlRates, [], 'all', 'omitnan');
if sdCtrl == 0, sdCtrl = 1; end

gain = (evtRates - muCtrl) ./ sdCtrl;
gain = gain(:);

end     % EOF
