function gain = evt_spkGain(muTimes, evtTimes, varargin)
% EVT_SPKGAIN Per-event multi-unit spike gain (z-score vs matched controls).
%
%   gain = EVT_SPKGAIN(muTimes, evtTimes, varargin)
%
%   SUMMARY:
%       Scores how strongly pooled multi-unit activity (MUA) increases during
%       each event: the in-event MUA rate z-scored against the rate over matched
%       control intervals (evt_ctrlTimes). A positive gain means the event
%       carries above-baseline spiking, the signature of a genuine population
%       event. Returns NaN per event when no MUA is available (criterion skipped).
%
%       Each event window is widened by a margin before the rate is counted:
%       ripple-associated spiking leads the detected start and lags the detected
%       end, because the band-pass envelope that sets the edges clips these
%       spiking shoulders. A small margin recovers the event's own recruitment;
%       controls are matched to the widened durations, so the z-score reference
%       stays unbiased. Empirically (lh100 NREM) the ripple-specific shoulder
%       sits within ~5 ms of each edge - beyond that the elevated rate is the
%       broad sharp-wave / clustering floor, which should not count per event.
%
%   INPUTS:
%       muTimes  - (Cell) {1 x 1} Pooled MUA spike times [s].
%       evtTimes - (Mat)  [N x 2] Event start/end times [s].
%       varargin - Parameter/Value:
%           'margin' - (Num) Symmetric widening [s] applied to every event
%                            before the rate and the matched controls.
%                            (Default: 0.005).
%
%   OUTPUTS:
%       gain     - (Vec)  [N x 1] MUA gain z-score per event; NaN when no MUA.
%
%   DEPENDENCIES:
%       evt_ctrlTimes, times2rate.
%
%   HISTORY:
%       Created: 05 Jul 2026 (extracted from ripp_qa's spike-gain block).
%       Updated: 260719 (add the 'margin' shoulder; controls match the widened
%                windows).

p = inputParser;
addRequired(p, 'muTimes', @iscell);
addRequired(p, 'evtTimes', @isnumeric);
addParameter(p, 'margin', 0.005, @(x) isnumeric(x) && isscalar(x) && x >= 0);
parse(p, muTimes, evtTimes, varargin{:});
margin = p.Results.margin;

N = size(evtTimes, 1);
gain = nan(N, 1);

if isempty(muTimes) || all(cellfun(@isempty, muTimes))
    return;                         % no MUA -> skip (all NaN)
end

% Widen each event by the margin (spiking leads the start and lags the end).
% Controls are matched to the same widened durations; evt_ctrlTimes merges any
% overlaps the widening creates and crops to the recording.
evtWide = [evtTimes(:, 1) - margin, evtTimes(:, 2) + margin];

% Matched control intervals (whole recording; mirrors the detection-time QA)
ctrlTimes = evt_ctrlTimes(evtWide);

evtRates  = times2rate(muTimes, 'winCalc', evtWide,   'binsize', Inf);
ctrlRates = times2rate(muTimes, 'winCalc', ctrlTimes, 'binsize', Inf);

muCtrl = mean(ctrlRates, 'all', 'omitnan');
sdCtrl = std(ctrlRates, [], 'all', 'omitnan');
if sdCtrl == 0
    sdCtrl = 1;
end

gain = (evtRates - muCtrl) ./ sdCtrl;
gain = gain(:);

end     % EOF
