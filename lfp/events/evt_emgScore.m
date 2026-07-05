function emgZ = evt_emgScore(emg, evtWins, fs, varargin)
% EVT_EMGSCORE Per-event EMG z-score against a baseline distribution.
%
%   emgZ = EVT_EMGSCORE(emg, evtWins, fs, varargin)
%
%   SUMMARY:
%       Measures muscle (EMG) contamination per event: the mean rectified EMG
%       over each event window, z-scored against the rectified-EMG distribution
%       over a baseline period. The caller chooses the event windows (the event
%       duration for ripples, a fixed peak window for point events like EDs) and
%       the baseline (e.g. NREM bouts, or the whole recording). A high z-score
%       means the event coincides with movement and should be rejected.
%
%   INPUTS:
%       emg      - (Vec)  EMG trace (same fs as the detection signal). [] -> NaN.
%       evtWins  - (Mat)  [N x 2] Event windows [start end] to average EMG over [s].
%       fs       - (Num)  Sampling frequency [Hz].
%       varargin - Parameter/Value pairs:
%           'baselineTimes' - (Mat) [M x 2] Intervals defining the baseline
%                                   distribution [s]. Empty -> whole recording.
%
%   OUTPUTS:
%       emgZ     - (Vec)  [N x 1] EMG z-score per event; NaN when no EMG is given.
%
%   DEPENDENCIES:
%       None.
%
%   HISTORY:
%       Created: 05 Jul 2026 (unifies the ripp_qa EMG block + ed_reject_emg).

p = inputParser;
addRequired(p, 'emg', @(x) isempty(x) || isnumeric(x));
addRequired(p, 'evtWins', @isnumeric);
addRequired(p, 'fs', @isnumeric);
addParameter(p, 'baselineTimes', [], @isnumeric);
parse(p, emg, evtWins, fs, varargin{:});
baselineTimes = p.Results.baselineTimes;

N = size(evtWins, 1);
emgZ = nan(N, 1);
if isempty(emg), return; end        % no EMG -> skip (all NaN)

emgAbs = abs(emg(:));
nSamp = numel(emgAbs);

% Mean rectified EMG over each event window
eventEmg = nan(N, 1);
for iEv = 1:N
    b1 = max(1, round(evtWins(iEv, 1) * fs));
    b2 = min(nSamp, round(evtWins(iEv, 2) * fs));
    if b2 >= b1
        eventEmg(iEv) = mean(emgAbs(b1 : b2));
    end
end

% Baseline distribution (whole recording if no baselineTimes)
if isempty(baselineTimes)
    baseAmp = emgAbs;
else
    baseMask = false(nSamp, 1);
    for iB = 1:size(baselineTimes, 1)
        b1 = max(1, round(baselineTimes(iB, 1) * fs));
        b2 = min(nSamp, round(baselineTimes(iB, 2) * fs));
        if b2 >= b1, baseMask(b1 : b2) = true; end
    end
    baseAmp = emgAbs(baseMask);
end
if isempty(baseAmp), return; end

baseSd = std(baseAmp, 'omitnan');
if baseSd == 0, baseSd = 1; end
emgZ = (eventEmg - mean(baseAmp, 'omitnan')) / baseSd;

end     % EOF
