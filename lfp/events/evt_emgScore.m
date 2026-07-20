function emgZ = evt_emgScore(emg, evtWins, fs, varargin)
% EVT_EMGSCORE Per-event EMG z-score against a baseline distribution.
%
%   emgZ = EVT_EMGSCORE(emg, evtWins, fs, varargin)
%
%   SUMMARY:
%       Measures muscle (EMG) contamination per event: the mean rectified EMG
%       over each event window, standardized against the rectified-EMG
%       distribution over a baseline period. The caller chooses the event
%       windows (the event duration for ripples, a fixed peak window for point
%       events like EDs) and the baseline (e.g. NREM bouts, or the whole
%       recording). A high score means the event coincides with movement and
%       should be rejected.
%
%       The comparison is made on the LOG of the amplitude, with a ROBUST
%       (median / MAD) baseline. EMG amplitude is multiplicative and heavy
%       tailed, so a linear mean/SD tracks the largest movement transients
%       rather than the resting bulk and the scale drifts between sessions -
%       which is what makes one threshold fail to transfer across mice. See the
%       note at the standardization step for the measured effect.
%
%       The baseline statistics depend only on the samples, not on the windows
%       being scored. Scoring a regular grid of windows against the same
%       baseline therefore yields the identical scale - which is how the viewer
%       draws a continuous trace comparable to the per-event score (the
%       'emgScore' transform in var_fetch).
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
%       emgZ     - (Vec)  [N x 1] Robust z per event (median / MAD of the
%                         baseline); NaN when no EMG is given.
%
%   DEPENDENCIES:
%       None.
%
%   HISTORY:
%       Created: 05 Jul 2026 (unifies the ripp_qa EMG block + ed_reject_emg).
%       260720 standardization moved to log amplitude with a median / MAD
%              baseline (was linear mean / SD). Chosen over a linear robust z
%              and over a window-matched baseline by measuring, on the MCU
%              cohort, the spread of the per-mouse WAKE/NREM boundary:
%              IQR/median 0.59 linear-robust, 0.37 window-matched, 0.34 log.
%              The optimal cut is now ~1.0 (0.79-1.39 over 15 mice). A stored
%              ripp.emg / ed.emgZ from an earlier run is on the old scale and is
%              NOT comparable - re-detect before reusing a threshold.

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
eventAmp = nan(N, 1);
for iEv = 1:N
    b1 = max(1, round(evtWins(iEv, 1) * fs));
    b2 = min(nSamp, round(evtWins(iEv, 2) * fs));
    if b2 >= b1
        eventAmp(iEv) = mean(emgAbs(b1 : b2));
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

% Compare on a LOG scale, robustly. EMG amplitude is multiplicative - gain,
% impedance and movement all scale it - so the rectified signal is strongly
% right skewed and mean/SD track the largest transients instead of the resting
% bulk. Taking the log makes that bulk near symmetric (the same reason
% AccuSleep stores log-RMS), and median/MAD then estimate its centre and width
% without the tail. Measured over the MCU cohort, this halves the spread of the
% WAKE/NREM decision boundary between mice (IQR/median 0.59 -> 0.34), which is
% what lets a single threshold transfer. Note the baseline statistics depend
% only on the samples, not on the event windows, so scoring any other set of
% windows against the same baseline lands on exactly this scale.
maxSamp = 2e6;                          % a baseline can reach 1e8 samples;
if numel(baseAmp) > maxSamp             % median/MAD are stable under a stride,
    baseAmp = baseAmp(1 : ceil(numel(baseAmp) / maxSamp) : end);   % and a fixed
end                                     % stride keeps the result reproducible

baseLog = log(baseAmp(baseAmp > 0));    % a zero sample carries no amplitude
if isempty(baseLog), return; end

baseLoc = median(baseLog);
baseScl = 1.4826 * median(abs(baseLog - baseLoc));
if ~isfinite(baseScl) || baseScl == 0
    baseScl = std(baseLog);             % degenerate bulk (flat / clipped)
end
if ~isfinite(baseScl) || baseScl == 0, baseScl = 1; end

evtLog = log(eventAmp);                 % a dead window -> NaN (metric absent),
evtLog(~isfinite(evtLog)) = NaN;        % never -Inf
emgZ = (evtLog - baseLoc) / baseScl;

end     % EOF
