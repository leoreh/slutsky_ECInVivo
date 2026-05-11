function ev = spontCa_detect(trace, fs, varargin)
% SPONTCA_DETECT Detects Ca transients in a single dF/F trace.
%
%   ev = SPONTCA_DETECT(TRACE, FS, ...) treats each event as a local
%   maximum followed by a decay. The event's START is the peak sample
%   itself; the STOP is the sample where the trace returns to thrBsl
%   while walking forward from the peak (or the next peak, whichever
%   comes first). No walk-back, no foot, no rise concept. Duration is
%   STOP - START, i.e. the decay length.
%
%   The input is assumed to be dF/F already. No rolling baseline.
%
%   PIPELINE:
%       1. Find local maxima ("rise stops"): samples i where
%          trace(i) > trace(i-1) AND trace(i) >= trace(i+1). Catches
%          impulse peaks and plateau onsets in one rule.
%       2. Keep only peaks with trace(i) > minAmp.
%       3. Greedy max-suppression in minIEI windows: keep the highest
%          peak per window, drop the rest.
%       4. Walk forward from each peak until trace <= thrBsl OR until
%          the next surviving peak. The walk-forward end is STOP.
%       5. Drop events whose stop - peak duration is below minDur.
%
%   INPUTS:
%       trace - (1 x nT) dF/F signal
%       fs    - (scalar) sampling rate (Hz)
%
%   OPTIONAL (Name-Value):
%       'minAmp' - (num) min peak amplitude (dF/F)              {0.08}
%       'minDur' - (num) min decay duration (stop - peak, s)    {0.4}
%       'minIEI' - (num) min peak-to-peak distance (s)          {1.0}
%       'thrBsl' - (num) absolute return threshold (dF/F)       {0.02}
%
%   OUTPUT:
%       ev struct with fields:
%         .start  (n x 1)  peak times (s)  -- this IS the peak
%         .stop   (n x 1)  decay-end times (s)
%         .amp    (n x 1)  peak amplitude (dF/F)
%         .dur    (n x 1)  stop - start (s)
%         .int    (n x 1)  integral over [peak, stop] (dF/F * s)
%
%   See also: SPONTCA_EVENTS, SPONTCA_COUPLE, SPONTCA_LOAD

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'trace', @(x) isnumeric(x) && isvector(x));
addRequired(p, 'fs',    @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'minAmp', 0.08, @isnumeric);
addParameter(p, 'minDur', 0.4,  @isnumeric);
addParameter(p, 'minIEI', 1.0,  @isnumeric);
addParameter(p, 'thrBsl', 0.02, @isnumeric);
parse(p, trace, fs, varargin{:});
P = p.Results;

trace = trace(:)';
nT = length(trace);
dt = 1 / fs;


%% ========================================================================
%  LOCAL MAXIMA (rise stops above minAmp)
%  ========================================================================

isUp   = false(1, nT);
isFlat = false(1, nT);
isUp(2:end)     = trace(2:end)   > trace(1:end-1);
isFlat(1:end-1) = trace(1:end-1) >= trace(2:end);
isRiseStop = isUp & isFlat & (trace > P.minAmp);
isRiseStop(isnan(trace)) = false;
peakIdx = find(isRiseStop);
pkVals  = trace(peakIdx);


%% ========================================================================
%  GREEDY MAX-SUPPRESSION (min peak distance)
%  ========================================================================

minIEISmp = max(1, round(P.minIEI * fs));
nEv = length(peakIdx);
if nEv > 1
    [~, order] = sort(pkVals, 'descend');
    keep = false(1, nEv);
    keep(order(1)) = true;
    for k = 2:nEv
        cand    = peakIdx(order(k));
        keptPos = peakIdx(keep);
        if all(abs(cand - keptPos) >= minIEISmp)
            keep(order(k)) = true;
        end
    end
    peakIdx = peakIdx(keep);
    [peakIdx, sortByTime] = sort(peakIdx);
    pkVals = pkVals(keep);
    pkVals = pkVals(sortByTime);
    nEv = length(peakIdx);
end


%% ========================================================================
%  WALK FORWARD FROM PEAK (peak -> stop)
%  ========================================================================
% Walk forward while the trace is above thrBsl. Bounded by the next
% surviving peak so adjacent events do not share samples.

stopIdx = zeros(1, nEv);
for iE = 1:nEv
    pk = peakIdx(iE);
    if iE < nEv
        rightBound = peakIdx(iE + 1) - 1;
    else
        rightBound = nT;
    end
    e = pk;
    while e < rightBound
        nxt = trace(e + 1);
        if isnan(nxt) || nxt <= P.thrBsl
            break;
        end
        e = e + 1;
    end
    stopIdx(iE) = e;
end


%% ========================================================================
%  MIN-DURATION FILTER
%  ========================================================================

durSmp  = stopIdx - peakIdx + 1;
keep    = durSmp >= max(2, round(P.minDur * fs));
peakIdx = peakIdx(keep);
stopIdx = stopIdx(keep);
pkVals  = pkVals(keep);
nEv     = length(peakIdx);


%% ========================================================================
%  INTEGRAL OVER [peak, stop]
%  ========================================================================

intg = nan(nEv, 1);
for iE = 1:nEv
    seg = trace(peakIdx(iE):stopIdx(iE));
    seg(isnan(seg)) = 0;
    intg(iE) = trapz(seg) * dt;
end


%% ========================================================================
%  ASSEMBLE OUTPUT
%  ========================================================================
% start == peak. dur = stop - start = decay length.

ev = struct();
ev.start = (peakIdx(:) - 1) * dt;
ev.stop  = (stopIdx(:) - 1) * dt;
ev.amp   = pkVals(:);
ev.dur   = ev.stop - ev.start;
ev.int   = intg;

end     % EOF
