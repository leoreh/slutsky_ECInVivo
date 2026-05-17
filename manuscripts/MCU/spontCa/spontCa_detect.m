function evTbl = spontCa_detect(trace, fs, varargin)
% SPONTCA_DETECT Detects Ca transients in a single dF/F trace.
%
%   evTbl = SPONTCA_DETECT(TRACE, FS, ...) detects events by their RISE:
%   a positive crossing of the smoothed derivative on the trace. Each event
%   is then validated by an amplitude check (peak >= minAmp) and by the
%   presence of at least one significantly negative derivative sample
%   after the peak. Stops walk forward until the derivative flattens
%   (three consecutive |d| < stopFlat samples). Amplitudes and integrals
%   are reported directly from the input trace.
%
%   Why derivative-based: at fs=3 Hz the rise of a Ca event is ~1 sample,
%   producing one large positive d/dt; the slow decay is many samples of
%   small negative d/dt. Plateau noise is symmetric and small in d/dt.
%   This is the signature that separates real events from plateau noise,
%   which amplitude thresholds on the raw trace cannot.
%
%   Note on baseline: this function does NOT subtract any rolling baseline
%   from the trace. The caller is expected to feed in a trace where slow
%   drift / plateaus have already been removed (e.g., the Pass B step in
%   spontCa_detectWrapper). amp / int are read directly off the trace.
%
%   INPUTS:
%       trace - (1 x nT) dF/F signal
%       fs    - (scalar) sampling rate (Hz)
%
%   OPTIONAL (Name-Value, user-facing):
%       'minAmp' - (num) min peak amplitude on trace (dF/F)              {0.05}
%       'minIEI' - (num) min peak-to-peak distance (s)                   {1.0}
%       'kNoise' - (num) rise-threshold multiplier of derivative noise   {3.5}
%       'minDur' - (num) min decay duration (stop - peak, s)             {0.4}
%
%   OUTPUT:
%       evTbl - table with one row per event, columns:
%         start (s)   peak time (start == peak by convention)
%         stop  (s)   decay-end time (where derivative flattens)
%         amp   (dF/F) peak amplitude on the input trace
%         dur   (s)   stop - start (decay length, not event span)
%         int   (dF/F * s) integral of trace over [peak, stop]
%
%   See also: SPONTCA_FINALIZE, SPONTCA_LOAD, SPONTCA_MANCUR

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'trace', @(x) isnumeric(x) && isvector(x));
addRequired(p, 'fs',    @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'minAmp', 0.05, @isnumeric);
addParameter(p, 'minIEI', 1.0,  @isnumeric);
addParameter(p, 'kNoise', 3.5,  @isnumeric);
addParameter(p, 'minDur', 0.4,  @isnumeric);
parse(p, trace, fs, varargin{:});
P = p.Results;

% Buried constants (do not expose; tune in source if you must).
pkLookAhead  = 3;      % samples to search for peak after derivative crossing
minDecaySmp  = 1;      % samples of strongly negative d/dt required after peak
kStop        = 1.0;    % multiplier on sigma_dt for "strongly negative" test
stopFlat     = 0.005;  % |d| below this counts as flat (dF/F per sample)
stopK        = 3;      % consecutive flat samples to terminate the walk-forward

trace = trace(:)';
nT = length(trace);
dt = 1 / fs;


%% ========================================================================
%  EARLY RETURN
%  ========================================================================

evTbl = emptyEvTbl();
if all(isnan(trace)) || nT < 4
    return;
end


%% ========================================================================
%  DERIVATIVE WITH LIGHT SMOOTHING
%  ========================================================================
% 3-sample boxcar on the trace stabilises the per-sample diff against
% high-frequency noise without distorting the peak time. Compute the
% derivative on the smoothed signal; amplitudes / integrals still use raw.

traceSm = movmean(trace, 3, 'omitnan');
d = [0, diff(traceSm)];   % per-sample dF/F change; same length as trace


%% ========================================================================
%  PER-CELL DERIVATIVE NOISE (MAD)
%  ========================================================================
% MAD on the derivative (not the trace) is robust to plateaus: slow
% baselines contribute near-zero d, real events contribute large positive
% spikes that the median rejects. Sigma_dt measures noise.

dValid = d(~isnan(d));
if isempty(dValid)
    return;
end
sigmaD = 1.4826 * median(abs(dValid - median(dValid)));
if sigmaD <= 0
    sigmaD = eps;
end
thrRise = P.kNoise * sigmaD;
thrNeg  = kStop * sigmaD;


%% ========================================================================
%  CANDIDATE RISES (leading edges of d > thrRise runs)
%  ========================================================================

isRise = d > thrRise;
isRise(isnan(d)) = false;
leadEdges = find(isRise & ~[false, isRise(1:end-1)]);

if isempty(leadEdges)
    return;
end


%% ========================================================================
%  LOCATE PEAK (walk forward from rise onset up to pkLookAhead samples)
%  ========================================================================

nCand   = length(leadEdges);
peakIdx = zeros(1, nCand);
for k = 1:nCand
    i0 = leadEdges(k);
    iLast = min(nT, i0 + pkLookAhead);
    bestIdx = i0;
    bestVal = trace(i0);
    for j = (i0 + 1):iLast
        if isnan(trace(j))
            break;
        end
        if trace(j) >= bestVal
            bestVal = trace(j);
            bestIdx = j;
        else
            break;
        end
    end
    peakIdx(k) = bestIdx;
end


%% ========================================================================
%  AMPLITUDE GATE
%  ========================================================================

ampLocal = trace(peakIdx);
keep     = ampLocal >= P.minAmp;
peakIdx  = peakIdx(keep);
ampLocal = ampLocal(keep);
nCand    = length(peakIdx);
if nCand == 0
    return;
end


%% ========================================================================
%  DECAY VALIDATION (at least one strongly negative sample after peak)
%  ========================================================================
% Rejects step-up-and-stay candidates: a real event has a falling phase;
% a plateau onset that just steps up to a new level does not, unless the
% sustained activity contains its own decay phases.

keep = false(1, nCand);
for k = 1:nCand
    rL = peakIdx(k) + 1;
    rR = min(nT, peakIdx(k) + pkLookAhead + 3);
    if rR < rL
        continue;
    end
    seg = d(rL:rR);
    if any(seg < -thrNeg)
        keep(k) = true;
    end
end
peakIdx  = peakIdx(keep);
ampLocal = ampLocal(keep);
nCand    = length(peakIdx);
if nCand == 0
    return;
end


%% ========================================================================
%  GREEDY IEI SUPPRESSION
%  ========================================================================
% Sort peaks by amplitude descending, keep highest, drop any within
% minIEI of a kept peak.

minIEISmp = max(1, round(P.minIEI * fs));
if nCand > 1
    [~, order] = sort(ampLocal, 'descend');
    keepFlag = false(1, nCand);
    keepFlag(order(1)) = true;
    for k = 2:nCand
        candPos = peakIdx(order(k));
        keptPos = peakIdx(keepFlag);
        if all(abs(candPos - keptPos) >= minIEISmp)
            keepFlag(order(k)) = true;
        end
    end
    peakIdx  = peakIdx(keepFlag);
    ampLocal = ampLocal(keepFlag);
    [peakIdx, sortByTime] = sort(peakIdx);
    ampLocal = ampLocal(sortByTime);
end
nEv = length(peakIdx);


%% ========================================================================
%  WALK FORWARD TO STOP (derivative flattens)
%  ========================================================================
% Bounded by the next surviving peak so events do not share samples.
% Increment a counter when |d| < stopFlat; reset when above. Stop when
% counter reaches stopK consecutive flat samples.

stopIdx = zeros(1, nEv);
for iE = 1:nEv
    pk = peakIdx(iE);
    if iE < nEv
        rightBound = peakIdx(iE + 1) - 1;
    else
        rightBound = nT;
    end
    flatCount = 0;
    e = pk;
    while e < rightBound
        e1 = e + 1;
        if isnan(d(e1))
            break;
        end
        if abs(d(e1)) < stopFlat
            flatCount = flatCount + 1;
        else
            flatCount = 0;
        end
        e = e1;
        if flatCount >= stopK
            break;
        end
    end
    stopIdx(iE) = e;
end


%% ========================================================================
%  MIN-DURATION FILTER
%  ========================================================================

durSmp = stopIdx - peakIdx + 1;
keep   = durSmp >= max(2, round(P.minDur * fs));
peakIdx  = peakIdx(keep);
stopIdx  = stopIdx(keep);
ampLocal = ampLocal(keep);
nEv      = length(peakIdx);


%% ========================================================================
%  INTEGRAL ON TRACE
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

starts = (peakIdx(:) - 1) * dt;     % start == peak time
stops  = (stopIdx(:) - 1) * dt;
amps   = ampLocal(:);                % amplitude at peak (on input trace)
durs   = stops - starts;
evTbl  = table(starts, stops, amps, durs, intg, ...
    'VariableNames', {'start', 'stop', 'amp', 'dur', 'int'});

end     % EOF


function tbl = emptyEvTbl()
% Canonical empty events table (zero rows, correct column types).
tbl = table( ...
    zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
    'VariableNames', {'start', 'stop', 'amp', 'dur', 'int'});
end
