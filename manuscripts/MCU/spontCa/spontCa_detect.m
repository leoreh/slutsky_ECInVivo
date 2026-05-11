function ev = spontCa_detect(trace, fs, varargin)
% SPONTCA_DETECT Detects Ca transients in a single dF/F trace.
%
%   ev = SPONTCA_DETECT(TRACE, FS, ...) detects events as rising flanks
%   in the trace. Each rising flank's last sample is the peak; the first
%   sample is the foot. No FINDPEAKS, no prominence, no trough-bounded
%   splitting.
%
%   The input is assumed to be dF/F already (centered near zero outside
%   events). No rolling baseline is subtracted. Noise is the robust SD of
%   the differenced signal (Allan-Pettersson) and feeds the height floor.
%
%   PIPELINE:
%       1. Identify rise stops: samples i where trace(i) > trace(i-1) and
%          trace(i) >= trace(i+1). This catches both impulse peaks (the
%          local max before a decay) and plateau onsets (the first sample
%          of a flat top after a rise).
%       2. Filter rise stops by absolute height: trace(i) > max(kThr *
%          noise, minAmp).
%       3. Enforce min peak distance (minIEI) via greedy max suppression.
%       4. For each surviving rise stop, walk back along the rising flank
%          while trace(s-1) < trace(s). This gives the foot of the rise -
%          the most recent sample whose value is strictly below every
%          subsequent sample up to the peak. The rise amplitude is
%          peak - trace(foot).
%       5. Drop events whose rise amplitude is below minRise. This is the
%          true rise of THIS event, not standard FINDPEAKS prominence
%          (which uses the higher flanking trough and breaks for
%          asymmetric events whose decays merge into plateaus).
%       6. Determine each event's stop: walk forward from the peak while
%          trace > thrBsl AND we have not reached the next event's foot.
%          This lets long decays extend through the plateau without
%          inventing an artificial split at the inter-peak trough.
%       7. Drop events shorter than minDur on the extended span.
%
%   INPUTS:
%       trace   - (1 x nT) dF/F signal
%       fs      - (scalar) sampling rate (Hz)
%
%   OPTIONAL (Name-Value):
%       'kThr'        - (num) height in noise SDs                 {3}
%       'minAmp'      - (num) absolute floor on peak amp (dF/F)   {0.05}
%       'minRise'     - (num) min rise above foot (dF/F)          {0.05}
%       'minRiseBnd'  - (num) min rise for an event to bound a
%                             neighbour's walk-forward and to absorb
%                             smaller overlapping events (dF/F)   {0.10}
%       'minDur'      - (num) min event duration on extended span (s) {0.4}
%       'minIEI'      - (num) min peak-to-peak distance (s)       {0.4}
%       'thrBsl'      - (num) absolute return threshold (dF/F)    {0.02}
%
%   OUTPUT:
%       ev struct with fields:
%         .start  (n x 1)  foot times (s)
%         .stop   (n x 1)  tail times (s)
%         .peak   (n x 1)  peak times (s)
%         .amp    (n x 1)  peak amplitude (dF/F)
%         .dur    (n x 1)  duration foot-to-tail (s)
%         .int    (n x 1)  integral over the extended span (dF/F * s)
%         .rise   (n x 1)  peak - trace(foot) (dF/F)
%         .noise  (scalar) robust noise std
%
%   See also: SPONTCA_EVENTS, SPONTCA_COUPLE, SPONTCA_LOAD

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'trace', @(x) isnumeric(x) && isvector(x));
addRequired(p, 'fs',    @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'kThr',    3,    @isnumeric);
addParameter(p, 'minAmp',  0.05, @isnumeric);
addParameter(p, 'minRise',    0.05, @isnumeric);
addParameter(p, 'minRiseBnd', 0.10, @isnumeric);
addParameter(p, 'minDur',     0.4,  @isnumeric);
addParameter(p, 'minIEI',     0.4,  @isnumeric);
addParameter(p, 'thrBsl',     0.02, @isnumeric);
parse(p, trace, fs, varargin{:});
P = p.Results;

trace = trace(:)';
nT = length(trace);
dt = 1 / fs;


%% ========================================================================
%  NOISE
%  ========================================================================

dx = diff(trace);
dx = dx(~isnan(dx));
if isempty(dx) || mad(dx, 1) == 0
    sNoise = max(eps, std(trace, 'omitnan'));
else
    sNoise = 1.4826 * mad(dx, 1) / sqrt(2);
end


%% ========================================================================
%  RISE STOPS (peak candidates)
%  ========================================================================
% trace(i) > trace(i-1) marks "the rise hit sample i".
% trace(i) >= trace(i+1) marks "sample i is not below the next sample".
% Together: the rise ended at sample i, either because i is an impulse
% peak (next sample is lower) or because i is the first sample of a
% plateau (next sample is equal).

minHeight = max(P.kThr * sNoise, P.minAmp);

isUp   = false(1, nT);
isFlat = false(1, nT);
isUp(2:end)   = trace(2:end)   > trace(1:end-1);
isFlat(1:end-1) = trace(1:end-1) >= trace(2:end);
isRiseStop = isUp & isFlat & (trace > minHeight);
isRiseStop(isnan(trace)) = false;
peakIdx = find(isRiseStop);
pkVals  = trace(peakIdx);


%% ========================================================================
%  MIN PEAK DISTANCE (greedy max suppression)
%  ========================================================================

minIEISmp = max(1, round(P.minIEI * fs));
nEv = length(peakIdx);
if nEv > 1
    [~, order] = sort(pkVals, 'descend');
    keep = false(1, nEv);
    keep(order(1)) = true;
    for k = 2:nEv
        cand = peakIdx(order(k));
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
%  FOOT (walk back along the rising flank)
%  ========================================================================

footIdx = zeros(1, nEv);
for iE = 1:nEv
    s = peakIdx(iE);
    while s > 1
        prev = trace(s - 1);
        if isnan(prev) || prev >= trace(s)
            break;
        end
        s = s - 1;
    end
    footIdx(iE) = s;
end


%% ========================================================================
%  RISE FILTER
%  ========================================================================

riseAmp = pkVals - trace(footIdx);
keep    = riseAmp >= P.minRise;
peakIdx = peakIdx(keep);
footIdx = footIdx(keep);
pkVals  = pkVals(keep);
riseAmp = riseAmp(keep);
nEv     = length(peakIdx);


%% ========================================================================
%  TAIL (walk forward from peak)
%  ========================================================================
% Stop when trace drops to or below thrBsl, OR we reach the foot of the
% next SIGNIFICANT event (rise >= minRiseBnd). Spurious "events" on a
% decay tail have rise just above minRise but well below minRiseBnd, so
% they do not bound the preceding real event's walk-forward.

isSig = riseAmp >= P.minRiseBnd;

stopIdx = zeros(1, nEv);
for iE = 1:nEv
    pk = peakIdx(iE);
    rightBound = nT;
    for jE = (iE + 1):nEv
        if isSig(jE)
            rightBound = footIdx(jE) - 1;
            break;
        end
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
%  ABSORB SPURIOUS EVENTS INSIDE SIGNIFICANT-EVENT SPANS
%  ========================================================================
% A non-significant event whose peak lies inside a significant event's
% extended span is absorbed (dropped). The significant event keeps its
% own properties; the spurious one is removed entirely.

sigIdx = find(isSig);
absorb = false(1, nEv);
for iE = 1:nEv
    if isSig(iE), continue; end
    for k = 1:length(sigIdx)
        jE = sigIdx(k);
        if peakIdx(iE) >= footIdx(jE) && peakIdx(iE) <= stopIdx(jE)
            absorb(iE) = true;
            break;
        end
    end
end
keep    = ~absorb;
peakIdx = peakIdx(keep);
footIdx = footIdx(keep);
stopIdx = stopIdx(keep);
pkVals  = pkVals(keep);
riseAmp = riseAmp(keep);
nEv     = length(peakIdx);


%% ========================================================================
%  MIN DUR
%  ========================================================================

durSmp = stopIdx - footIdx + 1;
keep   = durSmp >= max(2, round(P.minDur * fs));
peakIdx = peakIdx(keep);
footIdx = footIdx(keep);
stopIdx = stopIdx(keep);
pkVals  = pkVals(keep);
riseAmp = riseAmp(keep);
nEv     = length(peakIdx);


%% ========================================================================
%  INTEGRAL OVER EXTENDED SPAN
%  ========================================================================

intg = nan(nEv, 1);
for iE = 1:nEv
    seg = trace(footIdx(iE):stopIdx(iE));
    seg(isnan(seg)) = 0;
    intg(iE) = trapz(seg) * dt;
end


%% ========================================================================
%  ASSEMBLE OUTPUT
%  ========================================================================

ev = struct();
ev.start = (footIdx(:) - 1) * dt;
ev.stop  = (stopIdx(:) - 1) * dt;
ev.peak  = (peakIdx(:) - 1) * dt;
ev.amp   = pkVals(:);
ev.dur   = ev.stop - ev.start;
ev.int   = intg;
ev.rise  = riseAmp(:);
ev.noise = sNoise;

end     % EOF
