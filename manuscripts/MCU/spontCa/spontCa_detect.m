function ev = spontCa_detect(trace, fs, varargin)
% SPONTCA_DETECT Detects Ca transients in a single dF/F trace.
%
%   ev = SPONTCA_DETECT(TRACE, FS, ...) operates on one trace using
%   FINDPEAKS with a prominence (local-rise) threshold, then extends each
%   peak's start and stop via hysteresis to the first sample where the
%   trace drops below thrBsl * peakValue.
%
%   The input is assumed to be dF/F already (centered near zero outside
%   events). No rolling baseline is subtracted. Noise is the robust SD of
%   the differenced signal (Allan-Pettersson) and feeds the MinPeakHeight.
%
%   PIPELINE:
%       1. Peak picking: FINDPEAKS with
%             MinPeakHeight     = max(kThr * noise, minAmp)
%             MinPeakProminence = minRise
%             MinPeakDistance   = round(minIEI * fs)
%          Prominence is the standard "local rise" measure - the height
%          of a peak above the higher of its two flanking troughs. This
%          rejects fluctuations on decay tails and plateaus (their
%          flanking troughs sit at the same elevated level, so prominence
%          is tiny) while preserving real events sitting on top of an
%          elevation if their excursion is large enough.
%       2. Hysteresis extension. From each peak, walk left/right until
%          trace drops below thrBsl * peakValue. Walks are capped at the
%          neighbouring peaks so events do not bleed into each other.
%       3. Drop events shorter than minDur (on the extended span).
%
%   INPUTS:
%       trace   - (1 x nT) dF/F signal
%       fs      - (scalar) sampling rate (Hz)
%
%   OPTIONAL (Name-Value):
%       'kThr'    - (num) MinPeakHeight in noise SDs       {3}
%       'minAmp'  - (num) absolute floor on peak amp (dF/F) {0.05}
%       'minRise' - (num) MinPeakProminence (dF/F)         {0.05}
%       'minDur'  - (num) min event duration (s)           {0.4}
%       'minIEI'  - (num) min peak-to-peak distance (s)    {0.4}
%       'thrBsl'  - (num) return-to fraction of peak for start/stop {0.3}
%
%   OUTPUT:
%       ev struct with fields:
%         .start  (n x 1)  event start times (s) (hysteresis foot)
%         .stop   (n x 1)  event stop times (s)  (hysteresis tail)
%         .peak   (n x 1)  event peak times (s)
%         .amp    (n x 1)  peak amp (dF/F)
%         .dur    (n x 1)  duration (s)
%         .int    (n x 1)  integral (dF/F * s) over the extended span
%         .rise   (n x 1)  prominence from FINDPEAKS (dF/F)
%         .noise  (scalar) robust noise std
%
%   See also: SPONTCA_EVENTS, SPONTCA_COUPLE, SPONTCA_LOAD, FINDPEAKS

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'trace', @(x) isnumeric(x) && isvector(x));
addRequired(p, 'fs',    @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'kThr',    3,    @isnumeric);
addParameter(p, 'minAmp',  0.05, @isnumeric);
addParameter(p, 'minRise', 0.05, @isnumeric);
addParameter(p, 'minDur',  0.4,  @isnumeric);
addParameter(p, 'minIEI',  0.4,  @isnumeric);
addParameter(p, 'thrBsl',  0.3,  @isnumeric);
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
%  PEAK PICKING (FINDPEAKS with prominence)
%  ========================================================================

minHeight = max(P.kThr * sNoise, P.minAmp);
minDist   = max(1, round(P.minIEI * fs));

% FINDPEAKS doesn't tolerate NaN; substitute -Inf so they cannot be peaks.
traceClean = trace;
traceClean(isnan(traceClean)) = -Inf;

[pkVals, pkLocs, ~, proms] = findpeaks(traceClean, ...
    'MinPeakHeight',     minHeight, ...
    'MinPeakProminence', P.minRise, ...
    'MinPeakDistance',   minDist);

pkVals = pkVals(:);
pkLocs = pkLocs(:);
proms  = proms(:);
nEv    = length(pkLocs);


%% ========================================================================
%  HYSTERESIS EXTENSION
%  ========================================================================

startIdx = zeros(nEv, 1);
stopIdx  = zeros(nEv, 1);

for iE = 1:nEv
    pk = pkLocs(iE);
    returnThr = P.thrBsl * pkVals(iE);

    % Backward walk bound: the trough between the previous peak and this
    % one. Walks stop at the trough sample (inclusive on this side,
    % exclusive on the other), so adjacent events do not share samples.
    if iE > 1
        prevPk = pkLocs(iE - 1);
        [~, troughRel] = min(trace(prevPk:pk), [], 'omitnan');
        leftBound = prevPk + troughRel - 1;
    else
        leftBound = 1;
    end
    s = pk;
    while s > leftBound
        prev = trace(s - 1);
        if isnan(prev) || prev <= returnThr
            break;
        end
        s = s - 1;
    end
    startIdx(iE) = s;

    % Forward walk bound: one sample before the trough between this peak
    % and the next, so the trough belongs to the next event only.
    if iE < nEv
        nextPk = pkLocs(iE + 1);
        [~, troughRel] = min(trace(pk:nextPk), [], 'omitnan');
        rightBound = max(pk, pk + troughRel - 2);
    else
        rightBound = nT;
    end
    e = pk;
    while e < rightBound
        nxt = trace(e + 1);
        if isnan(nxt) || nxt <= returnThr
            break;
        end
        e = e + 1;
    end
    stopIdx(iE) = e;
end


%% ========================================================================
%  MIN-DURATION FILTER (on extended span)
%  ========================================================================

durSmp = stopIdx - startIdx + 1;
keep   = durSmp >= max(2, round(P.minDur * fs));
startIdx = startIdx(keep);
stopIdx  = stopIdx(keep);
pkLocs   = pkLocs(keep);
pkVals   = pkVals(keep);
proms    = proms(keep);
nEv      = length(pkLocs);


%% ========================================================================
%  INTEGRAL OVER EXTENDED SPAN
%  ========================================================================

intg = nan(nEv, 1);
for iE = 1:nEv
    seg = trace(startIdx(iE):stopIdx(iE));
    seg(isnan(seg)) = 0;
    intg(iE) = trapz(seg) * dt;
end


%% ========================================================================
%  ASSEMBLE OUTPUT
%  ========================================================================

ev = struct();
ev.start = (startIdx(:) - 1) * dt;
ev.stop  = (stopIdx(:)  - 1) * dt;
ev.peak  = (pkLocs(:)   - 1) * dt;
ev.amp   = pkVals;
ev.dur   = ev.stop - ev.start;
ev.int   = intg;
ev.rise  = proms;
ev.noise = sNoise;

end     % EOF
