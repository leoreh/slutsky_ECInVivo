function ev = spontCa_detect(trace, fs, varargin)
% SPONTCA_DETECT Detects Ca transients in a single dF/F trace.
%
%   ev = SPONTCA_DETECT(TRACE, FS, ...) picks local maxima with FINDPEAKS
%   (height + distance only), filters them by a LEFT-SIDE local rise gate
%   (median of a pre-peak window), then extends each peak's start and stop
%   via hysteresis to the first sample where the trace drops below an
%   absolute return threshold (close to baseline).
%
%   The input is assumed to be dF/F already (centered near zero outside
%   events). No rolling baseline is subtracted. Noise is the robust SD of
%   the differenced signal (Allan-Pettersson) and feeds the MinPeakHeight.
%
%   WHY LEFT-SIDE PROMINENCE
%   ------------------------
%   Standard FINDPEAKS prominence uses the HIGHER of the two flanking
%   troughs (walks both sides until a higher peak is found). For an
%   asymmetric event (fast rise, slow decay) followed by a plateau, the
%   right base is the plateau level - so a lone real event whose decay
%   merges into the plateau gets a tiny standard prominence and is
%   rejected. Conversely a noise bump on top of a plateau between two
%   big peaks shares a low left base far back at baseline, inflating
%   prominence if you used the lower base. Neither is what we want.
%   The local left-side rise (peak - median of a windowed pre-peak
%   baseline) is the correct measure: it depends only on the rise itself.
%
%   WHY ABSOLUTE RETURN THRESHOLD
%   -----------------------------
%   Hysteresis using thrBsl * peakValue scales the cutoff with peak
%   amplitude. A 0.35 peak stops at 0.105 (cuts the decay short while
%   still high) while a 0.08 peak stops at 0.024 (extends into noise).
%   An absolute thrBsl (dF/F) returns every event to the same baseline
%   regardless of peak height, which matches what "event over" means
%   biologically.
%
%   PIPELINE:
%       1. FINDPEAKS with MinPeakHeight = max(kThr * noise, minAmp) and
%          MinPeakDistance = round(minIEI * fs). No prominence filter.
%       2. Local left-rise filter. For each peak, look back to the window
%          [peak - tLead - tBase, peak - tLead] and compute the median.
%          Require peak - median >= minRise.
%       3. Hysteresis extension from peak to start and to stop until the
%          trace drops below the absolute thrBsl. Walks are bounded by the
%          troughs between consecutive peaks so adjacent events do not
%          share samples.
%       4. Drop events shorter than minDur on the extended span.
%
%   INPUTS:
%       trace   - (1 x nT) dF/F signal
%       fs      - (scalar) sampling rate (Hz)
%
%   OPTIONAL (Name-Value):
%       'kThr'    - (num) MinPeakHeight in noise SDs           {3}
%       'minAmp'  - (num) absolute floor on peak amp (dF/F)    {0.05}
%       'minRise' - (num) min local rise above pre-peak median {0.05}
%       'minDur'  - (num) min event duration (s)               {0.4}
%       'minIEI'  - (num) min peak-to-peak distance (s)        {0.4}
%       'tBase'   - (num) pre-peak baseline window length (s)  {1.5}
%       'tLead'   - (num) offset back from peak before window (s) {0.5}
%       'thrBsl'  - (num) absolute return threshold (dF/F)     {0.02}
%
%   OUTPUT:
%       ev struct with fields:
%         .start  (n x 1)  event start times (s) (hysteresis foot)
%         .stop   (n x 1)  event stop times (s)  (hysteresis tail)
%         .peak   (n x 1)  event peak times (s)
%         .amp    (n x 1)  peak amp (dF/F)
%         .dur    (n x 1)  duration (s)
%         .int    (n x 1)  integral (dF/F * s) over the extended span
%         .rise   (n x 1)  peak - pre-peak median (dF/F)
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
addParameter(p, 'tBase',   1.5,  @isnumeric);
addParameter(p, 'tLead',   0.5,  @isnumeric);
addParameter(p, 'thrBsl',  0.02, @isnumeric);
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
%  LOCAL MAXIMA (height + distance only)
%  ========================================================================

minHeight = max(P.kThr * sNoise, P.minAmp);
minDist   = max(1, round(P.minIEI * fs));

traceClean = trace;
traceClean(isnan(traceClean)) = -Inf;

[pkVals, pkLocs] = findpeaks(traceClean, ...
    'MinPeakHeight',   minHeight, ...
    'MinPeakDistance', minDist);
pkVals = pkVals(:);
pkLocs = pkLocs(:);


%% ========================================================================
%  LOCAL LEFT-RISE FILTER
%  ========================================================================

tBaseSmp = max(1, round(P.tBase * fs));
tLeadSmp = max(1, round(P.tLead * fs));

nEv  = length(pkLocs);
rise = nan(nEv, 1);
keep = false(nEv, 1);
for iE = 1:nEv
    pk   = pkLocs(iE);
    winR = pk - tLeadSmp;
    winL = pk - tLeadSmp - tBaseSmp + 1;
    winL = max(1, winL);
    if winR < 1
        continue;
    end
    base = median(trace(winL:winR), 'omitnan');
    if isnan(base)
        continue;
    end
    rise(iE) = pkVals(iE) - base;
    keep(iE) = rise(iE) >= P.minRise;
end

pkVals = pkVals(keep);
pkLocs = pkLocs(keep);
rise   = rise(keep);
nEv    = length(pkLocs);


%% ========================================================================
%  HYSTERESIS EXTENSION (absolute return threshold)
%  ========================================================================

startIdx = zeros(nEv, 1);
stopIdx  = zeros(nEv, 1);

for iE = 1:nEv
    pk = pkLocs(iE);

    % Backward walk bound: trough between previous peak and this one.
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
        if isnan(prev) || prev <= P.thrBsl
            break;
        end
        s = s - 1;
    end
    startIdx(iE) = s;

    % Forward walk bound: one sample before the trough so it belongs to
    % the next event only.
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
        if isnan(nxt) || nxt <= P.thrBsl
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
rise     = rise(keep);
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
ev.rise  = rise;
ev.noise = sNoise;

end     % EOF
