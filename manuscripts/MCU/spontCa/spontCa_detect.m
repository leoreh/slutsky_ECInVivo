function ev = spontCa_detect(trace, fs, varargin)
% SPONTCA_DETECT Detects Ca transients in a single dF/F trace.
%
%   ev = SPONTCA_DETECT(TRACE, FS, ...) operates on one trace using
%   threshold-based connected-component detection plus a local-rise gate.
%   Designed to be composed by SPONTCA_EVENTS or any caller that loops over
%   cells.
%
%   The input is assumed to be dF/F already (centered near zero outside
%   events). No rolling baseline is subtracted. Noise is the robust SD of
%   the differenced signal (Allan-Pettersson).
%
%   PIPELINE:
%       1. Absolute-amplitude threshold: trace > max(kThr * noise, minAmp).
%          Connected components define candidate events.
%       2. Merge components closer than minIEI.
%       3. Drop components shorter than minDur.
%       4. LOCAL-RISE GATE: for each candidate, locate the peak sample and
%          compute base = median(trace in [peak - tLead - tRise,
%          peak - tLead]). Require peak - base >= minRise. This rejects
%          fluctuations on decay tails and on elevated plateaus, which
%          have small or negative local rise even though their absolute
%          value is above threshold. tLead skips the rise itself (the peak
%          is not truly instantaneous); tRise is the pre-event baseline
%          window.
%
%   INPUTS:
%       trace   - (1 x nT) dF/F signal
%       fs      - (scalar) sampling rate (Hz)
%
%   OPTIONAL (Name-Value):
%       'kThr'    - (num) absolute threshold in noise SDs    {3}
%       'minAmp'  - (num) absolute floor on amplitude (dF/F) {0.02}
%       'minDur'  - (num) min event duration (s)             {0.4}
%       'minIEI'  - (num) min inter-event gap (s)            {0.4}
%       'minRise' - (num) min rise above local baseline (dF/F) {0.05}
%       'tLead'   - (num) offset back from peak before window (s) {0.5}
%       'tRise'   - (num) length of pre-peak baseline window (s)  {1.0}
%
%   OUTPUT:
%       ev struct with fields:
%         .start  (n x 1)  event start times (s)
%         .stop   (n x 1)  event stop times (s)
%         .amp    (n x 1)  peak amp (dF/F)
%         .dur    (n x 1)  duration (s)
%         .int    (n x 1)  integral (dF/F * s)
%         .rise   (n x 1)  peak - local baseline (dF/F)
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
addParameter(p, 'minAmp',  0.02, @isnumeric);
addParameter(p, 'minDur',  0.4,  @isnumeric);
addParameter(p, 'minIEI',  0.4,  @isnumeric);
addParameter(p, 'minRise', 0.05, @isnumeric);
addParameter(p, 'tLead',   0.5,  @isnumeric);
addParameter(p, 'tRise',   1.0,  @isnumeric);
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
%  CONNECTED-COMPONENT DETECTION
%  ========================================================================

thr  = max(P.kThr * sNoise, P.minAmp);
mask = trace > thr;
mask(isnan(trace)) = false;

edges = diff([0, mask, 0]);
startIdx = find(edges == 1);
stopIdx  = find(edges == -1) - 1;

% Greedy merge of events whose gap is shorter than minIEI
gapSmp = max(1, round(P.minIEI * fs));
mS = [];
mE = [];
for k = 1:length(startIdx)
    if ~isempty(mE) && (startIdx(k) - mE(end) - 1 < gapSmp)
        mE(end) = stopIdx(k);
    else
        mS(end+1) = startIdx(k);      %#ok<AGROW>
        mE(end+1) = stopIdx(k);       %#ok<AGROW>
    end
end
startIdx = mS;
stopIdx  = mE;

% Min duration filter
keep = (stopIdx - startIdx + 1) >= max(2, round(P.minDur * fs));
startIdx = startIdx(keep);
stopIdx  = stopIdx(keep);


%% ========================================================================
%  PER-EVENT METRICS (peak amp, integral, peak index)
%  ========================================================================

nEv0  = length(startIdx);
amp   = nan(nEv0, 1);
intg  = nan(nEv0, 1);
pkIdx = nan(nEv0, 1);

for iE = 1:nEv0
    seg = trace(startIdx(iE):stopIdx(iE));
    segNoNaN = seg;
    segNoNaN(isnan(segNoNaN)) = -inf;
    [pkVal, pkRel] = max(segNoNaN);
    pkIdx(iE) = startIdx(iE) + pkRel - 1;
    amp(iE)   = pkVal;
    segZ = seg;
    segZ(isnan(segZ)) = 0;
    intg(iE) = trapz(segZ) * dt;
end


%% ========================================================================
%  LOCAL-RISE GATE
%  ========================================================================

leadSmp = max(1, round(P.tLead * fs));
winSmp  = max(1, round(P.tRise * fs));

rise = nan(nEv0, 1);
keepRise = false(nEv0, 1);
for iE = 1:nEv0
    winR = pkIdx(iE) - leadSmp;
    winL = winR - winSmp + 1;
    winL = max(1, winL);
    if winR < 1
        continue;
    end
    base = median(trace(winL:winR), 'omitnan');
    if isnan(base)
        continue;
    end
    rise(iE) = amp(iE) - base;
    keepRise(iE) = rise(iE) >= P.minRise;
end

startIdx = startIdx(keepRise);
stopIdx  = stopIdx(keepRise);
amp      = amp(keepRise);
intg     = intg(keepRise);
rise     = rise(keepRise);


%% ========================================================================
%  ASSEMBLE OUTPUT
%  ========================================================================

ev = struct();
ev.start = (startIdx(:) - 1) * dt;
ev.stop  = (stopIdx(:)  - 1) * dt;
ev.dur   = ev.stop - ev.start;
ev.amp   = amp;
ev.int   = intg;
ev.rise  = rise;
ev.noise = sNoise;

end     % EOF
