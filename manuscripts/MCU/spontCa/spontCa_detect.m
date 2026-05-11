function ev = spontCa_detect(trace, fs, varargin)
% SPONTCA_DETECT Detects Ca transients in a single dF/F trace.
%
%   ev = SPONTCA_DETECT(TRACE, FS, ...) operates on one trace using
%   threshold-based connected-component detection. Designed to be composed
%   by SPONTCA_EVENTS or any caller that loops over cells.
%
%   The input is assumed to be dF/F already (centered near zero outside
%   events). No rolling baseline is subtracted. Noise is the robust SD of
%   the differenced signal (Allan-Pettersson), and the threshold is
%   max(kThr * noise, minAmp). Every above-threshold excursion is one
%   event - no FINDPEAKS, no prominence, no over-segmentation on slow
%   mito decays.
%
%   INPUTS:
%       trace   - (1 x nT) dF/F signal
%       fs      - (scalar) sampling rate (Hz)
%
%   OPTIONAL (Name-Value):
%       'kThr'   - (num) threshold in robust noise SDs   {3}
%       'minAmp' - (num) hard floor on amplitude (dF/F)  {0.02}
%       'minDur' - (num) min event duration (s)          {0.4}
%       'minIEI' - (num) min inter-event gap (s)         {0.4}
%
%   OUTPUT:
%       ev struct with fields:
%         .start  (n x 1)  event start times (s)
%         .stop   (n x 1)  event stop times (s)
%         .amp    (n x 1)  peak amp (dF/F)
%         .dur    (n x 1)  duration (s)
%         .int    (n x 1)  integral (dF/F * s)
%         .noise  (scalar) robust noise std
%
%   See also: SPONTCA_EVENTS, SPONTCA_LOAD

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'trace', @(x) isnumeric(x) && isvector(x));
addRequired(p, 'fs',    @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'kThr',   3,    @isnumeric);
addParameter(p, 'minAmp', 0.02, @isnumeric);
addParameter(p, 'minDur', 0.4,  @isnumeric);
addParameter(p, 'minIEI', 0.4,  @isnumeric);
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
%  PER-EVENT METRICS
%  ========================================================================

nEv = length(startIdx);
ev = struct();
ev.start = (startIdx(:) - 1) * dt;
ev.stop  = (stopIdx(:)  - 1) * dt;
ev.dur   = ev.stop - ev.start;
ev.amp   = nan(nEv, 1);
ev.int   = nan(nEv, 1);

for iE = 1:nEv
    seg = trace(startIdx(iE):stopIdx(iE));
    seg(isnan(seg)) = 0;
    ev.amp(iE) = max(seg);
    ev.int(iE) = trapz(seg) * dt;
end

ev.noise = sNoise;

end     % EOF
