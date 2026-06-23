function ed = ed_detect(sig, fs, varargin)
% ED_DETECT Detect epileptiform discharges via a moving z-score threshold.
%
%   ed = ED_DETECT(sig, fs, varargin)
%
%   SUMMARY:
%       Detects sharp epileptiform discharges (EDs) on a single signal using
%       a local (moving-window) z-score. Crossings of the z-threshold are
%       localised to the peak deflection, merged within a refractory window,
%       de-duplicated, and gated by a peak-to-peak amplitude criterion.
%       Cleaned port of the core of IED.detect_move_z (LdM), with every
%       constant exposed and the legacy 'negative' direction fixed.
%
%   INPUTS:
%       sig         - (Vec) Signal for detection (e.g. sSig.eeg).
%       fs          - (Num) Sampling frequency [Hz].
%       varargin    - Parameter/Value pairs:
%           'thr'      - (Num) Z-score threshold for the moving-z crossing
%                              and the amplitude gate. {7}
%           'thrDir'   - (Char) 'positive' | 'negative' | 'both'. {'both'}
%           'baseWin'  - (Num) Moving baseline window for mu/sigma [s]. {5}
%           'interDur' - (Num) Refractory / burst-merge window [s]. {0.025}
%           'ampWin'   - (Num) Half-window for the peak-to-peak gate [s]. {0.015}
%           'lowThr'   - (Num) Trough threshold (signal units) for the
%                              twin-peak merge rule. {0.2}
%           'minAmp'   - (Num) Optional absolute amplitude floor (signal
%                              units), applied on top of the z gate. {[]}
%
%   OUTPUTS:
%       ed          - (Struct) Partial detection result:
%           .pos     - (N x 1) Peak sample index into SIG (1-based).
%           .amp     - (N x 1) Peak-to-peak amplitude in +/-ampWin (signal units).
%           .ampZ    - (N x 1) Local moving z-score at the peak.
%           .info    - (Struct) Detection parameters used.
%
%   DEPENDENCIES:
%       peak2peak (Signal Processing Toolbox).
%
%   HISTORY:
%       Created: 22 Jun 2026

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'sig', @isnumeric);
addRequired(p, 'fs', @isnumeric);
addParameter(p, 'thr', 7, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'thrDir', 'both', @(x) any(strcmpi(x, {'positive', 'negative', 'both'})));
addParameter(p, 'baseWin', 5, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'interDur', 0.025, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'ampWin', 0.015, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'lowThr', 0.2, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'minAmp', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));

parse(p, sig, fs, varargin{:});
sig         = p.Results.sig(:);
fs          = p.Results.fs;
thr         = p.Results.thr;
thrDir      = lower(p.Results.thrDir);
baseWin     = p.Results.baseWin;
interDur    = round(p.Results.interDur * fs);   % samples
ampWin      = round(p.Results.ampWin * fs);     % samples
lowThr      = p.Results.lowThr;
minAmp      = p.Results.minAmp;

%% ========================================================================
%  SETUP
%  ========================================================================

ed = struct();
ed.pos = [];
ed.amp = [];
ed.ampZ = [];
ed.info.fs       = fs;
ed.info.thr      = thr;
ed.info.thrDir   = thrDir;
ed.info.baseWin  = baseWin;
ed.info.interDur = p.Results.interDur;
ed.info.ampWin   = p.Results.ampWin;
ed.info.lowThr   = lowThr;
ed.info.minAmp   = minAmp;
ed.info.nDetected = 0;

%% ========================================================================
%  MOVING Z-SCORE
%  ========================================================================

calcWin  = round(baseWin * fs);
sigMu    = movmean(sig, calcWin);
sigSigma = movstd(sig, calcWin);
sigSigma(sigSigma == 0) = eps;
localZ = (sig - sigMu) ./ sigSigma;
localZ(~isfinite(localZ)) = 0;

% Threshold. NOTE: the legacy 'negative' branch used (localZ < thr) which is
% almost always true; corrected here to (localZ < -thr).
switch thrDir
    case 'positive'
        thresholded = localZ > thr;
    case 'negative'
        thresholded = localZ < -thr;
    case 'both'
        thresholded = (localZ > thr) | (localZ < -thr);
end
thresholded = thresholded(:);
cross = find([0; diff(thresholded) > 0]);

% Return early if nothing crosses
if isempty(cross)
    fprintf('ed_detect: no discharges detected\n');
    return
end

%% ========================================================================
%  LOCALISE TO PEAK
%  ========================================================================

% Merge crossings closer than interDur/2 samples. This caps the maximum
% burst rate; the half-interDur keeps the later peak search valid.
ii = find(diff(cross) < interDur / 2);
while ~isempty(ii)
    cross(ii + 1) = [];
    ii = find(diff(cross) < interDur / 2);
end

% Adjust each crossing to the local max / min by magnitude (handles biphasic
% deflections) and drop crossings whose window runs off the signal edge.
nCross  = numel(cross);
peakVal = zeros(nCross, 1);
pos     = zeros(nCross, 1);
rmEdge  = false(nCross, 1);
for iC = 1:nCross
    if cross(iC) + interDur > numel(sig) || cross(iC) - interDur < 1
        rmEdge(iC) = true;
        continue
    end
    seg = sig(cross(iC) - interDur : cross(iC) + interDur);
    [vMax, pMax] = max(seg);
    [vMin, pMin] = min(seg);
    if abs(vMax) >= abs(vMin)
        peakVal(iC) = vMax; relPos = pMax;
    else
        peakVal(iC) = vMin; relPos = pMin;
    end
    pos(iC) = cross(iC) - interDur + relPos - 1;
end
peakVal(rmEdge) = [];
pos(rmEdge) = [];

% Overlapping windows can yield equal / out-of-order positions; collapse them
[pos, idxU] = uniquetol(pos, interDur, 'DataScale', 1);
peakVal = peakVal(idxU);

% Twin-peak merge: if the signal between two peaks never returns below lowThr,
% keep only the larger of the pair (a refractory rule by amplitude).
rmIdx = [];
for iP = 1:numel(pos) - 1
    if min(sig(pos(iP):pos(iP + 1))) > lowThr
        rmIdx = [rmIdx; iP + (peakVal(iP) > peakVal(iP + 1))]; %#ok<AGROW>
    end
end
rmIdx = unique(rmIdx);
pos(rmIdx) = [];
peakVal(rmIdx) = [];
fprintf('ed_detect: %d discharges after detection\n', numel(pos));

%% ========================================================================
%  AMPLITUDE GATE
%  ========================================================================

% Peak-to-peak amplitude in a narrow window must exceed the local threshold
% (thr local-SDs above the local mean), and an optional absolute floor.
thrAmp = thr .* sigSigma + sigMu;
nPos   = numel(pos);
amp    = zeros(nPos, 1);
rmAmp  = false(nPos, 1);
for iP = 1:nPos
    w1 = max(1, pos(iP) - ampWin);
    w2 = min(numel(sig), pos(iP) + ampWin);
    amp(iP) = peak2peak(sig(w1:w2));
    if amp(iP) < thrAmp(pos(iP))
        rmAmp(iP) = true;
    elseif ~isempty(minAmp) && amp(iP) < minAmp
        rmAmp(iP) = true;
    end
end
pos(rmAmp) = [];
amp(rmAmp) = [];
fprintf('ed_detect: %d discharges after amplitude gate\n', numel(pos));

%% ========================================================================
%  ORGANISE OUTPUT
%  ========================================================================

ed.pos  = pos(:);
ed.amp  = amp(:);
ed.ampZ = localZ(ed.pos);
ed.info.nDetected = numel(ed.pos);

end     % EOF
