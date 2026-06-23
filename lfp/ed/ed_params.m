function ed = ed_params(sig, ed, varargin)
% ED_PARAMS Per-event features for detected epileptiform discharges.
%
%   ed = ED_PARAMS(sig, ed, varargin)
%
%   SUMMARY:
%       Fills the per-event timing and shape features of the ed struct.
%       Ported from IED.@data/calculate_props + extract_discharges (LdM):
%       amplitude is measured peak-to-baseline (local moving mean), and the
%       half- and 10%-amplitude widths are found by threshold crossings on a
%       clip centred on the peak. The event extent (.times) is defined as
%       peakTime +/- dur/2 so it slots into ripp_states / rate code.
%
%   INPUTS:
%       sig         - (Vec) Same signal passed to ed_detect.
%       ed          - (Struct) Output of ed_detect (requires .pos, .info.fs).
%       varargin    - Parameter/Value pairs:
%           'marg'       - (Num) Half-window clipped around each peak [s]. {0.05}
%           'ampBaseWin' - (Num) Moving-mean window for the amplitude
%                                baseline [s]. {1}
%
%   OUTPUTS:
%       ed          - (Struct) With added fields:
%           .peakTime - (N x 1) Peak time [s] (relative; wrapper shifts to absolute).
%           .times    - (N x 2) Event start/end [s] = peakTime +/- dur/2.
%           .dur      - (N x 1) Half-amplitude width [ms].
%           .width10  - (N x 1) 10%-amplitude width [ms].
%           .amp      - (N x 1) Peak-to-peak amplitude (kept from ed_detect).
%
%   DEPENDENCIES:
%       None.
%
%   HISTORY:
%       Created: 22 Jun 2026

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'sig', @isnumeric);
addRequired(p, 'ed', @isstruct);
addParameter(p, 'marg', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'ampBaseWin', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);

parse(p, sig, ed, varargin{:});
sig        = p.Results.sig(:);
ed         = p.Results.ed;
marg       = p.Results.marg;
ampBaseWin = p.Results.ampBaseWin;

%% ========================================================================
%  SETUP
%  ========================================================================

fs      = ed.info.fs;
pos     = ed.pos(:);
nEvents = numel(pos);
nSig    = numel(sig);

% Peak time, relative to the (windowed) signal; ed_wrapper shifts to absolute
ed.peakTime = (pos - 1) / fs;

% Initialise outputs
ed.dur     = nan(nEvents, 1);
ed.width10 = nan(nEvents, 1);
if ~isfield(ed, 'amp') || numel(ed.amp) ~= nEvents
    ed.amp = nan(nEvents, 1);
end

% Baseline for the amplitude reference, and clip half-window
winMeans = movmean(sig, round(ampBaseWin * fs));
margSamp = floor(marg * fs);
detP     = margSamp + 1;            % index of the peak within a clip

%% ========================================================================
%  PER-EVENT FEATURES
%  ========================================================================

for iEv = 1:nEvents
    idx = (pos(iEv) - margSamp):(pos(iEv) + margSamp);
    if idx(1) < 1 || idx(end) > nSig
        continue                    % out-of-bound clip: leave as nan
    end
    seg   = sig(idx);
    peakV = seg(detP);
    baseV = winMeans(pos(iEv));
    ampBase = abs(peakV - baseV);

    halfLvl  = ampBase / 2;
    tenthLvl = ampBase / 10;
    if peakV > baseV
        ed.dur(iEv)     = width_from_thr(seg, detP, baseV + halfLvl,  @gt, fs);
        ed.width10(iEv) = width_from_thr(seg, detP, baseV + tenthLvl, @gt, fs);
    else
        ed.dur(iEv)     = width_from_thr(seg, detP, baseV - halfLvl,  @lt, fs);
        ed.width10(iEv) = width_from_thr(seg, detP, baseV - tenthLvl, @lt, fs);
    end
end

% Widths to ms
ed.dur     = ed.dur * 1000;
ed.width10 = ed.width10 * 1000;

%% ========================================================================
%  EVENT EXTENT
%  ========================================================================

% Half-amplitude width defines the shaded extent; events with no measurable
% width collapse to a degenerate [peakTime peakTime] interval (finite).
halfDurSec = (ed.dur / 1000) / 2;
halfDurSec(isnan(halfDurSec)) = 0;
ed.times = [ed.peakTime - halfDurSec, ed.peakTime + halfDurSec];

end     % EOF


%% ========================================================================
%  LOCAL: WIDTH FROM THRESHOLD
%  ========================================================================

function w = width_from_thr(seg, detP, crossVal, crossFun, fs)
% Width between the threshold crossings flanking the peak (in seconds).
right = crossFun(seg(detP:end), crossVal);
left  = crossFun(seg(1:detP),  crossVal);
lMark = find(diff(left),  1, 'last')  + 1;
rMark = find(diff(right), 1, 'first') + detP - 1;
if isempty(lMark) || isempty(rMark)
    w = nan;
else
    w = (rMark - lMark) / fs;
end
end
