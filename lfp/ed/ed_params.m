function ed = ed_params(edSig, ed)
% ED_PARAMS Per-event features for detected epileptiform discharges.
%
%   ed = ED_PARAMS(edSig, ed)
%
%   SUMMARY:
%       Measures the properties the ED pipeline reports and gates on. Two
%       reference frames, because they disagree exactly where it matters:
%       .ampG asks whether a deflection is large for this RECORDING, .ampZ
%       whether it stands out from its own NEIGHBOURHOOD. A modest deflection
%       inside an unusually quiet stretch scores high on one and low on the
%       other, and gating on both is worth more than either alone.
%
%       The neighbourhood is a RING - the samples from CORE out to RING on both
%       sides, with the event cut out - so a feature never compares an event to
%       itself. That is what the old moving mean / SD baseline got wrong.
%
%   INPUTS:
%       edSig - (Struct) ed_detect signal bundle (.lfp .filt .hf).
%       ed    - (Struct) partial detection result; requires .pos (peak sample),
%                        .bouts ([N x 2] candidate bounds) and .info.fs.
%
%   OUTPUTS:
%       ed    - (Struct) with the per-event fields added:
%           .peakTime - (N x 1) peak time [s], window-relative.
%           .times    - (N x 2) measured half-amplitude extent [s].
%           .amp      - (N x 1) peak-to-baseline amplitude of the RAW trace,
%                               SIGNED (signal units); the sign is the
%                               discharge polarity. Reported, never gated.
%           .ampG     - (N x 1) |amp| over the RECORDING's robust scale.
%                               Measured on the raw trace, not the band-passed
%                               one: most of a discharge's amplitude sits below
%                               the discharge band, and gating on the filtered
%                               peak instead costs most of the control-to-
%                               epileptic contrast (26x -> 16x, measured).
%           .ampZ     - (N x 1) band-passed peak over the RING's robust scale.
%           .hfRatio  - (N x 1) peak in the supra-physiological band over the
%                               peak in the discharge band. A discharge scores
%                               ~1; a glitch or movement transient is broadband
%                               and scores several times higher.
%           .dur      - (N x 1) half-amplitude width [ms].
%
%   DEPENDENCIES:
%       None.
%
%   HISTORY:
%       Created: 260622 (half- and 10%-amplitude widths, ported from IED).
%       Updated: 260720 (rebuilt on the ring; .width10 and .pol dropped - the
%                first was never read, the second is now the sign of .amp).

CORE = 0.025;               % half-window holding the event [s]
RING = 1.0;                 % outer half-window of its background [s]

fs   = ed.info.fs;
filt = edSig.filt;
hf   = edSig.hf;
raw  = edSig.lfp;
nSig = numel(filt);

pos = ed.pos(:);
nEv = numel(pos);
nCore = max(1, round(CORE * fs));
nRing = max(nCore + 1, round(RING * fs));

% one scale for the recording, behind .ampG, on a stride: a median and a MAD are
% stable under decimation, and a fixed stride keeps the value reproducible
sub  = raw(1 : max(1, floor(nSig / 2e6)) : end);
gScl = 1.4826 * median(abs(sub - median(sub)));
if ~isfinite(gScl) || gScl <= 0, gScl = 1; end

ed.peakTime = (pos - 1) / fs;
ed.times    = (ed.bouts - 1) / fs;      % refined below where a width is found
ed.amp      = nan(nEv, 1);
ed.ampG     = nan(nEv, 1);
ed.ampZ     = nan(nEv, 1);
ed.hfRatio  = nan(nEv, 1);
ed.dur      = nan(nEv, 1);

for iEv = 1 : nEv

    c1 = max(1, pos(iEv) - nCore);
    c2 = min(nSig, pos(iEv) + nCore);
    r1 = max(1, pos(iEv) - nRing);
    r2 = min(nSig, pos(iEv) + nRing);
    ring = [filt(r1 : c1 - 1); filt(c2 + 1 : r2)];
    if numel(ring) < nCore
        continue                        % too close to an edge to have a ring
    end

    [pkVal, iPk] = max(abs(filt(c1 : c2)));
    iPk = iPk + c1 - 1;
    ringScl = 1.4826 * median(abs(ring - median(ring)));

    ed.amp(iEv)  = raw(iPk) - median(raw(r1 : r2));
    ed.ampG(iEv) = abs(ed.amp(iEv)) / gScl;
    if ringScl > 0
        ed.ampZ(iEv) = pkVal / ringScl;
    end
    if pkVal > 0
        ed.hfRatio(iEv) = max(abs(hf(c1 : c2))) / pkVal;
    end

    % half-amplitude extent, both crossings found separately so the event is
    % free to be asymmetric; without them the candidate's own bounds stand
    [i1, i2] = halfWidth(filt, iPk, c1, c2);
    if ~isempty(i1)
        ed.times(iEv, :) = ([i1, i2] - 1) / fs;
        ed.dur(iEv)      = (i2 - i1) / fs * 1000;
    end
end

end     % EOF


% =========================================================================
%  LOCAL
% =========================================================================
function [i1, i2] = halfWidth(x, iPk, c1, c2)
% Samples where the deflection at IPK last / first crosses half its height.
% Searched over twice the core, so a discharge wider than the peak-search
% window is measured rather than truncated to it. No crossing -> empty.
half = x(iPk) / 2;
span = c2 - c1;
lo   = max(1, c1 - span);
hi   = min(numel(x), c2 + span);

if x(iPk) >= 0
    under = x(lo : hi) < half;
else
    under = x(lo : hi) > half;
end
rel = iPk - lo + 1;

i1 = find(under(1 : rel), 1, 'last');
i2 = find(under(rel : end), 1, 'first');
if isempty(i1) || isempty(i2)
    i1 = []; i2 = [];
    return
end
i1 = i1 + lo - 1;
i2 = i2 + rel - 2 + lo;

end     % halfWidth
