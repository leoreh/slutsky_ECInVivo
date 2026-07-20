function ed = ed_params(edSig, ed)
% ED_PARAMS Per-event features for detected epileptiform discharges.
%
%   ed = ED_PARAMS(edSig, ed)
%
%   SUMMARY:
%       Measures the three properties that define a discharge, plus the ones
%       worth reporting. Each gate asks one question, about a different way of
%       being wrong:
%
%       .fastZ  IS IT SHARP? The peak of the 60-150 Hz trace in units of that
%               band's own background. A discharge is a near-instantaneous
%               transient, so it is enormous here (35-45x); a sharp wave or any
%               other slow LFP deflection is a smooth excursion with almost no
%               energy this high (6-12x), however large it looks raw. This is
%               what separates a discharge from normal hippocampal activity.
%
%       .posZ   IS IT AN UPWARD DISCHARGE? The positive excursion of the RAW
%               trace above its own pre-event baseline, in units of the
%               recording's background. A discharge rises well above baseline
%               (~14x); a step artifact - the trace dropping and never
%               returning, which is the other thing that is sharp - has
%               essentially no positive component (~1.5x). Measured against a
%               baseline taken BEFORE the event, so a step cannot hide by
%               moving the local mean with it.
%
%       .isoZ   IS IT ALONE? The same peak over the robust scale of a RING -
%               the samples out to RING on both sides with the event itself cut
%               out. A discharge is a solitary spike on quiet background
%               (~47x); the other big confusion is a noisy or oscillatory
%               stretch in which the "event" is merely the largest of many
%               comparable peaks, and that scores low however sharp and
%               positive it is.
%
%       Each is blind to the others' failure mode, which is why all three are
%       needed and why three are enough: sharpness cannot see a step (AUC 0.68
%       against them), positivity cannot see a slow deflection, and neither
%       can see a noisy epoch. Measured on 47 discharges curated by hand in
%       three raMCU mice.
%
%   INPUTS:
%       edSig - (Struct) ed_detect signal bundle (.lfp .filt .sclFast .sclRaw).
%       ed    - (Struct) partial detection result; requires .pos (peak sample),
%                        .bouts ([N x 2] candidate bounds) and .info.fs.
%
%   OUTPUTS:
%       ed    - (Struct) with the per-event fields added:
%           .peakTime - (N x 1) peak time [s], window-relative.
%           .times    - (N x 2) measured half-amplitude extent [s].
%           .fastZ    - (N x 1) sharpness (see above).
%           .posZ     - (N x 1) positive discharge amplitude (see above).
%           .isoZ     - (N x 1) isolation from the local background (see above).
%           .amp      - (N x 1) peak-to-baseline amplitude of the raw trace,
%                               SIGNED (signal units). Reported, not gated;
%                               its sign is the dominant deflection's polarity.
%           .dur      - (N x 1) half-amplitude width [ms].
%
%   DEPENDENCIES:
%       None.
%
%   HISTORY:
%       Created: 260622 (half- and 10%-amplitude widths, ported from IED).
%       Updated: 260720 (rebuilt on an event-excluded ring).
%       Updated: 260721 (rebuilt again against 47 CURATED discharges. hfRatio
%                and ampG went with the 10-100 Hz detection band that produced
%                them - a band in which discharges and normal deflections are
%                not separable at all. The ring survived the move, remeasured
%                in the new band. See dev/ed_pipeline_rebuild.md.)

CORE = 0.020;               % half-window holding the event [s]
BASE = [-0.25 -0.10];       % pre-event baseline window [s]
RING = 1.0;                 % outer half-window of the background ring [s]

fs   = ed.info.fs;
raw  = edSig.lfp;
filt = edSig.filt;
nSig = numel(raw);

pos = ed.pos(:);
nEv = numel(pos);
nCore = max(1, round(CORE * fs));
nRing = max(nCore + 1, round(RING * fs));
bSamp = round(BASE * fs);

ed.peakTime = (pos - 1) / fs;
ed.times    = (ed.bouts - 1) / fs;      % refined below where a width is found
ed.fastZ    = nan(nEv, 1);
ed.posZ     = nan(nEv, 1);
ed.isoZ     = nan(nEv, 1);
ed.amp      = nan(nEv, 1);
ed.dur      = nan(nEv, 1);

for iEv = 1 : nEv

    p = pos(iEv);
    if p + bSamp(1) < 1 || p + nCore > nSig
        continue                        % no room for the baseline or the core
    end
    c1 = p - nCore;
    c2 = p + nCore;

    pkFast = max(abs(filt(c1 : c2)));
    ed.fastZ(iEv) = pkFast / edSig.sclFast;

    r1 = max(1, p - nRing);
    r2 = min(nSig, p + nRing);
    ring = [filt(r1 : c1 - 1); filt(c2 + 1 : r2)];
    if numel(ring) >= nCore
        sclRing = 1.4826 * median(abs(ring - median(ring)));
        if sclRing > 0, ed.isoZ(iEv) = pkFast / sclRing; end
    end

    base = median(raw(p + bSamp(1) : p + bSamp(2)));
    seg  = raw(c1 : c2) - base;
    ed.posZ(iEv) = max(seg) / edSig.sclRaw;

    [~, iPk] = max(abs(seg));
    ed.amp(iEv) = seg(iPk);

    % half-amplitude extent, both crossings measured separately so the event is
    % free to be asymmetric; without them the candidate's own bounds stand
    [i1, i2] = halfWidth(seg, iPk);
    if ~isempty(i1)
        ed.times(iEv, :) = ([i1, i2] + c1 - 2) / fs;
        ed.dur(iEv)      = (i2 - i1) / fs * 1000;
    end
end

end     % EOF


% =========================================================================
%  LOCAL
% =========================================================================
function [i1, i2] = halfWidth(seg, iPk)
% Samples where the deflection at IPK last / first crosses half its height.
% Indices are into SEG. No crossing inside the core returns empty, and the
% caller keeps the candidate's own bounds.
half = seg(iPk) / 2;
if seg(iPk) >= 0
    under = seg < half;
else
    under = seg > half;
end

i1 = find(under(1 : iPk), 1, 'last');
i2 = find(under(iPk : end), 1, 'first');
if isempty(i1) || isempty(i2)
    i1 = []; i2 = [];
    return
end
i2 = i2 + iPk - 1;

end     % halfWidth
