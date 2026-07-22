function [wv, shift, met] = evt_align(wv, tstamps, met, win)
% EVT_ALIGN Shift event waveforms so a chosen extremum sits at t = 0.
%
%   [wv, shift] = EVT_ALIGN(wv, tstamps, met, win)
%
%   SUMMARY:
%       Re-centres each snippet on its own trough (or peak) by an integer-
%       sample shift, for a cleaner average. Detection aligns to the largest
%       ABSOLUTE excursion, so in a biphasic event some rows are centred on
%       their positive peak and others on their negative trough - and averaging
%       that mix smears the mean and leaves a notch at t = 0. Fixing every row
%       to the same feature lines the events up.
%
%       A SHIFT, not a re-extraction. The snippet already spans wider than the
%       window worth showing, so the feature is found and the row is slid within
%       the samples already in hand - no second read of the LFP. Vacated samples
%       at the edge become NaN rather than WRAPPING (what circshift would do):
%       a wrap moves the far edge of the trace under the near edge, which is
%       wrong signal, whereas a small shift only ever exposes NaN out where the
%       display never reaches.
%
%       This is a DISPLAY choice and it assumes a feature of that sign exists.
%       'trough' on an event with no real negative component aligns to a noise
%       dip; the pipeline stays polarity-blind, and only a caller that has
%       looked at the data should ask for a fixed polarity (see ed_wvTbl).
%
%   INPUTS:
%       wv      - <mat>  [nEv x nSamp] per-event waveforms, one row per event.
%       tstamps - <vec>  [1 x nSamp] window time base (s), from evt_maps.
%       met     - <char> feature to centre on:
%                        'trough'   the minimum (negative-going event).
%                        'peak'     the maximum.
%                        'extremum' the largest absolute value, per event.
%                        'auto'     ONE choice for the whole set: the polarity
%                                   MOST events have. Each event's own biggest
%                                   swing (positive or negative) is one vote;
%                                   the majority wins, ties go to trough.
%                                   Polarity is a property of the electrode's
%                                   layer, so it is shared by a session's
%                                   events - decided once here, not per event,
%                                   which 'extremum' does and which lets noise
%                                   flip a near-flat event's sign. A VOTE, not
%                                   the sign of the mean waveform: a few large
%                                   opposite-going events inflate the mean's
%                                   other lobe and flip it (on raMCU4 the mean
%                                   reads peak while 28 of 35 events are
%                                   troughs, which the vote and the median
%                                   waveform both get right).
%                        'none'     leave the rows as they are.
%                        {'trough'}
%       win     - <num>  search half-window around t = 0 (s): the feature is
%                        taken inside +-win, so a deflection elsewhere in the
%                        snippet cannot capture the alignment. {0.01}
%
%   OUTPUTS:
%       wv      - <mat>  [nEv x nSamp] shifted; an all-NaN row is left as is.
%       shift   - <vec>  [nEv x 1] samples each row moved (0 for 'none' and for
%                        rows with no finite sample in the window).
%       met     - <char> the feature actually used, so 'auto' reports whether
%                        it resolved to 'peak' or 'trough'.
%
%   HISTORY:
%       260722 created for the per-mouse discharge average in ed_wvTbl.
%       260722b 'auto' picks the polarity from the average waveform.

if nargin < 3 || isempty(met), met = 'trough'; end
if nargin < 4 || isempty(win), win = 0.01; end
met = lower(met);

nEv = size(wv, 1);
shift = zeros(nEv, 1);
if strcmp(met, 'none'), return; end

t = tstamps(:)';
[~, iCtr] = min(abs(t));            % sample nearest t = 0
inWin = abs(t) <= win;

% resolve 'auto' ONCE: the polarity the majority of events actually have. Each
% event's biggest swing inside the window is one vote for peak or trough. A
% vote, not the sign of the MEAN waveform, because a few large opposite-going
% events inflate the mean's other lobe and flip it - on raMCU4 the mean reads
% peak while 28 of 35 events are troughs.
if strcmp(met, 'auto')
    seg = wv;
    seg(:, ~inWin) = NaN;
    [~, iMx] = max(abs(seg), [], 2);
    val = wv(sub2ind(size(wv), (1 : nEv)', iMx));
    val = val(isfinite(val));       % drop all-NaN edge rows
    met = 'trough';
    if sum(val > 0) > sum(val < 0), met = 'peak'; end
end

for iEv = 1 : nEv
    row = wv(iEv, :);
    seg = row;
    seg(~inWin) = NaN;              % restrict the search to +-win
    if all(isnan(seg)), continue; end

    switch met
        case 'trough',   [~, iPk] = min(seg);
        case 'peak',     [~, iPk] = max(seg);
        case 'extremum', [~, iPk] = max(abs(seg));
        otherwise
            error('evt_align:met', 'unknown met "%s"', met);
    end

    s = iCtr - iPk;                 % move iPk to the centre
    if s == 0, continue; end
    shift(iEv) = s;

    % output col j takes input col j - s; the rest is NaN, never wrapped
    out = nan(size(row));
    src = (1 : numel(row)) - s;
    ok  = src >= 1 & src <= numel(row);
    out(ok) = row(src(ok));
    wv(iEv, :) = out;
end

end     % EOF
