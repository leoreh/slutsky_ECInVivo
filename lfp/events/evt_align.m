function [wv, shift] = evt_align(wv, tstamps, met, win)
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
%                        'extremum' the largest absolute value.
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
%
%   HISTORY:
%       260722 created for the per-mouse discharge average in ed_wvTbl.

if nargin < 3 || isempty(met), met = 'trough'; end
if nargin < 4 || isempty(win), win = 0.01; end
met = lower(met);

nEv = size(wv, 1);
shift = zeros(nEv, 1);
if strcmp(met, 'none'), return; end

t = tstamps(:)';
[~, iCtr] = min(abs(t));            % sample nearest t = 0
inWin = abs(t) <= win;

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
