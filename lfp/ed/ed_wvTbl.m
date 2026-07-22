function [tbl, tstamps] = ed_wvTbl(basepaths, varargin)
% ED_WVTBL Accepted discharge waveforms across sessions, one row per event.
%
%   [tbl, tstamps] = ED_WVTBL(basepaths, varargin)
%
%   SUMMARY:
%       The waveform twin of ed_tbl. Where ed_tbl answers "how many, and
%       where", this answers "what did they look like" - one row per ACCEPTED
%       event, carrying its LFP snippet alongside the mouse and the state, so
%       guiTbl_xy can tile by mouse or genotype and overlay a median trace.
%       Its legend then reports the event count per group, which is the second
%       number worth having.
%
%       Waveforms are DETRENDED per event and not normalised. Every event sits
%       on its own slow deflection, and averaging raw snippets averages those
%       offsets too - the mean sags and its flanks follow whatever the drifts
%       happened to do. Amplitude is left alone because it is usually the thing
%       the average is reporting (see evt_detrend).
%
%       They are also ALIGNED (evt_align). Detection centres on the largest
%       ABSOLUTE excursion, so in a biphasic discharge some rows sit on their
%       positive peak and others on their negative trough - and the average of
%       that mix is smeared with a notch at t = 0. Fixing every row to the same
%       feature lines them up; it is a shift within the snippet, not a re-read.
%
%       The polarity is chosen PER SESSION, by default from the average
%       waveform ('auto'): a mouse whose mean discharge swings negative aligns
%       to the trough, one that swings positive to the peak. Polarity depends
%       on which layer the electrode sits in, so it is a per-mouse property and
%       this reads it off the data rather than assuming it. Force it with
%       'align', 'trough' | 'peak', or 'none' to keep detection's alignment.
%
%       A session with no accepted events contributes no rows, which is the
%       honest answer for a mouse that has no discharges - unlike ed_tbl, where
%       a zero-count row is exactly what the rate model needs.
%
%   INPUTS:
%       basepaths - (Cell) Session directories. {pwd}
%       varargin  - Parameter/Value:
%           'basenames' - (Cell) File stems. {folder names}
%           'detrend'   - (Char) evt_detrend method, or 'none'. {'edge'}
%           'align'     - (Char) evt_align feature: 'auto' | 'trough' |
%                                'peak' | 'extremum' | 'none'. {'auto'}
%           'alignWin'  - (Num)  evt_align search half-window (s). {0.01}
%
%   OUTPUTS:
%       tbl     - (Table) One row per accepted event:
%           .lfp      - (Num) [1 x nSamp] the waveform, as a matrix column.
%           .sbjID    - (Cat) mouse id.
%           .basename - (Cat) session stem.
%           .state    - (Cat) vigilance state at the peak, merged as in
%                       ed_tbl ('unscored' if the peak is in no bout).
%       tstamps - (Vec) [1 x nSamp] window time base (s), shared by all rows.
%
%   DEPENDENCIES:
%       evt_files, evt_detrend, evt_align, evt_stateMerge, get_mname.
%
%   HISTORY:
%       260721 created for the per-mouse waveform view in mcu_ed.
%       260722 alignment (evt_align), so the per-mouse averages stop smearing
%              where detection centred some events on the peak. Polarity is
%              'auto' per session - trough for a negative mouse, peak for a
%              positive one - since it depends on the recording layer.

p = inputParser;
addRequired(p, 'basepaths', @iscell);
addParameter(p, 'basenames', {}, @iscell);
addParameter(p, 'detrend', 'edge', @ischar);
addParameter(p, 'align', 'auto', @ischar);
addParameter(p, 'alignWin', 0.01, @isnumeric);
parse(p, basepaths, varargin{:});
basenames = p.Results.basenames;
flgDt     = p.Results.detrend;
flgAlign  = p.Results.align;
alignWin  = p.Results.alignWin;

if isempty(basenames)
    [~, basenames] = cellfun(@fileparts, basepaths, 'uni', false);
end
mnames = get_mname(basepaths);

tstamps = [];
rows = cell(numel(basepaths), 1);
for iPath = 1 : numel(basepaths)

    files = evt_files(basepaths{iPath}, basenames{iPath}, 'ed');
    if ~isfile(files.evt) || ~isfile(files.maps)
        warning('ed_wvTbl:noFile', 'no ed / edMaps for %s; skipping', ...
            basenames{iPath});
        continue
    end
    S = load(files.evt, 'ed');
    M = load(files.maps, 'edMaps');
    ed = S.ed;

    if isempty(tstamps)
        tstamps = M.edMaps.tstamps;
    elseif ~isequal(numel(M.edMaps.tstamps), numel(tstamps))
        % one table means one x axis; a session cut at another rate or width
        % would have to be resampled, and doing that silently would change the
        % waveforms being compared
        error('ed_wvTbl:tstamps', ...
            '%s has a %d-sample map, the others have %d', ...
            basenames{iPath}, numel(M.edMaps.tstamps), numel(tstamps));
    end

    acc = logical(ed.accepted(:));
    if ~any(acc), continue; end

    % detrend first, then align: evt_align finds the trough on a baseline it
    % can trust, and the shift does not disturb the flanks the detrend used
    wv = evt_detrend(double(M.edMaps.lfp(acc, :)), tstamps, flgDt);
    wv = evt_align(wv, tstamps, flgAlign, alignWin);

    sbj = mnames{iPath};
    if ~contains(basenames{iPath}, sbj), sbj = basenames{iPath}; end

    rows{iPath} = table(wv, ...
        repmat(categorical({sbj}), nnz(acc), 1), ...
        repmat(categorical(basenames(iPath)), nnz(acc), 1), ...
        evtState(ed, acc), ...
        'VariableNames', {'lfp', 'sbjID', 'basename', 'state'});
end

tbl = vertcat(rows{:});

end     % EOF


% =========================================================================
%  LOCAL
% =========================================================================
function s = evtState(ed, acc)
% The state at each accepted peak, with <undefined> promoted to its own level
% so those events get a group instead of vanishing from a categorical filter.
if ~isfield(ed, 'state') || isempty(ed.state)
    s = categorical(repmat({'unscored'}, nnz(acc), 1));
    return;
end
s = removecats(ed.state(acc));
s = s(:);
if any(isundefined(s))
    s = addcats(s, {'unscored'});
    s(isundefined(s)) = 'unscored';
end
s = evt_stateMerge(s);          % QWAKE -> WAKE, LSLEEP -> NREM

end     % evtState
