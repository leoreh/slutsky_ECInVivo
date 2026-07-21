function tbl = ed_tbl(basepaths, basenames)
% ED_TBL Discharge counts and rates per session, by vigilance state.
%
%   tbl = ED_TBL(basepaths, basenames)
%
%   SUMMARY:
%       Assembles the ED pipeline's answer across sessions: how many discharges
%       a mouse has, and how they distribute over vigilance states. One row per
%       session x state, following the repo's state-conditioning convention - a
%       state is a table ROW, never a trailing matrix dimension - so the result
%       drops straight into an LME with state fixed and sbjID grouping.
%
%       Counts are over ACCEPTED events, and exposure is the summed duration of
%       that state's scored bouts, so a rate is events per MINUTE OF THAT STATE.
%       A state that was scored but held no discharge still gets a row with
%       nEd = 0, which is the row a rate model needs and the one a naive count
%       silently drops.
%
%       The mouse is the unit. One row per session x state means each mouse
%       contributes one number per state, which is the level genotype was
%       assigned at - pooling events across mice would weight a long recording
%       or a busy mouse as if it were more animals.
%
%       States are MERGED for reporting (evt_stateMerge): QWAKE joins WAKE and
%       LSLEEP joins NREM. The scoring on disk is untouched. Counts and
%       exposures are pooled before the division, so a merged rate is the
%       pooled count over the pooled time.
%
%       Each session also gets a row with state 'ALL': every accepted event over
%       the whole analysis window. It is NOT the sum of the state rows - an
%       event whose peak falls in an unscored gap belongs to no state - and the
%       gap between the two is a direct read on how completely a session was
%       scored.
%
%   INPUTS:
%       basepaths - (Cell) Session directories. {pwd}
%       basenames - (Cell) File stems, for the flat layout where several
%                          recordings share one folder (the EA cohort) and the
%                          stem does not follow it. {folder names}
%
%   OUTPUT:
%       tbl - (Table) One row per session x state:
%           .sbjID    - (Cat) mouse id (get_mname).
%           .basename - (Cat) session stem.
%           .state    - (Cat) vigilance state (merged), or 'ALL'.
%           .nEd      - (Num) accepted discharges in that state.
%           .durState - (Num) scored duration of that state [min].
%           .edRate   - (Num) nEd / durState [events / min].
%
%   DEPENDENCIES:
%       basepaths2vars, evt_boutTimes, evt_files, evt_stateMerge,
%       get_mname.
%
%   HISTORY:
%       260720 created with the staged ED pipeline; the cross-session product
%              the pipeline exists to produce.
%       260721 rates moved from per hour to per MINUTE, with durState to match.
%              Discharges are counted in dozens per 24 h, so per-hour numbers
%              are small and per-minute ones smaller - but one unit throughout
%              beats two, and per minute is what the figures report.
%       260721b QWAKE and LSLEEP merged into WAKE and NREM for reporting.

if nargin < 1 || isempty(basepaths), basepaths = {pwd}; end
if nargin < 2 || isempty(basenames)
    [~, basenames] = cellfun(@fileparts, basepaths, 'uni', false);
end

mnames = get_mname(basepaths);
rows   = cell(numel(basepaths), 1);

for iPath = 1 : numel(basepaths)

    basepath = basepaths{iPath};
    basename = basenames{iPath};
    files = evt_files(basepath, basename, 'ed');
    if ~isfile(files.evt)
        warning('ed_tbl:noFile', 'no %s.ed.mat; skipping', basename);
        continue
    end

    S   = load(files.evt, 'ed');
    ed  = S.ed;
    acc = logical(ed.accepted(:));

    % get_mname reads the mouse out of the PATH, which only works when the
    % session has its own folder. In the flat layout the stem is the recording
    % id, so it stands in.
    sbj = mnames{iPath};
    if ~contains(basename, sbj), sbj = basename; end

    % Bouts are used only for their DURATIONS, which are invariant to the
    % window shift - events are matched to states through the stored ed.state
    % label, never by comparing times across frames.
    win = ed.info.win;
    if isinf(win(2)), winDur = Inf; else, winDur = win(2) - win(1); end
    v = basepaths2vars('basepaths', {basepath}, 'vars', {'sleep_states'});
    boutTimes = evt_boutTimes(v, win, winDur);

    nSt   = numel(boutTimes);           % 0 on a session that was never scored
    state = strings(nSt + 1, 1);
    nEd   = zeros(nSt + 1, 1);
    durSt = zeros(nSt + 1, 1);

    for iSt = 1 : nSt
        state(iSt) = v.ss.info.names{iSt};
        durSt(iSt) = sum(diff(boutTimes{iSt}, 1, 2)) / 60;
        nEd(iSt)   = sum(acc & ed.state == state(iSt));
    end
    state(end) = "ALL";
    nEd(end)   = sum(acc);
    durSt(end) = ed.info.sigDur / 60;

    % QWAKE -> WAKE, LSLEEP -> NREM. Counts and exposures are summed BEFORE the
    % division: a merged rate is the pooled count over the pooled time, not the
    % mean of the two rates, which would weight 20 min of LSLEEP as heavily as
    % 8 h of NREM.
    state = evt_stateMerge(state);
    [state, ~, iSt] = unique(state, 'stable');
    nEd   = accumarray(iSt, nEd);
    durSt = accumarray(iSt, durSt);

    keep = durSt > 0;                   % a state never scored gets no row
    rows{iPath} = table( ...
        repmat(categorical({sbj}), sum(keep), 1), ...
        repmat(categorical({basename}), sum(keep), 1), ...
        categorical(state(keep)), nEd(keep), durSt(keep), ...
        nEd(keep) ./ durSt(keep), 'VariableNames', ...
        {'sbjID', 'basename', 'state', 'nEd', 'durState', 'edRate'});
end

tbl = vertcat(rows{:});

end     % EOF
