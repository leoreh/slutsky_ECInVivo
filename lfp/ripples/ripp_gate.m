function [accepted, breakdown] = ripp_gate(ripp, qa)
% RIPP_GATE Apply a QA filter spec to detected ripples -> accepted mask.
%
%   [accepted, breakdown] = RIPP_GATE(ripp, qa)
%
%   SUMMARY:
%       The ripple QA engine: turns a filter spec - which vigilance states to
%       keep and per-metric [lo hi] ranges - into a per-event logical mask. Marks,
%       never removes: every detected event keeps its row; accepted just selects.
%       This is the single "filters -> mask" implementation. ripp_curate runs it
%       headless as the automatic gate and live as its GUI recomputes the
%       kept/removed split; ripp_screen applies it per method. It replaces the old
%       evt_qa call in the ripple path.
%
%       Permissive by construction, matching the pipeline's intent: an empty state
%       list keeps every state; a NaN metric (criterion unavailable - e.g. no MUA
%       gives NaN spkGain) passes rather than rejects. An event is accepted iff
%       its state is kept AND every metric is in range (or NaN).
%
%       Note on under-scored sessions: the state filter is categorical membership
%       on ripp.state, so events with an <undefined> state (peak outside every
%       scored bout) are NOT kept by a non-empty state list. On a session with
%       partial/absent sleep scoring, set qa.states = [] to gate on the metrics
%       alone rather than silently keeping unlabelled events.
%
%   INPUTS:
%       ripp - <struct> detected events with per-event fields aligned to ALL
%                       events: needs .peakTime, .state (for a state filter), and
%                       each metric named in qa.ranges (e.g. .emg .spkGain).
%       qa   - <struct> the filter spec (see ripp_methods '.qa'):
%           .states   - <vec>  AccuSleep state indices to keep (1=WAKE 2=QWAKE
%                              3=LSLEEP 4=NREM); [] = any state.
%           .unscored - <log>  also keep <undefined>-state events (peak in no
%                              scored bout). {false}
%           .ranges   - <struct> metric field -> [lo hi]; empty/absent = none.
%
%   OUTPUTS:
%       accepted  - <log> [N x 1] mask over all detected events.
%       breakdown - <struct> per-criterion masks for diagnostics: .state, one
%                            field per metric, and .accepted (their AND).
%
%   DEPENDENCIES:
%       as_loadConfig (resolves state indices to labels).
%
%   HISTORY:
%       260719 the ripple QA engine; absorbs evt_qa's role in the ripple path.

N = numel(ripp.peakTime);
accepted = true(N, 1);
breakdown = struct();

% state filter (categorical membership; an empty list keeps all states).
% qa.states may be AccuSleep indices (resolved to labels) or labels directly.
% qa.unscored (optional) also keeps events with an <undefined> state (peak in no
% scored bout); without it a non-empty state list drops those events.
stateOk = true(N, 1);
if isfield(qa, 'states') && ~isempty(qa.states)
    if isnumeric(qa.states)
        keepNames = ripp_stateNames(qa.states);
    else
        keepNames = cellstr(qa.states);
    end
    stateOk = ismember(ripp.state(:), keepNames);
end
if isfield(qa, 'unscored') && qa.unscored
    stateOk = stateOk | isundefined(ripp.state(:));
end
breakdown.state = stateOk;
accepted = accepted & stateOk;

% per-metric [lo hi] ranges (NaN passes; a missing field is a spec error)
if isfield(qa, 'ranges') && ~isempty(qa.ranges)
    flds = fieldnames(qa.ranges);
    for iFld = 1:numel(flds)
        fld = flds{iFld};
        if ~isfield(ripp, fld)
            error('ripp_gate:metric', ...
                'qa.ranges references missing per-event field "%s"', fld);
        end
        rng  = qa.ranges.(fld);
        vals = ripp.(fld)(:);
        metOk = isnan(vals) | (vals >= rng(1) & vals <= rng(2));
        breakdown.(fld) = metOk;
        accepted = accepted & metOk;
    end
end

breakdown.accepted = accepted;

end     % EOF


% =========================================================================
%  LOCAL
% =========================================================================
function names = ripp_stateNames(idx)
% Resolve AccuSleep state indices to their categorical labels (config order).
cfg = as_loadConfig([]);
allNames = cfg.names;
idx = idx(idx >= 1 & idx <= numel(allNames));
names = allNames(idx);
end
