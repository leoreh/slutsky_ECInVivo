function [accepted, breakdown] = evt_gate(evt, qa)
% EVT_GATE Apply a QA filter spec to detected events -> accepted mask.
%
%   [accepted, breakdown] = EVT_GATE(evt, qa)
%
%   SUMMARY:
%       The QA engine shared by both event pipelines: turns a filter spec -
%       which vigilance states to keep and per-metric [lo hi] ranges - into a
%       per-event logical mask. Marks, never removes: every detected event keeps
%       its row; accepted just selects. This is the single "filters -> mask"
%       implementation. ripp_curate and ed_curate run it headless as the
%       automatic gate and live as their GUIs recompute the kept/removed split;
%       ripp_screen applies it per method.
%
%       Permissive by construction, matching the pipelines' intent: an empty
%       state list keeps every state; a NaN metric (criterion unavailable - e.g.
%       no MUA gives NaN spkGain) passes rather than rejects. An event is
%       accepted iff its state is kept AND every metric is in range (or NaN).
%
%       Note on under-scored sessions: the state filter is categorical
%       membership on evt.state, so events with an <undefined> state (peak
%       outside every scored bout) are NOT kept by a non-empty state list. On a
%       session with partial or absent sleep scoring, either set qa.states = []
%       to gate on the metrics alone, or set qa.unscored = true to keep the
%       unlabelled events explicitly. A session with no .state field at all
%       (never state-labelled) passes the state criterion trivially.
%
%   INPUTS:
%       evt  - <struct> detected events with per-event fields aligned to ALL
%                       events: needs .peakTime, .state (for a state filter),
%                       and each metric named in qa.ranges (e.g. .emg .amp).
%       qa   - <struct> the filter spec (see ripp_methods / ed_methods '.qa'):
%           .states   - <vec|cellstr> AccuSleep state indices, or labels, to
%                              keep (1=WAKE 2=QWAKE 3=LSLEEP 4=NREM); [] = any.
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
%       260719 created as ripp_gate, the ripple QA engine.
%       260720 moved to lfp/events as evt_gate and shared with the ED pipeline
%              (a pure rename; the body was already event-agnostic). The .state
%              read gained an isfield guard, because an ED session may be
%              detected on a recording that was never sleep-scored.

N = numel(evt.peakTime);
accepted = true(N, 1);
breakdown = struct();

% state filter (categorical membership; an empty list keeps all states).
% qa.states may be AccuSleep indices (resolved to labels) or labels directly.
% qa.unscored (optional) also keeps events with an <undefined> state (peak in no
% scored bout); without it a non-empty state list drops those events.
stateOk = true(N, 1);
hasState = isfield(evt, 'state') && ~isempty(evt.state);
if isfield(qa, 'states') && ~isempty(qa.states) && hasState
    if isnumeric(qa.states)
        keepNames = evt_stateNames(qa.states);
    else
        keepNames = cellstr(qa.states);
    end
    stateOk = ismember(evt.state(:), keepNames);
end
if isfield(qa, 'unscored') && qa.unscored && hasState
    stateOk = stateOk | isundefined(evt.state(:));
end
breakdown.state = stateOk;
accepted = accepted & stateOk;

% per-metric [lo hi] ranges (NaN passes; a missing field is a spec error)
if isfield(qa, 'ranges') && ~isempty(qa.ranges)
    flds = fieldnames(qa.ranges);
    for iFld = 1:numel(flds)
        fld = flds{iFld};
        if ~isfield(evt, fld)
            error('evt_gate:metric', ...
                'qa.ranges references missing per-event field "%s"', fld);
        end
        rng  = qa.ranges.(fld);
        vals = evt.(fld)(:);
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
function names = evt_stateNames(idx)
% Resolve AccuSleep state indices to their categorical labels (config order).
cfg = as_loadConfig([]);
allNames = cfg.names;
idx = idx(idx >= 1 & idx <= numel(allNames));
names = allNames(idx);
end
