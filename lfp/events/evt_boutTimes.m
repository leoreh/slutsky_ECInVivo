function [boutTimes, vldTimes, nremTimes] = evt_boutTimes(v, win, sigDur)
% EVT_BOUTTIMES Prepare window-relative vigilance-state bout times.
%
%   [boutTimes, vldTimes, nremTimes] = EVT_BOUTTIMES(v, win, sigDur)
%
%   SUMMARY:
%       Shared sleep-state preparation for both event pipelines (ripples, ED).
%       Pulls the per-state bout times from the basepaths2vars struct v, shifts
%       them into the window-relative frame, and clips them to the window.
%       Derives two convenience sets from the clipped bouts: the valid-state
%       intervals for control matching (QWAKE + LSLEEP + NREM) and the NREM
%       intervals used as the ripple z-scoring baseline. Missing or partial
%       states degrade gracefully - any output that cannot be formed is returned
%       empty, which relaxes the criterion that consumes it downstream.
%
%   INPUTS:
%       v        - (Struct) basepaths2vars output. Reads the optional field
%                           .ss.bouts.times (cell, one [n x 2] matrix per state,
%                           s). Empty outputs when absent.
%       win      - (Vec)    Analysis window [start end] (s).
%       sigDur   - (Num)    Window duration (s); Inf for the full recording.
%
%   OUTPUTS:
%       boutTimes - (Cell)  Per-state bout times, window-relative and clipped.
%                           {} when states are absent.
%       vldTimes  - (Mat)   [m x 2] Valid-state intervals (states 2/3/4 stacked)
%                           for control matching. [] when < 4 states.
%       nremTimes - (Mat)   [k x 2] NREM intervals (state 4), the ripple
%                           z-scoring baseline. [] when < 4 states.
%
%   DEPENDENCIES:
%       None.
%
%   HISTORY:
%       Created: 260706 (absorbs the duplicated wrapper bout-time chains).

% Per-state bouts, shifted into the window frame and clipped to it
boutTimes = {};
if isfield(v, 'ss') && isfield(v.ss, 'bouts') && isfield(v.ss.bouts, 'times')
    boutTimes = v.ss.bouts.times;
    boutTimes = cellfun(@(x) x - win(1), boutTimes, 'uni', false);
    boutTimes = cellfun(@(x) x(x(:, 2) > 0 & x(:, 1) < sigDur, :), ...
        boutTimes, 'uni', false);
else
    warning('evt_boutTimes:noStates', ...
        'sleep_states not found; events left unlabelled.');
end

% Valid-state and NREM sets need the canonical 4-state layout; extract them
% from the already-clipped bouts so both stay inside the windowed signal.
vldTimes  = [];
nremTimes = [];
if numel(boutTimes) >= 4
    vldTimes  = vertcat(boutTimes{2}, boutTimes{3}, boutTimes{4});
    nremTimes = boutTimes{4};
end

end     % EOF
