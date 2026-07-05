function [idxGood, breakdown] = evt_qa(peakTime, varargin)
% EVT_QA Combine event quality-assurance criteria into one pass mask.
%
%   [idxGood, breakdown] = EVT_QA(peakTime, varargin)
%
%   SUMMARY:
%       A pure combiner shared by both event pipelines. It applies, in any
%       combination, two kinds of optional criteria and ANDs them:
%           1. Time inclusion / exclusion — keep events whose peak falls inside
%              'inTimes', drop those inside 'exTimes'.
%           2. Per-event metric thresholds — each metric column passes when its
%              value is within the matching [lo hi] range.
%       Every criterion is optional; a metric value of NaN passes (that event's
%       criterion is skipped), an empty range passes, and an empty 'inTimes'
%       keeps all events. Metrics are computed by the caller (e.g. evt_emgScore,
%       evt_spkGain) and thresholds are passed in — evt_qa makes no assumptions
%       about what a metric means.
%
%   INPUTS:
%       peakTime  - (Vec)  [N x 1] Event peak times [s].
%       varargin  - Parameter/Value pairs:
%           'inTimes' - (Mat)  [M x 2] Keep events whose peak is inside these
%                              intervals [s]. Empty -> keep all. (Default: []).
%           'exTimes' - (Mat)  [K x 2] Drop events whose peak is inside these
%                              intervals [s]. Empty -> drop none. (Default: []).
%           'metrics' - (Mat)  [N x j] Per-event metric values, one column each.
%           'ranges'  - (Cell) {1 x j} Each [lo hi]; the event passes column c
%                              if lo <= metric(:,c) <= hi (or the value is NaN,
%                              or the range is empty). (Default: {}).
%           'names'   - (Cell) {1 x j} Optional names for the breakdown fields.
%
%   OUTPUTS:
%       idxGood   - (Vec)    [N x 1] logical, true = passes every criterion.
%       breakdown - (Struct) Per-criterion masks: .inc, .exc, one per metric
%                            (named by 'names' or metric<c>), and .good.
%
%   DEPENDENCIES:
%       None.
%
%   HISTORY:
%       Created: 05 Jul 2026 (replaces ripp_qa + ed_reject_emg with one generic
%                combiner; metrics/thresholds are supplied by the caller).

% =========================================================================
%  ARGUMENTS
% =========================================================================
p = inputParser;
addRequired(p, 'peakTime', @isnumeric);
addParameter(p, 'inTimes', [], @isnumeric);
addParameter(p, 'exTimes', [], @isnumeric);
addParameter(p, 'metrics', [], @isnumeric);
addParameter(p, 'ranges', {}, @iscell);
addParameter(p, 'names', {}, @iscell);
parse(p, peakTime, varargin{:});

peakTime = p.Results.peakTime(:);
inTimes  = p.Results.inTimes;
exTimes  = p.Results.exTimes;
metrics  = p.Results.metrics;
ranges   = p.Results.ranges;
names    = p.Results.names;

N = numel(peakTime);

% =========================================================================
%  TIME INCLUSION / EXCLUSION
% =========================================================================
if isempty(inTimes)
    inc = true(N, 1);
else
    inc = in_any(peakTime, inTimes);
end

if isempty(exTimes)
    exc = false(N, 1);
else
    exc = in_any(peakTime, exTimes);
end

breakdown = struct('inc', inc, 'exc', exc);

% =========================================================================
%  METRIC THRESHOLDS
% =========================================================================
metPass = true(N, 1);
for c = 1:size(metrics, 2)
    vals = metrics(:, c);
    if numel(ranges) >= c && ~isempty(ranges{c})
        rng = ranges{c};
        pass = (vals >= rng(1) & vals <= rng(2)) | isnan(vals);
    else
        pass = true(N, 1);
    end
    metPass = metPass & pass;

    nm = sprintf('metric%d', c);
    if numel(names) >= c && ~isempty(names{c}), nm = names{c}; end
    breakdown.(nm) = pass;
end

% =========================================================================
%  COMBINE
% =========================================================================
idxGood = inc & ~exc & metPass;
breakdown.good = idxGood;

end     % EOF


% -------------------------------------------------------------------------
% Helper: true where t falls inside any [start end] row of iv
% -------------------------------------------------------------------------
function m = in_any(t, iv)
m = false(numel(t), 1);
for i = 1:size(iv, 1)
    m = m | (t >= iv(i, 1) & t <= iv(i, 2));
end
end
