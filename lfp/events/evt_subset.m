function evt = evt_subset(evt, mask)
% EVT_SUBSET Keep only the events selected by a logical mask.
%
%   evt = EVT_SUBSET(evt, mask)
%
%   SUMMARY:
%       Shared per-event field subsetting for both event pipelines. Applies the
%       pass mask to every per-event field of the event struct - each numeric,
%       logical, or categorical field whose first dimension equals numel(mask)
%       is kept row-wise; every other field (scalar structs such as .info,
%       char fields, coincidental sizes) is left untouched. Used after evt_qa to
%       drop the events that failed quality assurance.
%
%   INPUTS:
%       evt   - (Struct) Event struct with per-event fields ([N x k] rows).
%       mask  - (Vec)    [N x 1] logical, true = keep.
%
%   OUTPUTS:
%       evt   - (Struct) Same struct, per-event fields reduced to the kept rows.
%
%   DEPENDENCIES:
%       None.
%
%   HISTORY:
%       Created: 260706 (absorbs the duplicated wrapper QA-subset loops).

mask = logical(mask(:));
fn   = fieldnames(evt);
for iFld = 1:numel(fn)
    f = evt.(fn{iFld});
    if (isnumeric(f) || islogical(f) || iscategorical(f)) ...
            && size(f, 1) == numel(mask)
        evt.(fn{iFld}) = f(mask, :);
    end
end

end     % EOF
