function ev = spontCa_tbl2ev(evTbl, compartment)
% SPONTCA_TBL2EV  Inverse of SPONTCA_EV2TBL: pull one compartment's
% events from a long-format events table into the ev struct format
% (column vectors per field).
%
% USAGE
%   ev = spontCa_tbl2ev(evTbl, 'Cyto')
%
% INPUTS
%   evTbl       - table with columns {compartment, start, stop, amp,
%                 dur, int}; compartment is categorical.
%   compartment - char or string, one of {'Cyto','Mito'}.
%
% OUTPUT
%   ev          - struct with column-vector fields .start .stop .amp
%                 .dur .int. Empty fields if no rows match.
%
% See also: SPONTCA_EV2TBL, SPONTCA_DETECT

if isempty(evTbl)
    sub = evTbl;
else
    sub = evTbl(evTbl.compartment == char(compartment), :);
end
ev = struct( ...
    'start', sub.start(:), ...
    'stop',  sub.stop(:), ...
    'amp',   sub.amp(:), ...
    'dur',   sub.dur(:), ...
    'int',   sub.int(:));
end
