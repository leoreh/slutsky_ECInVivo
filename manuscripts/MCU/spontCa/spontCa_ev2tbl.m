function evTbl = spontCa_ev2tbl(ev, compartment)
% SPONTCA_EV2TBL  Convert an event struct (output of spontCa_detect) to
% a per-event table.
%
% Each row of the returned table is one event; columns are
%   {compartment, start, stop, amp, dur, int}
% Compartment is a categorical with levels {Cyto, Mito}. Times in s,
% amp in dF/F, int in dF/F * s.
%
% USAGE
%   evTbl = spontCa_ev2tbl(ev, 'Cyto')
%   evTbl = [spontCa_ev2tbl(evCyto, 'Cyto'); spontCa_ev2tbl(evMito, 'Mito')]
%
% INPUTS
%   ev          - struct with column-vector fields .start .stop .amp
%                 .dur .int (all same length, possibly empty)
%   compartment - char or string, one of {'Cyto','Mito'}.
%
% See also: SPONTCA_TBL2EV, SPONTCA_DETECT, SPONTCA_FINALIZE

cmpLevels = categorical({'Cyto', 'Mito'});
cmpVal    = categorical({char(compartment)}, categories(cmpLevels));
nE        = numel(ev.start);

if nE == 0
    evTbl = table( ...
        categorical(strings(0, 1), categories(cmpLevels)), ...
        zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
        zeros(0, 1), zeros(0, 1), ...
        'VariableNames', ...
        {'compartment', 'start', 'stop', 'amp', 'dur', 'int'});
    return;
end

evTbl = table( ...
    repmat(cmpVal, nE, 1), ...
    ev.start(:), ev.stop(:), ev.amp(:), ev.dur(:), ev.int(:), ...
    'VariableNames', {'compartment', 'start', 'stop', 'amp', 'dur', 'int'});
end
