function tbl = spk_byCond(fcn, bouts, lbls, varargin)
% SPK_BYCOND Runs a per-unit metric separately on each labelled interval set.
%
%   tbl = SPK_BYCOND(fcn, bouts, lbls, varargin)
%
%   SUMMARY:
%       Maps one metric function over a list of labelled interval sets and
%       stacks the results into a long table, one row per unit x condition.
%       The metric never learns what a condition is; it only ever sees one
%       [nBouts x 2] interval set and returns [nUnits x 1] per field. The
%       condition therefore lives in labelled rows instead of an unnamed
%       trailing matrix dimension, which is what lets two conditions (e.g.
%       vigilance state and time of day) be crossed without either one
%       having to know about the other.
%
%       Conditions are arbitrary. Vigilance states (ss.bouts.times), time
%       chunks (n2chunks), drug epochs - only bouts and lbls change.
%
%   INPUTS:
%       fcn      - (Fcn)   Handle taking one [nBouts x 2] interval set and
%                          returning a struct whose first field is [nUnits
%                          x ...]. Fields with a different row count (an
%                          info struct, say) are dropped. Close over the
%                          rest, e.g. @(b) spktimes_metrics(spktimes, b).
%       bouts    - (Cell)  {nCond x 1} of [nBouts x 2] interval sets [s].
%       lbls     - (Cell)  {nCond x 1} of condition names.
%       varargin - Parameter/Value pairs:
%           'colName' - (Char) Name of the condition column. {'state'}
%
%   OUTPUTS:
%       tbl      - (Table) [nUnits * nCond x nFlds + 2]. Carries uid (unit
%                          index within the session) and the condition
%                          column, so two such tables join on
%                          {'uid', colName}. Empty conditions contribute no
%                          rows.
%
%   DEPENDENCIES:
%       None.
%
%   HISTORY:
%       Created: 260719
%
%   See also: SPK_STATES, SPKTIMES_METRICS, BURST_STATS

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'fcn', @(x) isa(x, 'function_handle'));
addRequired(p, 'bouts', @iscell);
addRequired(p, 'lbls', @iscell);
addParameter(p, 'colName', 'state', @ischar);

parse(p, fcn, bouts, lbls, varargin{:});
colName = p.Results.colName;

nCond = numel(bouts);
if numel(lbls) ~= nCond
    error('bouts and lbls must have the same number of elements')
end


%% ========================================================================
%  MAP
%  ========================================================================

tblCond = cell(nCond, 1);

for iCond = 1 : nCond

    % a state with no bouts in this recording contributes no rows rather
    % than a block of nans, so exposure stays visible downstream
    if isempty(bouts{iCond})
        continue
    end

    % metric structs in this repo carry a scalar info field alongside the
    % per-unit ones; row count comes from the first field and anything that
    % does not match it is parameters, not data
    s = fcn(bouts{iCond});
    flds = fieldnames(s);
    nUnits = size(s.(flds{1}), 1);
    s = rmfield(s, flds(cellfun(@(f) size(s.(f), 1) ~= nUnits, flds)));

    tblCond{iCond} = struct2table(s);

    tblCond{iCond}.uid = (1 : nUnits)';
    tblCond{iCond}.(colName) = repmat(categorical(lbls(iCond)), nUnits, 1);
end


%% ========================================================================
%  STACK
%  ========================================================================

tblCond = tblCond(~cellfun(@isempty, tblCond));

if isempty(tblCond)
    tbl = table();
    return
end

tbl = vertcat(tblCond{:});
tbl = movevars(tbl, {'uid', colName}, 'Before', 1);

end     % EOF
