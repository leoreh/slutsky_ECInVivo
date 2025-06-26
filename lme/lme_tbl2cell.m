function [dataGrp, xVals, varLbls] = lme_tbl2cell(tbl, varsFxd, varRsp, varargin)
% LME_TBL2CELL Organizes table data into cell arrays for LME plotting.
%
% This function takes a data table and organizes it into cell arrays suitable
% for plotting functions. For single fixed effects, it creates a data matrix.
% For multiple fixed effects, it groups data by the specified group variable and creates
% separate matrices for each group.
%
% INPUTS:
%   tbl         (Required) Table: Data table used in LME analysis.
%   varsFxd     (Required) Cell Array: Names of fixed effects variables.
%   varRsp      (Required) Char: Name of the response variable.
%
% VARARGIN (Name-Value Pairs):
%   'flgCats'   (Optional) Logical: Force categorical ordering for all variables.
%                       Default is false (uses natural ordering).
%   'grpVar'    (Optional) Char: Name of the variable to use for grouping.
%                       If not specified, uses the first variable in varsFxd.
%
% OUTPUTS:
%   dataGrp     Cell Array/Matrix: For single variable, returns data matrix.
%               For multiple variables, returns cell array of data matrices
%               grouped by the specified group variable.
%   xVals       Cell Array/Char: For single variable, returns char array of levels.
%               For multiple variables, returns cell array of level names for
%               each group.
%   varLbls     String Array: For single variable, returns empty array.
%               For multiple variables, returns group labels.
%
% EXAMPLES:
%   % Single variable case
%   [dataGrp, xVals, varLbls] = lme_tbl2cell(dataTable, {'Group'}, 'Response')
%   % Returns: dataGrp = [nLevels x nSubjects] matrix
%   %          xVals = char array of group names
%   %          varLbls = []
%
%   % Two variable case with explicit group variable
%   [dataGrp, xVals, varLbls] = lme_tbl2cell(dataTable, {'Group', 'Time'}, 'Response', 'grpVar', 'Group')
%   % Returns: dataGrp = {[nTime x nSubjects_Group1], [nTime x nSubjects_Group2], ...}
%   %          xVals = {timeLevels, timeLevels, ...}
%   %          varLbls = ['Group1', 'Group2', ...]
%
% See also: lme_plot, lme_sigLines, lme_frml2vars

% 260524 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addParameter(p, 'flgCats', false, @islogical);
addParameter(p, 'grpVar', '', @ischar);
parse(p, varargin{:});

opts = p.Results;
varRsp = char(varRsp); % Ensure char format

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VALIDATE INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if all variables exist in the table
allVars = [varsFxd, {varRsp}];
missingVars = setdiff(allVars, tbl.Properties.VariableNames);
if ~isempty(missingVars)
    error('Variables not found in table: %s', strjoin(missingVars, ', '));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETERMINE GROUPING VARIABLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(varsFxd) >= 2
    % Multiple variables case - determine which is the group variable
    if ~isempty(opts.grpVar)
        % Use specified group variable
        if ~ismember(opts.grpVar, varsFxd)
            error('Specified grpVar "%s" not found in fixed effects variables: %s', ...
                opts.grpVar, strjoin(varsFxd, ', '));
        end
        grpVar = opts.grpVar;
        xVar = setdiff(varsFxd, grpVar, 'stable');
        xVar = xVar{1};
    else
        % Default: first variable is group, second is x-axis
        grpVar = varsFxd{1};
        xVar = varsFxd{2};
    end
else
    % Single variable case
    grpVar = '';
    xVar = varsFxd{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORGANIZE DATA BASED ON NUMBER OF VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(varsFxd) >= 2
    % Case: Multiple variables - group by grpVar, x-axis is xVar

    % Get unique values maintaining order
    if opts.flgCats || iscategorical(tbl.(grpVar))
        grpVals = categories(categorical(tbl.(grpVar)));
    else
        grpVals = unique(tbl.(grpVar), 'stable');
    end

    % Initialize outputs
    nGrps = length(grpVals);
    dataGrp = cell(nGrps, 1);
    xVals = cell(nGrps, 1);
    varLbls = strings(nGrps, 1);

    % Get data for each group
    for igrp = 1:nGrps
        % Get data for current group
        idx = tbl.(grpVar) == grpVals(igrp);
        grpTbl = tbl(idx, :);

        % Get unique x values using categorical order if possible
        if opts.flgCats || iscategorical(grpTbl.(xVar))
            xUnique = categories(categorical(grpTbl.(xVar)));
        else
            xUnique = unique(grpTbl.(xVar), 'stable');
        end
        nX = length(xUnique);

        % Organize data matrix
        dataMat = nan(nX, sum(idx));
        for ix = 1:nX
            xIdx = grpTbl.(xVar) == xUnique{ix};
            yVals = grpTbl.(varRsp)(xIdx);
            dataMat(ix, 1:length(yVals)) = yVals';
        end

        % Store results
        dataGrp{igrp} = dataMat;
        xVals{igrp} = string(xUnique);
        varLbls{igrp} = char(grpVals(igrp));
    end

else
    % Case: Single variable
    var1 = varsFxd{1};

    % Get unique values using categorical order if possible
    if opts.flgCats || iscategorical(tbl.(var1))
        xUnique = categories(categorical(tbl.(var1)));
    else
        xUnique = unique(tbl.(var1), 'stable');
    end
    nX = length(xUnique);

    % Initialize data matrix
    maxPts = max(histcounts(categorical(tbl.(var1))));
    dataMat = nan(nX, maxPts);

    % Fill data matrix
    for ix = 1:nX
        idx = tbl.(var1) == xUnique{ix};
        yVals = tbl.(varRsp)(idx);
        dataMat(ix, 1:length(yVals)) = yVals';
    end

    % Prepare outputs
    dataGrp = dataMat;
    xVals = char(xUnique);
    varLbls = [];
end

end 