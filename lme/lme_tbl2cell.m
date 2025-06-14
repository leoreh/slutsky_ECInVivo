function [dataGrp, xVals, varLbls] = lme_tbl2cell(tbl, varsFxd, varRsp, varargin)
% LME_TBL2CELL Organizes table data into cell arrays for LME plotting.
%
% This function takes a data table and organizes it into cell arrays suitable
% for plotting functions. For single fixed effects, it creates a data matrix.
% For multiple fixed effects, it groups data by the first variable and creates
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
%
% OUTPUTS:
%   dataGrp     Cell Array/Matrix: For single variable, returns data matrix.
%               For multiple variables, returns cell array of data matrices
%               grouped by first variable.
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
%   % Two variable case
%   [dataGrp, xVals, varLbls] = lme_tbl2cell(dataTable, {'Group', 'Time'}, 'Response')
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
% ORGANIZE DATA BASED ON NUMBER OF VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(varsFxd) >= 2
    % Case: Multiple variables - group by first var, x-axis is second var
    var1 = varsFxd{1};
    var2 = varsFxd{2};

    % Get unique values maintaining order
    if opts.flgCats || iscategorical(tbl.(var1))
        var1Vals = categories(categorical(tbl.(var1)));
    else
        var1Vals = unique(tbl.(var1), 'stable');
    end

    % Initialize outputs
    nGrps = length(var1Vals);
    dataGrp = cell(nGrps, 1);
    xVals = cell(nGrps, 1);
    varLbls = strings(nGrps, 1);

    % Get data for each group
    for igrp = 1:nGrps
        % Get data for current group
        idx = tbl.(var1) == var1Vals(igrp);
        grpTbl = tbl(idx, :);

        % Get unique x values using categorical order if possible
        if opts.flgCats || iscategorical(grpTbl.(var2))
            xUnique = categories(categorical(grpTbl.(var2)));
        else
            xUnique = unique(grpTbl.(var2), 'stable');
        end
        nX = length(xUnique);

        % Organize data matrix
        dataMat = nan(nX, sum(idx));
        for ix = 1:nX
            xIdx = grpTbl.(var2) == xUnique{ix};
            yVals = grpTbl.(varRsp)(xIdx);
            dataMat(ix, 1:length(yVals)) = yVals';
        end

        % Store results
        dataGrp{igrp} = dataMat;
        xVals{igrp} = string(xUnique);
        varLbls{igrp} = char(var1Vals(igrp));
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