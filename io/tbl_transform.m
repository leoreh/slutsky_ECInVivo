function tblOut = tbl_transform(tbl, varargin)
% TBL_TRANSFORM Applies log transformations, z-scoring, and normalization to table columns.
%
% SUMMARY:
% This function transforms table columns to improve their distributional
% properties for statistical analysis. It applies log transformation to
% highly skewed variables, z-scoring to standardize variables for
% regression analysis, and normalization relative to reference categories.
%
% INPUT (Required):
%   tbl         - Input table to be transformed.
%
% INPUT (Optional Key-Value Pairs):
%   varsExc     - Cell array of variable names to exclude from transformations {[]}.
%   varsInc     - Cell array of variable names to include in transformations {[]}.
%                  If provided, only these variables will be transformed (overrides varsExc).
%   flgZ        - Logical flag to apply z-scoring {true}.
%   flgLog      - Logical flag to apply log transformation to skewed variables {true}.
%   flgNorm     - Logical flag to apply normalization relative to reference category {false}.
%   skewThr     - Skewness threshold for log transformation {2}.
%   varsGrp     - Cell array of categorical variable names defining groups for
%                  separate transformation {[]}. If provided, transformations are
%                  applied separately within each group combination.
%   varNorm     - String. The name of the categorical variable whose reference
%                  category mean will be used for normalization (required if flgNorm=true).
%
% OUTPUT:
%   tblOut      - Transformed table with the same structure as input.
%
% EXAMPLE:
%   tblOut = tbl_transform(tbl, 'varsExc', {'UnitID', 'Group'}, 'flgZ', true);
%   tblOut = tbl_transform(tbl, 'varsGrp', {'Group', 'State'}, 'flgZ', true);
%   tblOut = tbl_transform(tbl, 'flgNorm', true, 'varNorm', 'Day', 'varsGrp', {'Group', 'State'});
%
% DEPENDENCIES:
%   None

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING & INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addRequired(p, 'tbl', @istable);
addParameter(p, 'varsExc', [], @(x) iscell(x) || isempty(x));
addParameter(p, 'varsInc', [], @(x) iscell(x) || isempty(x));
addParameter(p, 'varsGrp', [], @(x) iscell(x) || isempty(x));
addParameter(p, 'varNorm', '', @ischar);
addParameter(p, 'flgZ', false, @islogical);
addParameter(p, 'flgLog', false, @islogical);
addParameter(p, 'flgNorm', false, @islogical);
addParameter(p, 'skewThr', 2, @(x) isnumeric(x) && isscalar(x) && x > 0);

parse(p, tbl, varargin{:});

varsExc = p.Results.varsExc;
varsInc = p.Results.varsInc;
varsGrp = p.Results.varsGrp;
varNorm = p.Results.varNorm;
flgZ = p.Results.flgZ;
flgLog = p.Results.flgLog;
flgNorm = p.Results.flgNorm;
skewThr = p.Results.skewThr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get all variable names from the table
tblVars = tbl.Properties.VariableNames;
tblOut = tbl;

% Determine which variables to process
if ~isempty(varsInc)
    % Include only specified variables (takes precedence over varsExc)
    processVars = varsInc;
elseif ~isempty(varsExc)
    % Exclude specified variables
    excludeIdx = ismember(tblVars, varsExc);
    processVars = tblVars(~excludeIdx);
else
    % Process all variables
    processVars = tblVars;
end

% Filter to only numeric variables
numericIdx = cellfun(@(x) isnumeric(tbl.(x)), processVars);
processVars = processVars(numericIdx);

% Get reference category (assuming first category is reference)
if flgNorm
    catRef = categories(tblOut.(varNorm));
    catRef = catRef{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GROUP INDICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check if group-based transformation is requested
if ~isempty(varsGrp)
    % Ensure varsGrp is cellstr
    if isstring(varsGrp)
        varsGrp = cellstr(varsGrp);
    end

    % Find unique combinations of grouping variables
    uGrps = unique(tblOut(:, varsGrp), 'rows');

    % Find all group indices in advance and store in cell array
    idxGrps = cell(height(uGrps), 1);
    for iGrp = 1:height(uGrps)
        uRow = uGrps(iGrp, :);

        % Find rows matching the current unique group combination
        idxGrp = true(height(tblOut), 1);
        for iVar = 1:length(varsGrp)
            idxGrp = idxGrp & (tblOut.(varsGrp{iVar}) == uRow.(varsGrp{iVar}));
        end
        idxGrps{iGrp} = idxGrp;
    end
else
    % Single group (all data)
    idxGrps = {true(height(tblOut), 1)};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% APPLY TRANSFORMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iVar = 1:length(processVars)
    varName = processVars{iVar};

    % Apply transformations to each group
    for iGrp = 1:length(idxGrps)
        idxGrp = idxGrps{iGrp};
        varData = tblOut.(varName)(idxGrp);

        if flgLog
            % Check for non-negativity before log transform
            if min(varData) >= 0
                s = skewness(varData);
                if s > skewThr
                    % Apply log10 transformation with offset only to zero values
                    offset = min(varData(varData > 0)) / 2;
                    varData(varData == 0) = offset;
                    varData = log10(varData);
                    fprintf('Log-transforming variable: %s (skewness=%.2f)\n', varName, s);
                end
            end
        end

        % Apply normalization (after log transform, before z-scoring)
        if flgNorm
            % Find rows within this group that match the reference category of varNorm
            idxRef = tblOut.(varNorm) == catRef;
            idxRefGroup = idxGrp & idxRef;

            % Calculate the mean for the reference category within this group
            refMean = mean(tblOut.(varName)(idxRefGroup), 'omitnan');

            % Check for issues with reference mean
            if isnan(refMean) || refMean == 0
                if refMean == 0
                    refMean = eps; % Avoid division by zero, result will be large
                end
            end

            % Calculate percentage difference from reference
            % For values above reference, we want them to be >100%
            % For values below reference, we want them to be <100%
            varDiff = (varData - refMean) / abs(refMean);
            aboveRef = varData > refMean;
            normData = 100 + (aboveRef .* abs(varDiff) - ~aboveRef .* abs(varDiff)) * 100;

            % Update varData with normalized values
            varData = normData;
        end

        % Apply z-scoring
        if flgZ
            varData = zscore(varData);
        end

        tblOut.(varName)(idxGrp) = varData;
    end
end

end     % EOF