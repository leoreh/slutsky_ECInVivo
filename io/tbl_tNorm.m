function tblOut = tbl_tNorm(tbl, varargin)
% TBL_TNORM Normalizes vector columns in a table to a specific window.
%
% SUMMARY:
% This function transforms table columns (vectors/matrices) by normalizing
% them based on statistics calculated from a specific window. It supports
% both row-wise and group-wise normalization using MATLAB's native
% NORMALIZE function.
%
% INPUT (Required):
%   tbl         - Input table to be transformed.
%
% INPUT (Optional Key-Value Pairs):
%   varsInc     - Cell array of variable names to include in normalization {[]}.
%                  If provided, only these variables will be processed.
%                  Variables must be numeric matrices.
%   varsGrp     - Cell array of categorical variable names defining groups for
%                  separate transformation {[]}. If provided, normalization
%                  statistics (e.g., mean, std) are calculated per group.
%   winNorm     - 2-element numeric vector [Start, End] defining the
%                  normalization window indices {[]}. If empty, the entire
%                  vector is used.
%   Method      - String defining the normalization method {'percentage'}.
%                  Options:
%                  'percentage' - Divide by mean ( x / mean(win) )
%                  'zscore'     - Subtract mean, divide by std ( (x-mu)/sigma )
%                  'center'     - Subtract mean ( x - mu )
%                  'scale'      - Divide by std ( x / sigma )
%                  'range'      - Scale to range ( (x-min)/(max-min) )
%
% OUTPUT:
%   tblOut      - Transformed table with the same structure as input.
%
% EXAMPLE:
%   tbl = tbl_tNorm(tbl, 'varsInc', {'LFP'}, 'winNorm', [1, 100]);
%   tbl = tbl_tNorm(tbl, 'varsInc', {'FR'}, 'varsGrp', {'Genotype'}, 'Method', 'zscore');
%
% DEPENDENCIES:
%   None
%
%   See also: TBL_TRANSFORM, NORMALIZE

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'tbl', @istable);
addParameter(p, 'varsInc', [], @(x) iscell(x) || ischar(x) || isstring(x) || isempty(x));
addParameter(p, 'varsGrp', [], @(x) iscell(x) || ischar(x) || isstring(x) || isempty(x));
addParameter(p, 'winNorm', [], @(x) isnumeric(x) && (isempty(x) || numel(x)==2));
addParameter(p, 'Method', 'percentage', @(x) ischar(x) || isstring(x));

parse(p, tbl, varargin{:});

varsInc = p.Results.varsInc;
varsGrp = p.Results.varsGrp;
winNorm = p.Results.winNorm;
method = p.Results.Method;

% Handle unmatched arguments for normalize
unmatched = p.Unmatched;
normArgs = reshape([fieldnames(unmatched), struct2cell(unmatched)]', 1, []);

% Standardize cell arrays
if ischar(varsInc) || isstring(varsInc), varsInc = cellstr(varsInc); end
if ischar(varsGrp) || isstring(varsGrp), varsGrp = cellstr(varsGrp); end


%% ========================================================================
%  INITIALIZATION
%  ========================================================================

tblOut = tbl;
tblVars = tbl.Properties.VariableNames;

% Determine variables to process
if ~isempty(varsInc)
    processVars = varsInc;
else
    processVars = tblVars;
end

% Filter to only valid variables (present in table)
idxExist = ismember(processVars, tblVars);
if any(~idxExist)
    warning('Some variables in varsInc were not found in the table.');
    processVars = processVars(idxExist);
end

% Filter to only numeric variables
validIdx = false(size(processVars));
for iVar = 1:numel(processVars)
    v = processVars{iVar};
    val = tbl.(v);
    if isnumeric(val)
        validIdx(iVar) = true;
    end
end
processVars = processVars(validIdx);


%% ========================================================================
%  GROUP INDICES
%  ========================================================================

% Check if group-based transformation is requested
if ~isempty(varsGrp)
    % Find unique combinations of grouping variables
    uGrps = unique(tblOut(:, varsGrp), 'rows');

    % Find all group indices in advance
    idxGrps = cell(height(uGrps), 1);
    for iGrp = 1:height(uGrps)
        uRow = uGrps(iGrp, :);

        % Find rows matching the current unique group combination
        idxGrp = true(height(tblOut), 1);
        for iVar = 1:length(varsGrp)
            v = varsGrp{iVar};
            colData = tblOut.(v);
            val = uRow.(v);

            % Handle string vs numeric comparison
            if isstring(colData) || iscategorical(colData)
                % For string/categorical
                idxGrp = idxGrp & (colData == val);
            else
                % Numeric
                idxGrp = idxGrp & (colData == val);
            end
        end
        idxGrps{iGrp} = idxGrp;
    end
else
    idxGrps = {};
end


%% ========================================================================
%  TRANSFORM
%  ========================================================================

for iVar = 1:length(processVars)
    varName = processVars{iVar};
    varData = tblOut.(varName);

    % Determine Window Indices
    wS = 1; wE = inf;
    if ~isempty(winNorm)
        wS = max(1, winNorm(1));
        wE = winNorm(2);
    end

    % ---------------------------------------
    % Case 1: Row-wise Normalization (Default)
    % ---------------------------------------
    if isempty(idxGrps)

        % Vectorized Matrix Operation
        nC = size(varData, 2);
        currE = min(nC, wE);

        if currE >= wS
            winMat = varData(:, wS:currE);

            % Calculate Row-wise Stats
            mu = mean(winMat, 2, 'omitnan');
            sig = std(winMat, 0, 2, 'omitnan');
            mn = min(winMat, [], 2, 'omitnan');
            mx = max(winMat, [], 2, 'omitnan');

            [C, S] = getNormParams(method, mu, sig, mx, mn);

            varData = normalize(varData, 2, 'center', C, 'scale', S, normArgs{:});
        end

        tblOut.(varName) = varData;

    else
        % ---------------------------------------
        % Case 2: Group-wise Normalization
        % ---------------------------------------

        for iGrp = 1:length(idxGrps)
            idxGrp = idxGrps{iGrp};
            grpData = varData(idxGrp, :);

            % Extract Pooled Window Data
            nC = size(grpData, 2);
            currE = min(nC, wE);
            allWinVals = [];

            if currE >= wS
                subMat = grpData(:, wS:currE);
                allWinVals = subMat(:);
            end

            % Calculate Group Stats
            [mu, sig, mx, mn] = getStats(allWinVals);

            % Determine Normalization Params
            [C, S] = getNormParams(method, mu, sig, mx, mn);

            % Apply Normalize to Group Rows
            varData(idxGrp, :) = normalize(varData(idxGrp, :), 2, 'center', C, 'scale', S, normArgs{:});
        end

        tblOut.(varName) = varData;
    end
end

end

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function [mu, sig, mx, mn] = getStats(vals)
% Calculates basic stats ignoring NaNs
vals = vals(~isnan(vals));
if isempty(vals)
    mu = NaN; sig = NaN; mx = NaN; mn = NaN;
else
    mu = mean(vals);
    sig = std(vals);
    mx = max(vals);
    mn = min(vals);
end
end

function [C, S] = getNormParams(method, mu, sig, mx, mn)
% Map method to Center/Scale parameters
switch lower(method)
    case 'percentage'
        C = 0;
        S = mu;
    case 'zscore'
        C = mu;
        S = sig;
    case 'center'
        C = mu;
        S = 1;
    case 'scale'
        C = 0;
        S = sig;
    case 'range'
        C = mn;
        S = mx - mn;
    otherwise
        C = mu;
        S = sig;
end

S(S == 0) = eps;
end
