function tblOut = tbl_transform(tbl, varargin)
% TBL_TRANSFORM Applies log/logit transformations, z-scoring, and normalization.
%
%   TBLOUT = TBL_TRANSFORM(TBL, ...) transforms table columns to improve their
%   distributional properties for statistical analysis. It applies log or
%   logit transformations to highly skewed variables, z-scoring to
%   standardize variables, and normalization relative to a reference
%   category.
%
%   INPUTS:
%       tbl         - (table) Input data table.
%       varargin    - (param/value) Optional parameters:
%
%           SELECTION:
%           'varsInc' : (cell) Variables to transform {[]}. Overrides varsExc.
%           'varsExc' : (cell) Variables to exclude from transform {[]}.
%           'varsGrp' : (cell) Categorical variables defining groups for
%                       separate transformation operations {[]}.
%
%           TRANSFORMATION:
%           'flgZ'    : (logical) Z-score variables {true}.
%           'logBase' : (mix) Transformation mode {[]}:
%                       - [] or empty: No transformation.
%                       - 'logit': Logit transform log(y/(1-y)).
%                       - 'e': Natural Log transform log(y).
%                       - (numeric): Log transform with specific base.
%           'skewThr' : (numeric) Skewness threshold for auto-log {2}.
%           'flg0'    : (logical) Add offset if zero-inflated {false}.
%
%           NORMALIZATION:
%           'varNorm' : (char) Grouping variable for normalization. If
%                       provided, normalizes values relative to the mean of
%                       the first category (Reference).
%
%           MISC:
%           'verbose' : (logical) Print actions {false}.
%
%   OUTPUTS:
%       tblOut      - (table) Transformed table.
%
%   EXAMPLES:
%       % Z-score only
%       tbl = tbl_transform(tbl, 'flgZ', true);
%
%       % Log10 transform skewed variables
%       tbl = tbl_transform(tbl, 'logBase', 10, 'flgZ', false);
%
%       % Logit transform (for proportions)
%       tbl = tbl_transform(tbl, 'logBase', 'logit');
%
%       % Normalize 'FR' relative to 'Day' baseline
%       tbl = tbl_transform(tbl, 'varNorm', 'Day', 'varsInc', {'FR'});
%
%   See also: LME_ANALYSE, LME_TBLPREP

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);

% Selection
addParameter(p, 'varsInc', [], @(x) iscell(x) || isempty(x));
addParameter(p, 'varsExc', [], @(x) iscell(x) || isempty(x));
addParameter(p, 'varsGrp', [], @(x) iscell(x) || isempty(x));

% Transformation
addParameter(p, 'flgZ', false, @islogical);
addParameter(p, 'logBase', [], @(x) isempty(x) || isnumeric(x) || ischar(x));
addParameter(p, 'skewThr', 2, @isnumeric);
addParameter(p, 'flg0', false, @islogical);

% Normalization
addParameter(p, 'varNorm', '', @ischar);

% Misc
addParameter(p, 'verbose', false, @islogical);

parse(p, tbl, varargin{:});

varsInc = p.Results.varsInc;
varsExc = p.Results.varsExc;
varsGrp = p.Results.varsGrp;
flgZ = p.Results.flgZ;
logBase = p.Results.logBase;
skewThr = p.Results.skewThr;
flg0 = p.Results.flg0;
varNorm = p.Results.varNorm;
verbose = p.Results.verbose;


%% ========================================================================
%  INITIALIZATION
%  ========================================================================

% Get all variable names
tblVars = tbl.Properties.VariableNames;
tblOut = tbl;

% Determine variables to process
if ~isempty(varsInc)
    processVars = varsInc;
elseif ~isempty(varsExc)
    excludeIdx = ismember(tblVars, varsExc);
    processVars = tblVars(~excludeIdx);
else
    processVars = tblVars;
end

% Filter to numeric only
numericIdx = cellfun(@(x) isnumeric(tbl.(x)), processVars);
processVars = processVars(numericIdx);

% Check Normalization
flgNorm = ~isempty(varNorm);
catRef = [];

if flgNorm
    if ismember(varNorm, varsGrp)
        error('varNorm "%s" cannot be included in varsGrp', varNorm);
    end
    % Get reference category (first one)
    catRef = categories(tblOut.(varNorm));
    catRef = catRef{1};
end

if flgNorm && flgZ
    warning('Z-scoring will overwrite the scale set by normalization.');
end


%% ========================================================================
%  GROUP INDICES
%  ========================================================================

if ~isempty(varsGrp)
    % Ensure cellstr
    if isstring(varsGrp), varsGrp = cellstr(varsGrp); end

    uGrps = unique(tblOut(:, varsGrp), 'rows');
    idxGrps = cell(height(uGrps), 1);

    for iGrp = 1:height(uGrps)
        uRow = uGrps(iGrp, :);
        idxGrp = true(height(tblOut), 1);
        for iVar = 1:length(varsGrp)
            idxGrp = idxGrp & (tblOut.(varsGrp{iVar}) == uRow.(varsGrp{iVar}));
        end
        idxGrps{iGrp} = idxGrp;
    end
else
    idxGrps = {true(height(tblOut), 1)};
end


%% ========================================================================
%  TRANSFORM
%  ========================================================================

for iVar = 1:length(processVars)
    varName = processVars{iVar};
    varData = tblOut.(varName);

    % --- Global Analysis (Log/Logit/Offset) ---

    % Determine Transformation Mode
    isLogit = ischar(logBase) && strcmpi(logBase, 'logit');
    isLog = (isnumeric(logBase) && ~isempty(logBase)) || ...
        (ischar(logBase) && strcmpi(logBase, 'e'));

    % Logit Transform
    if isLogit
        varData = max(eps, min(1 - eps, varData));
        varData = log(varData ./ (1 - varData));
        if verbose
            fprintf('[%s] Applying Logit\n', varName);
        end
    end

    % Offset (Zero-Inflation)
    % Add offset if requested OR if we are about to log-transform
    if flg0 || isLog
        if any(varData(:) == 0) && all(varData(:) >= 0)
            c = min(varData(varData > 0)) / 2;
            varData = varData + c;
            if verbose
                fprintf('[%s] Adding Offset %.4f\n', varName, c);
            end
        end
    end

    % Log Transform (Conditional on Skewness)
    if isLog
        % Check skewness on pooled data
        s = skewness(varData(~isnan(varData)));

        if s > skewThr
            if ischar(logBase) && strcmpi(logBase, 'e')
                varData = log(varData);
                baseStr = 'Ln';
            else
                varData = log(varData) ./ log(logBase);
                baseStr = sprintf('Log%.1f', logBase);
            end
            if verbose
                fprintf('[%s] Applying %s (Skew=%.2f)\n', varName, baseStr, s);
            end
        end
    end

    % Write back global changes
    tblOut.(varName) = varData;

    % --- Group-wise Operations (Norm & Z) ---

    for iGrp = 1:length(idxGrps)
        idxGrp = idxGrps{iGrp};
        dataGrp = tblOut.(varName)(idxGrp);

        % Normalization
        if flgNorm
            idxRef = tblOut.(varNorm) == catRef;
            idxRefGroup = idxGrp & idxRef;

            refMean = mean(tblOut.(varName)(idxRefGroup), 'all', 'omitnan');

            if refMean == 0, refMean = eps; end

            if isnan(refMean)
                warning('[%s] Group %d RefMean is NaN. Skipping norm.\n', varName, iGrp);
            else
                % Use "Number Line" normalization
                varDiff = (dataGrp - refMean) / abs(refMean);
                dataGrp = 100 + varDiff * 100;
            end
        end

        % Z-Scoring
        if flgZ
            dataGrp = (dataGrp - mean(dataGrp, 'omitnan')) ./ ...
                std(dataGrp, 'omitnan');
        end

        % Update Table
        tblOut.(varName)(idxGrp) = dataGrp;
    end
end

end     % EOF
