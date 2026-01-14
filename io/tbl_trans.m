function [tblOut, transParams] = tbl_trans(tbl, varargin)
% TBL_TRANS Applies log/logit transformations, z-scoring, and normalization.
%
%   [TBLOUT, PARAMS] = TBL_TRANS(TBL, ...) transforms table columns to
%   improve their distributional properties for statistical modeling.
%
%   MODES:
%       1. FIT MODE: Calculates statistics (Mean, SD, Skewness, Reference)
%          from TBL, applies transformations, and returns the parameters in
%          PARAMS.
%       2. APPLY MODE: Uses parameters from a previous run (passed via
%          'template') to transform TBL using those exact statistics.
%          Prevents data leakage between Training and Test sets.
%       3. INVERSE MODE: Reverses the transformations to recover original values.
%
%   INPUTS:
%       tbl         - (table) Input data table.
%       varargin    - (param/value) Optional parameters:
%
%           MODE SELECTION:
%           'template': (struct) Output PARAMS from a previous run. If
%                       provided, the function runs in APPLY or INVERSE mode.
%           'flgInv'  : (logical) If true (and template provided), reverses
%                       transformations. {false}
%
%           VARIABLE SELECTION (Fit Mode Only):
%           'varsInc' : (cell) Variables to transform {[]}. Overrides varsExc.
%           'varsExc' : (cell) Variables to exclude from transform {[]}.
%           'varsGrp' : (cell/char) Categorical variables defining groups for
%                       separate Z-scoring/Normalization {[]}.
%
%           TRANSFORMATION SETTINGS (Fit Mode Only):
%           'flgZ'    : (logical) Z-score variables {true}.
%                       NOTE: If true, 'varNorm' is ignored (Conflict).
%           'logBase' : (mix) Transformation mode {[]}:
%                       - []/empty: No log.
%                       - 'logit': Logit transform log(y/(1-y)).
%                       - 'e': Natural Log transform log(y).
%                       - (numeric): Log transform with specific base.
%           'skewThr' : (numeric) Skewness threshold for auto-log {2}.
%           'flg0'    : (logical) Force offset for zero-values {false}.
%                       Note: Offset is auto-added if Log+Zeros detected.
%
%           NORMALIZATION (Fit Mode Only):
%           'varNorm' : (char) Grouping variable for normalization.
%                       Normalizes values relative to the mean of the first
%                       category (Reference). {''}
%
%           MISC:
%           'verbose' : (logical) Print actions {false}.
%
%   OUTPUTS:
%       tblOut      - (table) Transformed table.
%       transParams - (struct) Structure containing applied settings and
%                     statistics, suitable for use as 'template'.
%
%   EXAMPLES:
%       % 1. FIT: Transform Training Data
%       [tblTrn, P] = tbl_trans(dataTrn, 'flgZ', true, 'logBase', 10);
%
%       % 2. APPLY: Transform Test Data (using Training stats)
%       tblTst = tbl_trans(dataTst, 'template', P);
%
%       % 3. INVERSE: Recover Original Data
%       tblRec = tbl_trans(tblTst, 'template', P, 'flgInv', true);
%
%   See also: LME_ANALYSE, LME_ABLATION

%% ========================================================================
%  INPUT PARSING
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);

% Mode A: Apply Template
addParameter(p, 'template', [], @(x) isempty(x) || isstruct(x));

% Mode B: Fit New Parameters
addParameter(p, 'varsInc', [], @(x) iscell(x) || isempty(x));
addParameter(p, 'varsExc', [], @(x) iscell(x) || isempty(x));
addParameter(p, 'varsGrp', [], @(x) iscell(x) || ischar(x) || isstring(x) || isempty(x));
addParameter(p, 'flgZ', false, @islogical);
addParameter(p, 'logBase', [], @(x) isempty(x) || isnumeric(x) || ischar(x));
addParameter(p, 'skewThr', 2, @isnumeric);
addParameter(p, 'flg0', false, @islogical);
addParameter(p, 'varNorm', '', @ischar);
addParameter(p, 'verbose', false, @islogical);

% Mode C: Inverse Transform (Back to Original)
addParameter(p, 'flgInv', false, @islogical);

parse(p, tbl, varargin{:});

tmpl    = p.Results.template;
flgInv  = p.Results.flgInv;
verbose = p.Results.verbose;
varsInc = p.Results.varsInc;
varsExc = p.Results.varsExc;


%% ========================================================================
%  CENTRALIZED VARIABLE SELECTION (procVars)
%  ========================================================================

tblVars = tbl.Properties.VariableNames;

if isempty(tmpl)
    % --- FIT MODE CANDIDATES ---
    % Candidates are all numeric, non-categorical variables in the table.
    isNum = varfun(@(x) isnumeric(x) & ~iscategorical(x), tbl, 'OutputFormat', 'uniform');
    candidates = tblVars(isNum);
else
    % --- APPLY / INVERSE MODE CANDIDATES ---
    % Candidates are variables present in both the Template and the Table.
    varsTmpl = fieldnames(tmpl.varsTrans);
    candidates = varsTmpl(ismember(varsTmpl, tblVars));
end

% --- FILTERING (varsInc / varsExc) ---
if ~isempty(varsInc)
    % Whitelist: Only process variables in varsInc (that are also candidates)
    varsTrans = candidates(ismember(candidates, varsInc));
elseif ~isempty(varsExc)
    % Blacklist: Process candidates NOT in varsExc
    varsTrans = candidates(~ismember(candidates, varsExc));
else
    % Default: Process all candidates
    varsTrans = candidates;
end

% Ensure cell array of strings
if ischar(varsTrans), varsTrans = {varsTrans}; end


%% ========================================================================
%  ROUTE LOGIC (FIT vs APPLY vs INVERSE)
%  ========================================================================

if isempty(tmpl)
    % --- FIT MODE ---
    if verbose, fprintf('[TBL_TRANS] Fit Mode \n'); end
    % Note: varsTrans is passed explicitly
    [tblOut, transParams] = trans_fit(tbl, p.Results, varsTrans);
else
    if flgInv
        % --- INVERSE MODE ---
        if verbose, fprintf('[TBL_TRANS] Inverse Mode \n'); end
        [tblOut, transParams] = trans_inv(tbl, tmpl, p.Results, varsTrans);
    else
        % --- APPLY MODE ---
        if verbose, fprintf('[TBL_TRANS] Apply Mode \n'); end
        [tblOut, transParams] = trans_apply(tbl, tmpl, p.Results, varsTrans);
    end
end

end     % EOF


%% ========================================================================
%  CORE: FIT PARAMETERS & TRANSFORM
%  ========================================================================
function [tblOut, params] = trans_fit(tbl, args, varsTrans)

% Parse Settings
varsGrp = args.varsGrp;
if ischar(varsGrp) || isstring(varsGrp)
    varsGrp = cellstr(varsGrp);
end

flgZ    = args.flgZ;
argLogBase = args.logBase; % Arg
skewThr = args.skewThr;
flg0    = args.flg0;
varNorm = args.varNorm;
verbose = args.verbose;

% Conflict Resolution: Z-Score vs Norm
if flgZ && ~isempty(varNorm)
    warning('tbl_trans:Conflict', 'Z-scoring is enabled. Normalization (varNorm) will be IGNORED.');
    varNorm = ''; % Disable Norm
end

% Normalization Ref
flgNorm = ~isempty(varNorm);
catRef = [];
if flgNorm
    if ismember(varNorm, varsGrp)
        error('varNorm "%s" cannot be included in varsGrp', varNorm);
    end
    cats = categories(tbl.(varNorm));
    catRef = cats{1}; % First category is reference
end

% Output Structure
params = struct();
params.varsTrans = struct(); % Stores individual var params
params.varsGrp = varsGrp;
params.varNorm = varNorm;
params.catRef  = catRef;
tblOut = tbl;

% Identify Groups (for Stats Calculation)
if ~isempty(varsGrp)
    [uGrps, ~, grpIdx] = unique(tblOut(:, varsGrp), 'rows');
else
    uGrps = table();
    grpIdx = ones(height(tblOut), 1);
end

% Process Each Variable
for iVar = 1:numel(varsTrans)
    var  = varsTrans{iVar};
    data = tblOut.(var);

    % --- DETERMINE LOG/OFFSET (Global) ---
    % Default: Use argument logBase
    currLogBase = argLogBase;
    doOffset    = flg0;
    c           = 0; % Offset value

    % Logic:
    % - If 'logit': Always apply
    % - If 'e' or numeric: Check Skewness. If skew > thr, KEEP currLogBase. If not, set to [].
    % - If [] (none): Do nothing.

    doLogit = ischar(currLogBase) && strcmpi(currLogBase, 'logit');
    doLog   = (isnumeric(currLogBase) && ~isempty(currLogBase)) || ...
        (ischar(currLogBase) && strcmpi(currLogBase, 'e'));

    % Logit Checks
    if doLogit
        N = length(data);
        data = (data * (N - 1) + 0.5) / N;
        data = log(data ./ (1 - data));
        if verbose, fprintf('[%s] Applied Logit.\n', var); end
    elseif doLog
        % Log Checks (Skewness)
        s = skewness(data(~isnan(data)));
        if s > skewThr
            doOffset = true; % Check offset if log enabled
            if verbose, fprintf('[%s] Skew %.2f > %.1f. Log enabled.\n', var, s, skewThr); end
        else
            % Disable log if skew not met
            currLogBase = [];
            doLog = false;
        end
    end

    % Calculate Offset (Before Log)
    if doOffset || doLog
        if any(data(:) == 0) && all(data(:) >= 0)
            c = min(data(data > 0)) / 2; % Half of min non-zero
            data = data + c;
            if verbose, fprintf('[%s] Added Offset %.4f.\n', var, c); end
        end
    end

    % Apply Log
    if doLog
        if ischar(currLogBase) && strcmpi(currLogBase, 'e')
            data = log(data);
        elseif isnumeric(currLogBase) && ~isempty(currLogBase)
            data = log(data) ./ log(currLogBase);
        end
    end

    % Update Table (Transformation applied)
    tblOut.(var) = data;

    % --- CALC GROUP STATS ---
    nGrps = max(grpIdx);
    stats = table();
    if ~isempty(varsGrp), stats = uGrps; end

    stats.Mean    = nan(nGrps, 1);
    stats.SD      = nan(nGrps, 1);
    stats.RefMean = nan(nGrps, 1);

    for iGrp = 1:nGrps
        idxK = (grpIdx == iGrp);
        vals = data(idxK);

        stats.Mean(iGrp) = mean(vals, 'omitnan');
        stats.SD(iGrp)   = std(vals, 'omitnan');

        if flgNorm
            % Find Reference Subgroup within this Group
            idxRef = (tbl.(varNorm) == catRef);
            valsRef = data(idxK & idxRef);
            stats.RefMean(iGrp) = mean(valsRef, 'omitnan');
        end
    end

    % --- APPLY GROUP TRANS (Z-Score OR Norm) ---
    for iGrp = 1:nGrps
        idxK = (grpIdx == iGrp);
        vals = tblOut.(var)(idxK);

        % Normalization
        if flgNorm
            rm = stats.RefMean(iGrp);
            if ~isnan(rm) && rm ~= 0
                vals = 100 + ((vals - rm) / abs(rm) * 100);
            elseif rm == 0
                vals = 100 + vals;
            end
        end

        % Z-Score
        if flgZ
            mu = stats.Mean(iGrp);
            sigma = stats.SD(iGrp);

            if sigma ~= 0
                vals = (vals - mu) / sigma;
            else
                vals = vals - mu;
            end
        end

        tblOut.(var)(idxK) = vals;
    end

    % Save Parameters
    pVar = struct();
    pVar.logBase  = currLogBase; % Now sole source of truth
    pVar.offset   = c;
    pVar.flgNorm  = flgNorm;
    pVar.flgZ     = flgZ;
    pVar.stats    = stats;

    params.varsTrans.(var) = pVar;
end

end


%% ========================================================================
%  CORE: APPLY TRANSFORMATION
%  ========================================================================
function [tblOut, params] = trans_apply(tbl, tmpl, args, varsTrans) %#ok<INUSD>

params = tmpl; % Stats from template
tblOut = tbl;

% Match Groups
idxGrp = match_grp(tbl, tmpl);

for iVar = 1:numel(varsTrans)
    var = varsTrans{iVar};
    pVar = tmpl.varsTrans.(var);
    data = tblOut.(var);

    lb = pVar.logBase;

    % --- GLOBAL TRANS ---
    % Logit using Smithson & Verkuilen (2006) Squeeze transformation
    if ischar(lb) && strcmpi(lb, 'logit')
        N = length(data);
        data = (data * (N - 1) + 0.5) / N;
        data = log(data ./ (1 - data));
    end

    % Offset
    if pVar.offset > 0
        data = data + pVar.offset;
    end

    % Log
    if (isnumeric(lb) && ~isempty(lb)) || (ischar(lb) && strcmpi(lb, 'e'))
        if ischar(lb) % 'e'
            data = log(data);
        else
            data = log(data) ./ log(lb);
        end
    end

    % --- GROUP TRANS ---
    stats = pVar.stats;

    vecMean = nan(height(tbl), 1);
    vecSD   = nan(height(tbl), 1);
    vecRef  = nan(height(tbl), 1);

    validRows = (idxGrp > 0);
    validIdx  = idxGrp(validRows);

    vecMean(validRows) = stats.Mean(validIdx);
    vecSD(validRows)   = stats.SD(validIdx);

    % Apply Norm
    if pVar.flgNorm
        vecRef(validRows) = stats.RefMean(validIdx);
        isOK = ~isnan(vecRef) & vecRef ~= 0;

        data(~isOK) = NaN;
        data(isOK)  = 100 + ((data(isOK) - vecRef(isOK)) ./ abs(vecRef(isOK)) * 100);
    end

    % Apply Z
    if pVar.flgZ
        isOK = ~isnan(vecMean) & ~isnan(vecSD);
        data(~isOK) = NaN;

        isConst = (vecSD == 0);
        data(isOK & ~isConst) = (data(isOK & ~isConst) - vecMean(isOK & ~isConst)) ./ vecSD(isOK & ~isConst);
        data(isOK & isConst)  = data(isOK & isConst) - vecMean(isOK & isConst);
    end

    tblOut.(var) = data;
end
end


%% ========================================================================
%  CORE: INVERSE TRANSFORM
%  ========================================================================
function [tblOut, params] = trans_inv(tbl, tmpl, args, varsTrans)

verbose = args.verbose;

params = tmpl;
tblOut = tbl;

% Match Groups
idxGrp = match_grp(tbl, tmpl);

for iVar = 1:numel(varsTrans)
    var = varsTrans{iVar};
    pVar = tmpl.varsTrans.(var);
    data = tblOut.(var);

    stats = pVar.stats;
    vecMean = stats.Mean(idxGrp);
    vecSD   = stats.SD(idxGrp);

    % --- Step 1: Reverse Group Transforms (Z or Norm) ---

    % Reverse Z-Score
    if pVar.flgZ
        isConst = (vecSD == 0);
        data(~isConst) = data(~isConst) .* vecSD(~isConst) + vecMean(~isConst);
        data(isConst)  = data(isConst) + vecMean(isConst);
    end

    % Reverse Norm
    if pVar.flgNorm
        vecRef = stats.RefMean(idxGrp);
        isOK = ~isnan(vecRef) & vecRef ~= 0;
        isZero = (vecRef == 0);

        data(isOK)   = ((data(isOK) - 100) ./ 100 .* abs(vecRef(isOK))) + vecRef(isOK);
        data(isZero) = data(isZero) - 100;
        data(~isOK & ~isZero) = NaN;
    end

    % --- Step 2: Reverse Global Transforms ---

    lb = pVar.logBase;

    % Reverse Log
    if (isnumeric(lb) && ~isempty(lb)) || (ischar(lb) && strcmpi(lb, 'e'))
        if ischar(lb)
            data = exp(data);
        else
            data = lb .^ data;
        end
    end

    % Reverse Offset (Must happen after exp if forward was log(x+c))
    if pVar.offset > 0
        data = data - pVar.offset;
    end

    % Reverse Logit
    if ischar(lb) && strcmpi(lb, 'logit')
        % 1. Sigmoid: y = log(p / (1-p)) -> p = 1 / (1 + exp(-y))
        data = 1 ./ (1 + exp(-data));

        % 2. Reverse "Squeeze" (Smithson & Verkuilen)
        % Forward: y' = (y(N-1) + 0.5) / N
        % Inverse: y  = (y' * N - 0.5) / (N - 1)
        N = length(data);
        if N > 1
            data = (data * N - 0.5) / (N - 1);
        end
    end

    tblOut.(var) = data;
    if verbose, fprintf('[%s] Applied Inverse Transform.\n', var); end
end

end


%% ========================================================================
%  HELPER: MATCH GROUPS
%  ========================================================================
function idxGrp = match_grp(tbl, tmpl)
% Calculates group indices for Apply/Inverse modes by matching against template.

varsGrp = tmpl.varsGrp;
tblVars = tbl.Properties.VariableNames;

if ~isempty(varsGrp)
    % Ensure grouping vars exist in table
    if all(ismember(varsGrp, tblVars))
        fns = fieldnames(tmpl.varsTrans);
        if isempty(fns)
            % Should rarely happen if we have groups
            warning('tbl_trans:NoVars', 'Template contains no variable parameters to retrieve groups from.');
            idxGrp = ones(height(tbl), 1);
        else
            % Retrieve unique groups from the statistics of the first variable
            pFirst = tmpl.varsTrans.(fns{1});
            uGrps  = pFirst.stats(:, varsGrp);

            % Match rows
            [~, idxGrp] = ismember(tbl(:, varsGrp), uGrps, 'rows');

            % Handle unknown groups
            if any(idxGrp == 0)
                error('tbl_trans:NewGroups', ...
                    '%d rows have groups not seen in Template. Using default (NaN or 1).', ...
                    sum(idxGrp==0));
            end
        end
    else
        warning('tbl_trans:MissingGroups', 'Grouping variables missing in table. Using default group.');
        idxGrp = ones(height(tbl), 1);
    end
else
    idxGrp = ones(height(tbl), 1);
end

end
