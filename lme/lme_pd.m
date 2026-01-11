function [pdRes, hFig] = lme_pd(mdl, vars, varargin)
% LME_PD Plots Partial Dependence (Marginal Effects) for LME/GLME models.
%
%   [PDRES, HFIG] = LME_PD(MDL, VARS, ...)
%   Calculates and visualizes the marginal effect of one or two predictors
%   on the response, while holding all other predictors constant (at Mean
%   or Mode).
%
%   INPUTS:
%       mdl         - (Required) Fitted LME/GLME object.
%       vars        - (Required) Cell array of 1 or 2 variable names to vary.
%                     e.g. {'Group'} or {'Group', 'pBspk'}.
%       ...         - (Optional) Name-Value Pairs:
%                     'hAx'          - Axis handle to plot into.
%                     'nGrid'        - Number of grid points for continuous vars (default 100).
%                     'confLvl'      - Confidence Level (default 0.95).
%                     'clr'          - Color matrix or 'auto'.
%                     'transParams'  - (struct) Transformation parameters from
%                                      LME_ANALYSE / TBL_TRANSFORM. If provided,
%                                      predictions and grid values are back-
%                                      transformed to original units.
%
%   OUTPUTS:
%       pdRes       - Table containing the grid, predictions, and CIs.
%                     (In original units if transParams is provided).
%       hFig        - Figure handle.
%
% =========================================================================

%% ========================================================================
%  ARGUMENT PARSING
%  ========================================================================

p = inputParser;
addRequired(p, 'mdl');
addRequired(p, 'vars', @(x) ischar(x) || iscell(x));
addParameter(p, 'hAx', [], @(x) isempty(x) || isgraphics(x));
addParameter(p, 'nGrid', 100, @isnumeric);
addParameter(p, 'confLvl', 0.95, @isnumeric);
addParameter(p, 'clr', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'transParams', [], @(x) isempty(x) || isstruct(x));

parse(p, mdl, vars, varargin{:});

hAx     = p.Results.hAx;
nGrid   = p.Results.nGrid;
confLvl = p.Results.confLvl;
clr     = p.Results.clr;
transP  = p.Results.transParams;

% Ensure vars is cell
if ischar(vars), vars = {vars}; end
if length(vars) > 2
    error('LME_PD:MaxVars', 'Can only visualize up to 2 variables.');
end


%% ========================================================================
%  PREDICT
%  ========================================================================

% 1. Generate Grid (Transformed Space)
%    We use the model's data, so this grid is in the same space as training.
[tblGrid, varyInfo] = get_grid(mdl, vars, nGrid);

% 2. Predict (Transformed Space)
%    Use 'Conditional', false to compute marginal (population-level) effects.
%    This sets random effects to 0 (population mean).
alpha = 1 - confLvl;
[yPred, yCI] = predict(mdl, tblGrid, 'Conditional', false, ...
    'Prediction', 'curve', 'Alpha', alpha);

% 3. Back-Transform (to Original Space)
%    If transParams are provided, we invert the transformations.
%    We use 'tbl_transform' in 'inverse' mode.

plotGrid = tblGrid; % Copy for plotting
yPredPlot = yPred;
yCIPlot   = yCI;

if ~isempty(transP)

    % A. PREDICTORS
    % table `plotGrid` already has all predictors and grouping vars.
    plotGrid = tbl_transform(plotGrid, 'template', transP, 'flgInv', true);

    % B. RESPONSE
    % predictions are a vector. to use tbl_transform, we must wrap them
    % in a table. Crucially, we must include the GROUPING variables
    % from `plotGrid` so that `tbl_transform` can look up the correct
    % group-specific stats (if Z-scoring/Norm was grouped).

    respName = mdl.ResponseName;
    if isfield(transP.varList, respName)

        % Identify grouping vars needed
        varsGrp = transP.varsGrp;

        % Create Temp Table: [Groups | Response | CI_L | CI_U]
        % We actually handle Main Prediction and Bounds separately or together?
        % Simpler to do one by one or create 3 columns if they are same var?
        % tbl_transform expects specific variable names.

        % Strategy: Create a table with 'respName' and reuse it 3 times
        % or just do it once if we trust the loop.

        % Construct base table with Groups
        tblRspBase = table();
        if ~isempty(varsGrp)
            % Copy grouping columns from grid
            tblRspBase = plotGrid(:, varsGrp);
        else
            % Dummy col to ensure height matches
            tblRspBase = table((1:height(plotGrid))');
        end

        % Helper to invert a vector treated as 'respName'
        batch_inv = @(vec) local_inv_vec(vec, respName, tblRspBase, transP);

        yPredPlot    = batch_inv(yPredPlot);
        yCIPlot(:,1) = batch_inv(yCIPlot(:,1));
        yCIPlot(:,2) = batch_inv(yCIPlot(:,2));
    end
end


% Store results
pdRes = plotGrid;
pdRes.Pred  = yPredPlot;
pdRes.Lower = yCIPlot(:, 1);
pdRes.Upper = yCIPlot(:, 2);


%% ========================================================================
%  PLOT
%  ========================================================================

if isempty(hAx)
    hFig = figure('Color', 'w', 'Name', 'Partial Dependence');
    hAx = axes('Parent', hFig);
else
    hFig = ancestor(hAx, 'figure');
end
hold(hAx, 'on');

% Identify X and Grouping
varsPlot = varyInfo.vars;
isNum    = varyInfo.isNum;

if isscalar(varsPlot)
    % 1D Case
    xName = varsPlot{1};
    grpName = [];
    isXNum = isNum(1);

    % Dummy Grouping
    gIds = ones(height(pdRes), 1);
    gLbls = {'All'};
    nGrps = 1;
else
    % 2D Case
    var1 = varsPlot{1};
    var2 = varsPlot{2};

    % Heuristic: If one is Categorical, it becomes the GROUPING variable.
    % If both Continuous, the second one is Grouping.
    if ~isNum(1) && isNum(2)
        xName = var2; grpName = var1; isXNum = true;
    else
        xName = var1; grpName = var2; isXNum = isNum(1);
    end

    [gIds, gLbls] = findgroups(pdRes.(grpName));
    nGrps = length(gLbls);
end

% Colors
if isempty(clr)
    try
        cfg = mcu_cfg();
        if nGrps == 2, clr = cfg.clr.grp; else, clr = lines(nGrps); end
    catch
        clr = lines(nGrps);
    end
end
if size(clr,1) < nGrps, clr = lines(nGrps); end


% Plot Loop
for iG = 1:nGrps
    idx = (gIds == iG);
    subTbl = pdRes(idx, :);

    % Extract Data
    xData = subTbl.(xName);
    y     = subTbl.Pred;
    yL    = subTbl.Lower;
    yU    = subTbl.Upper;

    c = clr(iG, :);
    dispName = string(gLbls(iG));

    if isXNum
        % --- CONTINUOUS X: Line + Shaded CI ---
        % Sort for line plotting
        [xData, idxSort] = sort(xData);
        y = y(idxSort); yL = yL(idxSort); yU = yU(idxSort);

        % Shaded CI (Manual Fill)
        xFill = [xData; flipud(xData)];
        yFill = [yL; flipud(yU)];

        % Use alpha blending
        fill(hAx, xFill, yFill, c, 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');

        % Mean Line
        plot(hAx, xData, y, '-', 'Color', c, 'LineWidth', 2, ...
            'DisplayName', dispName);

    else
        % --- CATEGORICAL X: Errorbar ---
        xCats = double(categorical(xData));
        uCats = unique(xCats);

        % Jitter/Offset if multiple groups
        offset = 0;
        if nGrps > 1
            offset = (iG - mean(1:nGrps)) * 0.15;
        end
        xPlot = xCats + offset;

        errorbar(hAx, xPlot, y, y-yL, yU-y, 'o-', ...
            'Color', c, 'LineWidth', 1.5, ...
            'MarkerFaceColor', 'w', 'CapSize', 8, ...
            'DisplayName', dispName);

        xticks(hAx, uCats);
        xticklabels(hAx, categories(categorical(xData)));
    end
end

% Labels
xlabel(hAx, xName, 'Interpreter', 'none');
ylabel(hAx, mdl.ResponseName, 'Interpreter', 'none');
grid(hAx, 'on');

if ~isempty(grpName)
    lgd = legend(hAx, 'Location', 'best', 'Interpreter', 'none');
    title(lgd, grpName);
end
title(hAx, 'Partial Dependence', 'FontWeight', 'normal');

end


%% ========================================================================
%  HELPER: GENERATE GRID
%  ========================================================================
function [tblGrid, varyInfo] = get_grid(mdl, varsVary, nGrid)
% GET_GRID Generates a synthetic "prediction grid" for marginal effects.
%
%   [GRIDTBL, VARYINFO] = GET_GRID(MDL, VARSVARY, NGRID)
%
%   'gridTbl' is a synthetic dataset used to probe the fitted model. It
%   isolates the effect of specific predictors (varsVary) by holding all
%   other predictors constant.
%
%   HOW IT WORKS:
%   1. Identify all predictors in the model.
%   2. For the variables of interest (varsVary):
%      - Create a range of values (e.g., 100 points from min to max).
%   3. For all OTHER variables (Fixed Effects):
%      - Fix them to a single representative value.
%      - Continuous Vars -> Average (Mean).
%      - Categorical Vars -> Most Common (Mode).
%   4. Generate 'gridTbl' by creating every possible combination.

% Identify predictors
allVars = mdl.Variables.Properties.VariableNames;
respName = mdl.ResponseName;
predVars = setdiff(allVars, respName);

% Initialize
gridArgs = {};
varyInfo.vars = varsVary;
varyInfo.isNum = false(size(varsVary));
data = mdl.Variables;

for iVar = 1:length(predVars)
    name = predVars{iVar};
    raw = data.(name);
    isVarying = ismember(name, varsVary);

    if isVarying
        % Varying: Generate Grid
        if isnumeric(raw) && ~iscategorical(raw)
            vals = linspace(min(raw), max(raw), nGrid)';
            varyInfo.isNum(strcmp(varsVary, name)) = true;
        else
            vals = unique(raw);
        end
        gridArgs{end+1} = vals; %#ok<AGROW>
    else
        % Fixed: Set to Mean/Mode
        if isnumeric(raw) && ~iscategorical(raw)
            vals = mean(raw, 'omitnan');
        else
            if iscell(raw)
                t = tabulate(raw);
                [~, idx] = max([t{:,2}]);
                vals = {t{idx,1}}; % Wrap for ndgrid
            else
                vals = mode(raw);
            end
        end
        gridArgs{end+1} = vals; %#ok<AGROW>
    end
end

[outGrids{1:length(predVars)}] = ndgrid(gridArgs{:});
tblArgs = cellfun(@(x) x(:), outGrids, 'UniformOutput', false);
tblGrid = table(tblArgs{:}, 'VariableNames', predVars);

end


%% ========================================================================
%  HELPER: VECTOR INVERSE (WRAPPER)
%  ========================================================================
function vecOut = local_inv_vec(vecIn, varName, tblBase, transP)
% Wraps a vector into a table to use tbl_transform

T = tblBase;
T.(varName) = vecIn; % Add data column

% Transform
T_Out = tbl_transform(T, 'template', transP, 'flgInv', true);

% Extract
vecOut = T_Out.(varName);

end


%% ========================================================================
%  NOTES
%  ========================================================================
%  1. VARIABLES & DIMENSIONS
%     This function currently supports only 1 or 2 varying variables.
%     Why? Because marginal effect plots are designed to show relationships
%     on a 2D plane (X vs Y).
%       - 1 Var:  Curve (X vs Y)
%       - 2 Vars: Interaction (Separate Curves for each Level, or Surface)
%     Expanding to 3+ variables requires 4D visualization (e.g., small
%     multiples or sliders), which significantly complicates the code and
%     interpretation. If you need a 3-way interaction (e.g., Group * Drug * Time),
%     it is best to filter your dataset or loop over one variable yourself
%     and generate multiple plots.
%
%  2. NEGATIVE VALUES & TRANSFORMATIONS
%     If you see values like -0.7 to 3.2 for a probability 'pBspk', it is
%     likely because the variable was Log (or Logit) transformed and/or
%     Z-scored during the preprocessing logic in 'lme_analyse'. The model
%     estimates coefficients in this *transformed* space.
%     To see the plot in natural units (0 to 1), you must back-transform
%     the X-axis and Y-axis. This function now supports a 'transParams'
%     inputs (available from 'lmeInfo'). If provided, the grid and
%     predictions will be automatically mapped back to the original scale.
%
%     Does this distort the CI?
%     No. Transforming the CI endpoints (Lower/Upper) is a statistically
%     valid way to obtain the CI on the original scale, provided the
%     transformation is monotonic (like Log or Exp).
%
%  3. INTERPRETATION
%     These plots show the *Marginal Effect* (Partial Dependence). They
%     represent the predicted outcome for an "average" unit as we vary
%     only the target variable(s). All other predictors are held fixed at
%     their global mean or mode.
%  ========================================================================
