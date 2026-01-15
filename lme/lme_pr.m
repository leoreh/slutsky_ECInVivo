function [hFig, tblRes] = lme_pr(mdl, vars, varargin)
% LME_PR Plots Partial Residuals (Added Variable Plot) for LME models.
%
%   [HFIG, TBLRES] = LME_PR(MDL, VARS, ...)
%   Visualizes the partial residuals of a model to check for non-linearity
%   or validate interaction terms.
%
%   METHOD:
%       1. Create a "Reduced Model" by removing the specified variable(s)
%          and their interactions from the original formula.
%       2. Calculate the Residuals of this reduced model.
%       3. Plot these Residuals against the variable of interest.
%
%   INPUTS:
%       mdl         - (Required) Fitted LME/GLME object.
%       vars        - (Required) Cell array of 1 or 2 variable names to
%                     visualize/remove.
%                     - 1 Var: Plot residuals vs Var.
%                     - 2 Vars: Plot residuals vs Var1, grouped by Var2.
%       ...         - (Optional) Name-Value Pairs:
%                     'hAx'          - Axis handle to plot into.
%                     'clr'          - Color matrix or 'auto'.
%                     'transParams'  - (struct) Transformation parameters from
%                                      LME_ANALYSE / TBL_TRANS. If provided,
%                                      predictors are back-transformed.
%                     'verbose' - (logical) Print progress.
%
%   OUTPUTS:
%       hFig        - Figure handle.
%       tblRes      - Table containing the data and calculated residuals.
%
%   See also: LME_PD, LME_FIT, LME_FRML2RMV

%% ========================================================================
%  ARGUMENT PARSING
%  ========================================================================

p = inputParser;
addRequired(p, 'mdl');
addRequired(p, 'vars', @(x) ischar(x) || iscell(x));
addParameter(p, 'hAx', [], @(x) isempty(x) || isgraphics(x));
addParameter(p, 'clr', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'transParams', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'verbose', true, @islogical);

parse(p, mdl, vars, varargin{:});

hAx         = p.Results.hAx;
clr         = p.Results.clr;
transParams = p.Results.transParams;
verbose     = p.Results.verbose;

% Ensure vars is cell
if ischar(vars), vars = {vars}; end
if length(vars) > 2
    error('[LME_PR] Can only visualize up to 2 variables.');
end

%% ========================================================================
%  FIT REDUCED MODEL
%  ========================================================================

if verbose, fprintf('[LME_PR] Fitting Reduced Model...\n'); end

% Original Formula & Data
origFrml = char(mdl.Formula);
tblRes   = mdl.Variables;

% Determine Reduced Formula
% We iterate through 'vars' and remove each one from the formula.
redFrml = origFrml;
for iVar = 1:length(vars)
    redFrml = lme_frml2rmv(redFrml, vars{iVar});
end

if verbose
    fprintf('   Original: %s\n', origFrml);
    fprintf('   Reduced:  %s\n', redFrml);
end

% Determine Distribution
dist = mdl.Distribution;
if isobject(dist) % Handle newer MATLAB versions/GLME objects
    dist = char(dist.Name);
end

% Fit Reduced Model
mdlRed = lme_fit(tblRes, redFrml, 'dist', dist);

% Calculate Residuals
% Residual = Observed - Predicted(Reduced)
tblRes.Resid = residuals(mdlRed, 'ResidualType', 'Raw');


%% ========================================================================
%  BACK TRANSFORM
%  ========================================================================

tblFit = tblRes; % Keep transformed data for fitting the "Trend Line"

if ~isempty(transParams)
    if verbose, fprintf('[LME_PR] Back-transforming predictors\n'); end
    tblPlot = tbl_trans(tblRes, 'template', transParams, 'flgInv', true);
else
    tblPlot = tblRes;
end


%% ========================================================================
%  PLOT
%  ========================================================================

if isempty(hAx)
    hFig = figure('Color', 'w', 'Name', 'Partial Residual Plot');
    hAx = axes('Parent', hFig);
else
    hFig = ancestor(hAx, 'figure');
end
hold(hAx, 'on');

% Identify X and Grouping
varsPlot = vars;
isNum    = cellfun(@(x) isnumeric(tblPlot.(x)) && ~iscategorical(tblPlot.(x)), varsPlot);

if isscalar(varsPlot)
    % 1D Case
    xName = varsPlot{1};
    grpName = [];
    isXNum = isNum(1);

    gIds = ones(height(tblPlot), 1);
    gLbls = {'All'};
    nGrps = 1;
else
    % 2D Case
    var1 = varsPlot{1};
    var2 = varsPlot{2};

    % Heuristic: If one is Categorical, it becomes the Grouping variable.
    % If both Continuous, the second one is Grouping.
    if ~isNum(1) && isNum(2)
        xName = var2; grpName = var1; isXNum = true;
    else
        xName = var1; grpName = var2; isXNum = isNum(1);
    end

    [gIds, gLbls] = findgroups(tblPlot.(grpName));
    nGrps = length(gLbls);
end

% Colors
if isempty(clr)
    cfg = mcu_cfg();
    if nGrps == 2
        clr = cfg.clr.grp;
    else
        clr = lines(nGrps);
    end
end

% Plot Loop
for iG = 1:nGrps
    idx = (gIds == iG);
    subFit  = tblFit(idx, :);  % For Fitting (Transformed)
    subPlot = tblPlot(idx, :); % For Plotting (Original)

    % Extract Data
    xFit  = subFit.(xName);
    xPlot = subPlot.(xName);
    yData = subPlot.Resid; % Residuals are same in both (created after fit)

    c = clr(iG, :);
    dispName = string(gLbls(iG));

    % Pre-process for Scatter
    if ~isnumeric(xPlot) && ~islogical(xPlot) && ~iscategorical(xPlot)
        xPlot = categorical(xPlot);
    end

    % Scatter Points (Back-Transformed X)
    scatter(hAx, xPlot, yData, 20, c, 'filled', ...
        'MarkerFaceAlpha', 0.4, 'HandleVisibility', 'off');

    % Fit Line (in Transformed Space)
    if isXNum
        % Fit Linear Trend on TRANSFORMED data (y ~ x_trans + c)
        lm = fitlm(xFit, yData);

        % Generate Grid in Transformed Space
        xGridTrans = linspace(min(xFit), max(xFit), 100)';
        yPred = predict(lm, xGridTrans);

        % Back-Transform X-Grid for Plotting
        if ~isempty(transParams) && isfield(transParams.varsTrans, xName)
            % Create dummy table for inversion
            tmp = table(xGridTrans, 'VariableNames', {xName});
            tmp = tbl_trans(tmp, 'template', transParams, 'varsInc', {xName}, 'flgInv', true);
            xGridOrign = tmp.(xName);
        else
            xGridOrign = xGridTrans;
        end

        plot(hAx, xGridOrign, yPred, '-', 'Color', c, 'LineWidth', 2, ...
            'DisplayName', dispName);
    else
        % Categorical X: Plot Mean of Residuals
        xCats = double(categorical(xPlot));
        uCats = unique(xCats);

        % Calculate Means
        yMeans = splitapply(@mean, yData, findgroups(xPlot));

        % Plot connecting line
        plot(hAx, uCats, yMeans, 'o-', 'Color', c, 'LineWidth', 2, ...
            'MarkerFaceColor', 'w', 'DisplayName', dispName);

        xticklabels(hAx, categories(categorical(xPlot)));
    end
end

% Labels & Aesthetics
xlabel(hAx, xName, 'Interpreter', 'none');
ylabel(hAx, sprintf('Residuals (%s | Reduced)', mdl.ResponseName), 'Interpreter', 'none');
title(hAx, 'Partial Residual Plot', 'FontWeight', 'normal');
grid(hAx, 'on');

if ~isempty(grpName)
    lgd = legend(hAx, 'Location', 'best', 'Interpreter', 'none');
    title(lgd, grpName);
end

if ~isempty(transParams) && isfield(transParams.varsTrans, xName) ...
        && ~isempty(transParams.varsTrans.(xName).logBase)
    set(hAx, 'XScale', 'log')
end

end         % EOF
