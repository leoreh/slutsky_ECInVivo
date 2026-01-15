function hFig = lme_plotRes(mdl, varargin)
% LME_PLOTRES Diagnostic visualizations for LME/GLME models.
%
%   LME_PLOTRES(MDL) generates a set of diagnostic plots to evaluate model
%   assumptions, specifically targeting heteroscedasticity and normality of
%   residuals and random effects.
%
%   The function produces the following plots:
%       - Scale-Location: Square root of absolute standardized residuals vs.
%         fitted values. Used to check for homoscedasticity.
%       - BLUP Q-Q: Normal probability plots for the Best Linear Unbiased
%         Predictors (Random Effects). Checked per grouping variable/term.
%       - Residual Q-Q (LME Only): Normal probability plot of Standardized
%         Conditional Residuals to check for normality of error terms.
%
%   INPUTS:
%       mdl         - (object) Fitted LinearMixedModel or GeneralizedLinearMixedModel.
%       varargin    - (param/value) Optional parameters:
%                     'Name'    : (char) Figure name (default: 'LME Diagnostics').
%
%   OUTPUTS:
%       hFig        - (handle) Handle to the created figure.
%
%   See also: PLOTRESIDUALS, RANDOMEFFECTS, LME_ANALYSE

%% ========================================================================
%  SETUP
%  ========================================================================

p = inputParser;
addRequired(p, 'mdl', @(x) isa(x, 'LinearMixedModel') || isa(x, 'GeneralizedLinearMixedModel'));
addParameter(p, 'Name', 'LME Diagnostics', @ischar);
parse(p, mdl, varargin{:});

figName = p.Results.Name;

% Determine Model Type
isLME = isa(mdl, 'LinearMixedModel');

hFig = figure('Color', 'w', 'Name', figName);
t = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');


%% ========================================================================
%  CALCULATE RESIDUALS (STANDARDIZED)
%  ========================================================================
%  We use Standardized residuals for both Scale-Location and Q-Q plots.
%  Note: GLME objects in older MATLAB versions might not support 'Standardized'
%  explicitly in plotResiduals, but the residuals method supports it.

try
    resStd = residuals(mdl, 'ResidualType', 'Standardized');
catch
    warning('Standardized residuals not supported for this model type. Falling back to Pearson.');
    resStd = residuals(mdl, 'ResidualType', 'Pearson');
end
yFit = fitted(mdl);


%% ========================================================================
%  SCALE-LOCATION
%  ========================================================================
%  Y: Sqrt(|Standardized Residuals|)
%  X: Fitted Values
%  Trend: Smoothed line (Loess)

nexttile;
hold on;

% Transform
yScale = sqrt(abs(resStd));

% Scatter
scatter(yFit, yScale, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.4);

% Trend Line (Smoothing spline/Lowess)
% Sort for plotting line
[xSort, idxSort] = sort(yFit);
ySort = yScale(idxSort);

% Simple moving average or robust smooth if available
% Using 'rloess' reduces outlier influence
ySmooth = smooth(xSort, ySort, 0.5, 'rloess');
plot(xSort, ySmooth, 'r-', 'LineWidth', 2);

title('Scale-Location');
xlabel('Fitted Values');
ylabel('\surd|Standardized Residuals|');
grid on;
axis tight;


%% ========================================================================
%  BLUP Q-Q PLOTS (Random Effects)
%  ========================================================================
%  Check normality of random effects (b)

[b, bStats] = randomEffects(mdl);
% bStats contains 'Name' (Grouping Var) and 'Group' (Term Name)

% Identify unique terms (e.g., (Intercept) in Subject, Days in Subject)
% Concatenate Name and Group to find unique combinations
[~, ~, idxTerm] = unique(strcat(bStats.Name, bStats.Group));
termIDs = unique(idxTerm);

for iRand = 1:length(termIDs)
    mask = (idxTerm == termIDs(iRand));
    bSubset = b(mask);

    % Term Label
    currentName = bStats.Name(mask);
    currentGroup = bStats.Group(mask);

    groupName = char(currentName(1));
    termName = char(currentGroup(1));

    plotTitle = sprintf('BLUP Q-Q: %s | %s', termName, groupName);

    nexttile;
    % Use qqplot (Statistical Toolbox)
    hQ = qqplot(bSubset);

    title(plotTitle);
    xlabel('Theoretical Quantiles');
    ylabel('Sample Quantiles');
    grid on;

    % Tweak appearance to match style
    if ~isempty(hQ)
        set(hQ(1), 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'Marker', 'o', 'MarkerSize', 4);
        set(hQ(2), 'Color', 'r', 'LineWidth', 1.5); % The reference line
        set(hQ(3), 'Color', 'r', 'LineStyle', '--'); % Bounds if present
    end
end


%% ========================================================================
%  RESIDUAL Q-Q PLOT (LME Only)
%  ========================================================================
%  Standardized residuals vs Normal distribution.
%  Concept checks for Normality of Epsilon.

if isLME
    nexttile;

    % Use qqplot on the pre-calculated standardized residuals
    hQ = qqplot(resStd);

    % Customize title/labels
    title('Residual Q-Q (Standardized)');
    xlabel('Theoretical Quantiles');
    ylabel('Sample Quantiles');
    grid on;

    % Tweak appearance
    if ~isempty(hQ)
        set(hQ(1), 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none', 'Marker', 'o', 'MarkerSize', 4);
        set(hQ(2), 'Color', 'r', 'LineWidth', 1.5);
        set(hQ(3), 'Color', 'r', 'LineStyle', '--');
    end
end


% Improve Layout
set(findall(hFig, '-property', 'FontSize'), 'FontSize', 10);
title(t, figName, 'FontSize', 12, 'FontWeight', 'bold');

end     % EOF
