function [tblRes, hFig] = lme_pr(mdl, varX, varargin)
% LME_PR Plots Partial Residuals (CPR) or Partial Regression (AVP) for LME models.
%
%   [HFIG, TBLRES] = LME_PR(MDL, VARX, ...)
%
%   This function visualizes the marginal relationship between a target predictor
%   (VARX) and the response, after accounting for all other variables in the model.
%
%   MODES:
%       'residual'   (Default) Component-Plus-Residual Plot (CPR).
%                    Plots (Residuals + Component) vs Raw Predictor.
%                    * Interpretation: Functional Form (Linearity).
%                    * Units: X-axis is Raw.
%
%       'regression' Added Variable Plot (AVP).
%                    Plots (Y Residuals) vs (X Residuals).
%                    * Interpretation: Partial Correlation & Effect Size.
%                    * Units: X-axis is Standardized (if lme_analyse was used).
%                      The slope = The Model Coefficient (Beta).
%
%   INPUTS:
%       mdl         - (Required) Fitted LME object.
%       varX        - (Required) Name of the target predictor (char).
%       varargin    - (Optional) Key-Value pairs:
%                     'varGrp'      : (char) Variable to group by (coloring).
%                     'flgMode'     : 'residual' | 'regression'.
%                     'distX'       : Distribution for X model (AVP only). Default: 'Normal'.
%                     'transParams' : (struct) Output from LME_ANALYSE/TBL_TRANS.
%                                     If provided in 'regression' mode, X-axis
%                                     is un-standardized to Model Units.
%                     'hAx', 'clr'  : Graphics handles and colors.
%
%   EXAMPLE (Spike Allocation):
%       % Does Burst Recovery (dBrst) predict Single Recovery (dSngl),
%       % controlling for baseline firing (fr) and brightness (pBspk)?
%
%       frml = 'dSngl_rel ~ (dBrst_rel + pBspk + fr) * Group + (1|Name)';
%       mdl = lme_analyse(tbl, frml, 'dist', 'normal');
%
%       % Visualize the Partial Regression (AVP)
%       lme_pr(mdl, 'dBrst_rel', 'flgMode', 'regression', 'varGrp', 'Group');
%
%       % The slopes of the pr plot will be identicle to the estimate of the
%       % model for dBrst_rel (the fit must be linear and not orthogonal.
%
%   OUTPUTS:
%       tblRes      - Table with data and residuals.
%       hFig        - Figure handle.
%
%   See also: LME_ANALYSE, LME_FIT, PLOT_LINEREG

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'mdl');
addRequired(p, 'varX', @ischar);

% Optional
addParameter(p, 'varGrp', '', @ischar);
addParameter(p, 'flgMode', 'residual', @(x) any(validatestring(x, {'residual', 'regression'})));
addParameter(p, 'transParams', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'distX', 'Normal', @ischar);
addParameter(p, 'hAx', [], @(x) isempty(x) || isgraphics(x));
addParameter(p, 'clr', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'verbose', true, @islogical);

parse(p, mdl, varX, varargin{:});

varGrp      = p.Results.varGrp;
flgMode     = p.Results.flgMode;
distX       = p.Results.distX;
transParams = p.Results.transParams;
hAx         = p.Results.hAx;
clr         = p.Results.clr;
verbose     = p.Results.verbose;


%% ========================================================================
%  PREPARATION
%  ========================================================================

% Data & Groups
tblMdl = mdl.Variables;
if ~isempty(varGrp) && ismember(varGrp, tblMdl.Properties.VariableNames)
    [grpIds, grpNames] = findgroups(tblMdl.(varGrp));
    nGrps = length(grpNames);
else
    grpIds = ones(height(tblMdl), 1);
    grpNames = {'All'};
    nGrps = 1;
end

% Graphic Handles
if isempty(hAx)
    hFig = figure('Color', 'w', 'Name', ['Partial ' flgMode]);
    hAx = axes('Parent', hFig);
else
    hFig = ancestor(hAx, 'figure');
end
hold(hAx, 'on');

% Colors
if isempty(clr)
    cfg = mcu_cfg();
    if nGrps == 2
        clr = cfg.clr.grp;
    else
        clr = lines(nGrps);
    end
end


%% ========================================================================
%  CALCULATE RESIDUALS (Y)
%  ========================================================================
%  Regress Y against everything EXCEPT the target predictor X.

if verbose
    fprintf('[LME_PR] Mode: %s. Fitting Reduced Y Model...\n', flgMode);
end

% Construct Reduced Formula (Y ~ Rest)
frml = char(mdl.Formula);

% Remove the target predictor (X) from the formula
frmlY = lme_frml2rmv(frml, varX);

% Fit Reduced Y Model
distY = 'Normal';
if isa(mdl, 'GeneralizedLinearMixedModel')
    distY = mdl.Distribution;
end

mdlRedY = lme_fit(tblMdl, frmlY, 'dist', distY);
tblMdl.ResidY = residuals(mdlRedY, 'ResidualType', 'Raw');


%% ========================================================================
%  CALCULATE RESIDUALS (X) - [REGRESSION MODE ONLY]
%  ========================================================================
%  Regress X against everything else to isolate its unique variance.

if strcmpi(flgMode, 'regression')
    if verbose
        fprintf('[LME_PR] Fitting Reduced X Model (%s)...\n', distX);
    end

    % Construct Reduced Formula (X ~ Rest)
    rhs = extractAfter(frmlY, '~');
    frmlX = [varX ' ~ ' rhs];

    % Fit Reduced X Model
    mdlRedX = lme_fit(tblMdl, frmlX, 'dist', distX);
    tblMdl.ResidX = residuals(mdlRedX, 'ResidualType', 'Raw');

    % Un-standardize X to Model Units
    xUnits = 'Standardized (Z-Score)';
    if ~isempty(transParams)

        pVar = transParams.varsTrans.(varX);
        if pVar.flgZ
            % Multiply by SD (Intercept is 0 for residuals)
            % This restores the mechanistic slope (Beta_raw = Beta_std / SD)
            sigma = pVar.stats.SD(1); % Assume global Z-scoring (nGrps=1)
            if sigma ~= 0
                tblMdl.ResidX = tblMdl.ResidX .* sigma;
                xUnits = 'Model Units (Centered)';
                if verbose
                    fprintf('[LME_PR] Un-standardizing X-axis by SD=%.4f\n', sigma);
                end
            end
        end
    end

    varX_Plot = 'ResidX';

else
    % 'residual' mode: Plot against Raw X
    varX_Plot = varX;
end


%% ========================================================================
%  PLOT
%  ========================================================================

lims = max(abs([tblMdl.(varX_Plot); tblMdl.ResidY])) * 1.1;

% Determine Regression Type
regType = 'ortho';
if strcmpi(flgMode, 'regression')
    if exist('xUnits', 'var') && strcmpi(xUnits, 'Standardized (Z-Score)')
        regType = 'linear';
    end
end

% Grouping Variable
gArg = [];
if ~isempty(varGrp) && ismember(varGrp, tblMdl.Properties.VariableNames)
    gArg = tblMdl.(varGrp);
end

% Plot via plot_scat
plot_scat(tblMdl, varX_Plot, 'ResidY', ...
    'g', gArg, ...
    'c', clr, ...
    'fitType', regType, ...
    'alpha', 0.4, ...
    'sz', 20, ...
    'hAx', hAx, ...
    'flgStats', true); % Stats will trigger fit text in plot_scat -> plot_lineReg
hold(hAx, 'on');

% Reference Lines
xlim(hAx, [-lims, lims]); ylim(hAx, [-lims, lims]);
plot(hAx, [-lims, lims], [-lims, lims], '--k', 'HandleVisibility', 'off'); % Identity
xline(hAx, 0, '-', 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
yline(hAx, 0, '-', 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');

% Decoration
grid(hAx, 'on');
box(hAx, 'on');

if strcmpi(flgMode, 'regression')
    if exist('xUnits', 'var')
        lblX = sprintf('%s | %s', varX, xUnits);
    else
        lblX = [varX ' | Reduced (Residuals)'];
    end
    xlabel(hAx, lblX, 'Interpreter', 'none');
    ylabel(hAx, [mdl.ResponseName ' | Reduced (Residuals)'], 'Interpreter', 'none');
    title(hAx, 'Partial Regression (AVP)');
else
    xlabel(hAx, varX, 'Interpreter', 'none');
    ylabel(hAx, [mdl.ResponseName ' | Reduced (Residuals)'], 'Interpreter', 'none');
    title(hAx, 'Partial Residuals (CPR)');
end

if nGrps > 1
    legend(hAx, 'Location', 'best', 'Interpreter', 'none');
end

tblRes = tblMdl;

end     % EOF
