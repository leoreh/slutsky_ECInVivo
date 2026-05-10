

%% ========================================================================
%  DEMONSTRATION: MARGINALIZATION vs RESIDUALIZATION
%  ========================================================================
%  When a model has multiple predictors, a raw scatter of X vs Y shows the
%  MARGINAL (unadjusted) relationship — conflating the direct effect of X
%  with indirect effects mediated by correlated covariates. To isolate the
%  PARTIAL (adjusted) effect, there are two strategies:
%
%  -----------------------------------------------------------------------
%  (A) MARGINALIZATION  — "Ask the model"
%  -----------------------------------------------------------------------
%  Keep the fitted model intact. Generate predictions while VARYING the
%  target predictor and HOLDING nuisance covariates constant. The model
%  does the controlling internally through its estimated coefficients.
%  Output: predicted CURVES with confidence intervals.
%
%   A1. LSMEANS / MEM  (Conditional Expectation at the Mean)
%       Fix each covariate to its MEAN (continuous) or MODE (categorical).
%       Predict for a grid of target predictor values.
%       Question answered: "What does a TYPICAL neuron look like?"
%       Use case: Publication figures. Clean predicted curves with CIs.
%       Function: lme_lsmeans
%
%   A2. PDP / AME  (Average Marginal Effect)
%       For each target value, predict using EVERY neuron's actual covariate
%       values, then average the predictions. Preserves the observed
%       covariate distribution rather than collapsing it to a single point.
%       Question answered: "What does the POPULATION look like on average?"
%       Use case: When covariate distributions are skewed or non-linear
%                 interactions make the "mean subject" unrepresentative.
%       Function: plotPartialDependence (MATLAB built-in, no CIs)
%
%   Note: In linear models without interactions between the target predictor
%   and covariates, LSMEANS and PDP give identical point estimates. They
%   diverge in non-linear models or when interactions are present.
%
%  -----------------------------------------------------------------------
%  (B) RESIDUALIZATION  — "Clean the data"
%  -----------------------------------------------------------------------
%  Remove the effect of nuisance covariates from the data itself, then
%  examine the leftover (residual) structure. The controlling happens by
%  subtracting predicted effects from the raw observations.
%  Output: scatter of ADJUSTED data points (no CIs).
%
%   B1. Partial Residual / CPR  (Component-Plus-Residual)
%       Y-axis: Residuals from a reduced model that excludes the target
%               predictor (i.e., the part of Y unexplained by covariates).
%       X-axis: Raw target predictor (original scale preserved).
%       Question answered: "Is the assumed functional form correct?"
%       Use case: Diagnostics — detecting non-linearity, heteroscedasticity.
%                 Keeps X in its natural units for intuitive reading.
%       Function: lme_pr(..., 'flgMode', 'residual')
%
%   B2. Partial Regression / AVP  (Added Variable Plot)
%       Y-axis: Residuals(Y | covariates) — Y cleaned of covariates.
%       X-axis: Residuals(X | covariates) — X cleaned of covariates.
%       Both axes centered at 0. The slope equals the model coefficient (Beta).
%       Question answered: "What is the unique partial correlation, net of
%                           all shared variance with covariates?"
%       Use case: Effect validation — confirming coefficient magnitude,
%                 detecting influential outliers and leverage points,
%                 assessing collinearity (narrow x-spread = high collinearity).
%       Function: lme_pr(..., 'flgMode', 'regression')
%
%  -----------------------------------------------------------------------
%  AXIS TRANSFORMS: WHEN AND WHY
%  -----------------------------------------------------------------------
%  Variables like firing rates (log-normal) and proportions (bounded [0,1])
%  often need transforms. But there are TWO distinct reasons to transform,
%  and confusing them causes the "curving fit line" problem:
%
%  1. MODELING TRANSFORM (statistical, before fitting):
%     Changes the functional form the model assumes.
%     - Response: lme_analyse('dist','log-normal') applies log(Y) so the
%       model assumes a multiplicative (log-linear) relationship.
%     - Predictors: lme_analyse auto-applies log10 if skewness > 2.
%     These are stored in transParams and define "model space." The model's
%     coefficients, residuals, and predictions all live in this space.
%
%  2. DISPLAY TRANSFORM (cosmetic, after analysis):
%     Stretches the axis for visual clarity (e.g., set(gca,'YScale','log')).
%     Doesn't change any statistics or model behavior.
%
%  THE PITFALL:
%     If you plot raw Y and fit a linear line, then set the Y axis to log,
%     the fit CURVES — because the fit was computed in raw space but
%     displayed on a log axis. Fix: transform the data BEFORE fitting so
%     the fit and display agree, or don't show a fit line at all.
%
%  AUTO-LOG IN lme_analyse (common surprise):
%     Even with 'flgStnd', false (which disables z-scoring), lme_analyse
%     still calls tbl_trans('logBase', 10, 'skewThr', 2) on ALL numeric
%     predictors. If a predictor's skewness > 2, it gets log10-transformed.
%     This means pBurst (a right-skewed proportion) enters the model as
%     log10(pBurst), NOT raw pBurst. The Y-axis residuals and X-axis
%     predictor values in lme_pr all live in this auto-transformed space.
%
%  WHICH PANELS NEED WHAT:
%     Raw scatter  — Data is in natural units. Apply display transforms
%                    (logit for proportions, log for rates) AND transform
%                    the data BEFORE passing to plot_scat, so the fit line
%                    matches the axis. This is purely cosmetic.
%     LSMEANS      — lme_lsmeans back-transforms predictions to natural
%                    units (Hz). Use log Y axis for readability. No logit
%                    on X: the model sweeps a grid over pBurst, and the
%                    curve naturally captures the shape. No mismatch.
%     CPR          — Y-axis residuals are in model space (log units,
%                    centered at 0). The X-axis in lme_pr shows model-
%                    space X (auto-log10), but for display we back-
%                    transform to raw and apply logit (natural for
%                    proportions). The fit line is recomputed in this
%                    display space. This is a display choice — the model
%                    coefficient comes from the AVP panel, not CPR.
%     AVP          — Both axes are residualized (model space, centered at
%                    0, can be negative). No transforms possible. Note:
%                    lme_pr calls plot_lineEq which forces symmetric axes
%                    (x matches y range) — override xlim afterward.
%
%  -----------------------------------------------------------------------
%  EXAMPLE BELOW
%  -----------------------------------------------------------------------
%  Question: How does baseline burstiness (pBurst) predict steady-state
%  single-spike firing (ss_frSingle), controlling for steady-state burst
%  firing (ss_frBurst), across genotypes?
%  -----------------------------------------------------------------------

% Fit the model
frml = 'ss_frSingle ~ (pBurst + ss_frBurst) * genotype + (1|sbjID)';
dist = 'log-normal';
[mdlDemo, ~, infoDemo, ~] = lme_analyse(tbl, frml, ...
    'dist', dist, 'flgStnd', false, 'verbose', false);

% Figure
hFig = figure('Color', 'w', 'Name', 'Marginalization vs Residualization');
hTile = tiledlayout(hFig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
cfg = mcu_cfg();

% Display-transformed data for Panel 1 (see AXIS TRANSFORMS above).
% logit(pBurst) spreads the proportion away from [0,1] bounds.
% log(ss_frSingle) handles the log-normal distribution.
% Both transforms are applied to the DATA so the fit line stays straight.
tbl.log_ss_frSingle = log(tbl.ss_frSingle);

% --- Panel 1: RAW (Unadjusted) ---
% No model, no covariate adjustment. The slope here reflects the MARGINAL
% relationship, which conflates the direct effect of pBurst with indirect
% effects through correlated covariates (ss_frBurst). Axes are display-
% transformed (logit X, log Y) for readability, with the data transformed
% before fitting so the fit line and display agree.
hAx = nexttile;
plot_scat(tbl, 'pBurst_trans', 'log_ss_frSingle', ...
    'g', tbl.genotype, 'c', cfg.clr.grp, ...
    'fitType', 'linear', 'alpha', 0.4, 'sz', 20, ...
    'hAx', hAx, 'flgStats', true);
xlabel(hAx, 'pBurst (logit)', 'Interpreter', 'none');
ylabel(hAx, 'ss\_frSingle (log)', 'Interpreter', 'none');
title(hAx, 'Raw (Unadjusted)');

% --- Panel 2: MARGINALIZATION (LSMEANS) ---
% Fix ss_frBurst = mean. Sweep pBurst across its range. Predict
% ss_frSingle for a hypothetical "average" neuron in each genotype.
% lme_lsmeans back-transforms from model space to natural units (Hz).
% Log Y axis is a display transform for readability — predictions are
% already in Hz. No logit on X: the model generates a smooth grid over
% pBurst values, and the curve naturally handles the nonlinear mapping.
hAx = nexttile;
lme_lsmeans(mdlDemo, {'pBurst', 'genotype'}, ...
    'transParams', infoDemo.transParams, 'hAx', hAx);
title(hAx, 'Marginalization (LSMEANS)');
set(hAx, 'YScale', 'log');

% --- Panel 3: RESIDUALIZATION — Partial Residual (CPR) ---
% Y = residuals from a reduced model that excludes pBurst. These residuals
% are in MODEL SPACE (log units for log-normal), centered near 0.
% X = logit(pBurst) for display. We compute residuals manually (same logic
% as lme_pr 'residual' mode) and plot against tbl.pBurst_trans (logit)
% instead of the model-space pBurst (auto-log10 by lme_analyse).
% This is a display choice — logit is natural for proportions. The
% residuals (Y-axis) remain in model space regardless of the X-axis.
tblMdl = mdlDemo.Variables;
assert(height(tblMdl) == height(tbl), ...
    'Row mismatch between tbl and mdl.Variables (NaN exclusions?)');
frmlRed = lme_frml2rmv(char(mdlDemo.Formula), 'pBurst');
mdlRed = lme_fit(tblMdl, frmlRed, 'dist', 'Normal');
residY_cpr = residuals(mdlRed, 'ResidualType', 'Raw');
hAx = nexttile;
plot_scat([], tbl.pBurst_trans, residY_cpr, ...
    'g', tblMdl.genotype, 'c', cfg.clr.grp, ...
    'fitType', 'linear', 'alpha', 0.4, 'sz', 20, ...
    'hAx', hAx, 'flgStats', true);
xlabel(hAx, 'pBurst (logit)', 'Interpreter', 'none');
ylabel(hAx, [mdlDemo.ResponseName ' | Reduced (Residuals)'], 'Interpreter', 'none');
title(hAx, 'Residualization (CPR)');

% --- Panel 4: RESIDUALIZATION — Partial Regression (AVP) ---
% Both axes are residualized (model space). Y = residuals(Y | covariates),
% X = residuals(pBurst | covariates). The residualized pBurst is no longer
% a proportion — it's a centered deviation around 0, so logit doesn't apply.
% Linear axes only. Slope = model coefficient (Beta).
% Note: lme_pr calls plot_lineEq, which forces symmetric axis limits based
% on the max absolute value across BOTH axes. Since Y-residuals span a
% wider range than X-residuals, the x-axis gets stretched far beyond the
% data. We override xlim afterward to fit the actual X range.
hAx = nexttile;
[tblAVP, ~] = lme_pr(mdlDemo, 'pBurst', ...
    'flgMode', 'regression', 'hAx', hAx, ...
    'varGrp', 'genotype', 'transParams', infoDemo.transParams);
title(hAx, 'Residualization (AVP)');
xRange = [min(tblAVP.ResidX), max(tblAVP.ResidX)];
pad = 0.1 * diff(xRange);
xlim(hAx, [xRange(1) - pad, xRange(2) + pad]);
