function [yPred, info] = mea_lmCens(tbl, frml, varargin)
% MEA_LMCENS Censored Regression (Tobit) Imputation for MEA Data.
%
%   [YPRED, INFO] = MEA_LMCENS(TBL, FRML, ...) performs Tobit censored
%   regression to impute values for units that fall below a detection
%   threshold (Left Censoring).
%
%   Process:
%       1. Prepares table (Log/Z-Score) using LME_ANALYSE -> TBL_TRANS.
%       2. Fits Tobit Model (FITLMCENS) using provided censoring mask.
%       3. Imputes censored values using conditional expectation in Transformed Space.
%       4. Back-transforms predictions to Original Space.
%
%   INPUTS:
%       tbl         - (table) Data table.
%       frml        - (char/string) Model formula.
%       varargin    - (param/value)
%                     'censVar'            : (char) Name of logical variable in TBL indicating censored observations.
%                     'dist'               : (char) 'log-normal' (default) or 'Normal'.
%                     'flgPlot'            : (logical) Plot diagnostics.
%
%   OUTPUTS:
%       yPred       - (double) Rehabilitated response vector (Original scale).
%       info        - (struct) Model object and stats.
%
%   See also: FITLMCENS, TBL_TRANS, MEA_FRRCV

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addRequired(p, 'frml', @(x) ischar(x) || isstring(x));
addParameter(p, 'censVar', '', @ischar);
addParameter(p, 'dist', 'log-normal', @ischar);
addParameter(p, 'flgPlot', true, @islogical);

parse(p, tbl, frml, varargin{:});

censVar = p.Results.censVar;
dist = p.Results.dist;
flgPlot = p.Results.flgPlot;
frml = char(frml);

if isempty(censVar) || ~ismember(censVar, tbl.Properties.VariableNames)
    error('Must provide valid "censVar" (logical column name) for censoring.');
end

%% ========================================================================
%  FITTING
%  ========================================================================

% FORMULA
[~, varResp, ~, ~] = lme_frml2vars(frml);

% EXTRACT CENSORING MASK (Raw)
% Note: We extract this before transformation to ensure alignment.
% Assumes lme_analyse does not drop rows (no NaNs in used variables).
maskCens =  logical(tbl.(censVar));

% TRANSFORM
% We use lme_analyse to leverage the standard transformation pipeline.
% : This fits a standard LME as a side-effect, which we ignore.
[~, ~, lmeInfo, tblMdl] = lme_analyse(tbl, frml, 'dist', dist, 'verbose', false);

yTrans = tblMdl.(varResp);

% Check Threshold Uniformity in Transformed Space
% (Ideally, all censored values should be identical in Log-Space)
valsCens = yTrans(maskCens);
if std(valsCens) > 1e-6
    warning('Transformed censored values are not identical. Check clamping logic.');
end
thrCens = median(valsCens); % Robust estimate of L in transformed space

% GIT
fprintf('[MEA_LMCENS] Fitting Tobit (N=%d, %.1f%% Censored)\n', ...
    length(yTrans), 100*mean(maskCens));

% Fit Model
mdl = fitlmcens(tblMdl, varResp, frml, "Censoring", maskCens);

%% ========================================================================
%  IMPUTATION (Transformed Space)
%  ========================================================================

% Predict Latent Mean (Mu)
mu = predict(mdl, tblMdl);

sigma = mdl.Sigma;

% Impute Censored Values
% E[y | y < L] = mu - sigma * (phi(alpha) / Phi(alpha))
% alpha = (L - mu) / sigma
yPredTrans = yTrans;

idxCen = maskCens;
if any(idxCen)
    muCen = mu(idxCen);
    alpha = (thrCens - muCen) / sigma;

    pdfA = normpdf(alpha);
    cdfA = normcdf(alpha);
    cdfA(cdfA < realmin) = realmin;

    lambda = pdfA ./ cdfA;
    yPredTrans(idxCen) = muCen - sigma * lambda;
end

% Construct Info
info.mdl = mdl;
info.transParams = lmeInfo.transParams;
info.thrCens = thrCens;

%% ========================================================================
%  INVERSE TRANSFORMATION (Back to Original Space)
%  ========================================================================

% Use tbl_trans in INVERSE mode
% This undoes Log/Logit/Z-score on the response and predictors
% We only care about the response.
% Create a temp table with the IMPUTED response for back-transformation
tblMdl.(varResp) = yPredTrans;
tblOut = tbl_trans(tblMdl, 'template', lmeInfo.transParams, ...
    'flgInv', true, 'verbose', false);
yPred = tblOut.(varResp);

%% ========================================================================
%  PLOTTING
%  ========================================================================

% Calculate Back-Transformed Predicted Latent (Mu) for Plotting
tblMdl.(varResp) = mu;
tblOutMu = tbl_trans(tblMdl, 'template', lmeInfo.transParams, 'flgInv', true, 'verbose', false);
muOrig = tblOutMu.(varResp);

%% ========================================================================
%  PLOTTING
%  ========================================================================

if flgPlot
    yRaw = tbl.(varResp);
    figure('Color', 'w', 'Name', 'mea_lmeCens', 'Position', [100 100 1400 400]);

    % Plot 1: Transformed Space Fit
    subplot(1,3,1); hold on;
    scatter(mu(~idxCen), yTrans(~idxCen), 15, 'k', 'filled', 'DisplayName', 'Observed');
    scatter(mu(idxCen), yTrans(idxCen), 15, 'r', 'DisplayName', 'Censored (Clamped)');
    scatter(mu(idxCen), yPredTrans(idxCen), 15, 'b', 'filled', 'DisplayName', 'Imputed (Exp)');
    yline(thrCens, '--k', 'DisplayName', 'Threshold');
    refline(1,0);
    xlabel('Predicted Latent (Transformed)'); ylabel('Response (Transformed)');
    title('Tobit Imputation (Log-Space)'); legend('Location','best'); grid on;

    % Plot 2: Original Space Fit (Back-Transformed)
    subplot(1,3,2); hold on;
    scatter(muOrig(~idxCen), yPred(~idxCen), 15, 'k', 'filled', 'DisplayName', 'Observed');
    scatter(muOrig(idxCen), yRaw(idxCen), 15, 'r', 'DisplayName', 'Censored (Clamped)');
    scatter(muOrig(idxCen), yPred(idxCen), 15, 'b', 'filled', 'DisplayName', 'Imputed (Exp)');

    % Approx Threshold line in original space (Back-transform threshold?)
    % Since thrCens is median in log space, we can roughly estimate:
    if length(unique(yRaw(maskCens))) == 1
        yline(unique(yRaw(maskCens)), '--k', 'DisplayName', 'Min Limit');
    end

    refline(1,0);
    xlabel('Predicted Latent (Original)'); ylabel('Response (Original)');
    title('Rehabilitation (Original Scale)'); legend('Location','best'); grid on;

    % Plot 3: Final Distribution (Original Space)
    subplot(1,3,3); hold on;
    histogram(yRaw, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2, 'DisplayName', 'Original (Zero/Clamped)');
    histogram(yPred, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2, 'DisplayName', 'Rehabilitated');
    title('Distribution (Original Scale)'); legend('Location','best'); grid on;
end

end         % EOF
