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

    figure('Color', 'w', 'Name', 'mea_lmeCens', 'Position', [50 50 1400 500]);
    tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    % --- PLOT 1: Transformed Space (Fit Comparison) ---
    nexttile; hold on;
    sObs = scatter(mdl.Fitted(~maskCens), yTrans(~maskCens), 20, 'k', 'filled', ...
        'MarkerFaceAlpha', 0.3, 'DisplayName', 'Obs');
    sCens = scatter(mdl.Fitted(maskCens), yTrans(maskCens), 20, 'r', 'filled', ...
        'MarkerFaceAlpha', 0.3, 'DisplayName', 'Cens');
    sImp = scatter(mdl.Fitted(maskCens), yPredTrans(maskCens), 20, 'b', 'filled', ...
        'MarkerFaceAlpha', 0.4, 'DisplayName', 'Imp');

    refline(1, 0);
    hY = yline(thrCens, 'k--', 'DisplayName', 'Thr');

    ylabel('Response (Transformed)');
    xlabel('Fitted (Latent)');
    title('Tobit Fit (Transformed)');
    legend([sObs, sCens, sImp, hY], 'Location', 'best');
    grid on; axis square;

    % --- PLOT 2: Original Space (The Solution) ---
    nexttile; hold on;
    
    % Raw data (response and first predictor)
    yRaw = tbl.(varResp);
    varsPred = lme_frml2vars(frml);
    xName = varsPred{1}; 
    xVals = tbl.(xName);

    % Scatters
    sObs = scatter(xVals(~maskCens), yRaw(~maskCens), 25, 'k', 'filled', ...
        'MarkerFaceAlpha', 0.4, 'DisplayName', 'Obs');
    sCens = scatter(xVals(maskCens), yRaw(maskCens), 25, 'r', 'filled', ...
        'MarkerFaceAlpha', 0.4, 'DisplayName', 'Cens');
    sImp = scatter(xVals(maskCens), yPred(maskCens), 20, 'b', 'filled', ...
        'MarkerFaceAlpha', 0.4, 'DisplayName', 'Imp');
    
    % Regression Lines
    % - OLS on Observed Only
    % - OLS on Original (Clamped)
    % - OLS on Imputed
    mdlObs = fitlm(xVals(~maskCens), yRaw(~maskCens));
    mdlOrig = fitlm(xVals, yRaw);
    mdlImp = fitlm(xVals, yPred);

    % To plot lines properly, we predict on a sorted range
    xRange = linspace(min(xVals), max(xVals), 100)';

    % Predict OLS Obs
    yLine = predict(mdlObs, xRange);
    p1 = plot(xRange, yLine, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Fit: Obs');

    % Predict OLS All
    yLine = predict(mdlOrig, xRange);
    p2 = plot(xRange, yLine, 'r:', 'LineWidth', 1.5, 'DisplayName', 'Fit: Orig');

    % Predict Tobit (Re-transform)
    yLine = predict(mdlImp, xRange);
    p3 = plot(xRange, yLine, 'b:', 'LineWidth', 1.5, 'DisplayName', 'Fit: Imp');

    xlabel(xName);
    ylabel(varResp);
    legend([sObs, sCens, sImp, p1, p2, p3], 'Location', 'best', 'NumColumns', 2);
    title('Censored Imputation');
    grid on; axis square;

    % --- PLOT 3: Distributions (Log-Log) ---
    nexttile; hold on;

    % Determine log-scale
    edges = 50;
    if contains(lower(dist), 'log')
        edgeMin = (log10(min(yPred)));
        edgeMax = (log10(max(yPred)));
        edges = logspace(edgeMin, edgeMax, edges);
        set(gca, 'XScale', 'log');
    end

    histogram(yRaw, edges, 'Normalization', 'pdf', 'DisplayStyle', 'bar', ...
        'EdgeColor', 'none', 'FaceColor', 'r', 'FaceAlpha', 0.2, 'LineWidth', 2, 'DisplayName', 'Raw');
    histogram(yRaw, edges, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', ...
        'EdgeColor', 'r', 'LineWidth', 1, 'HandleVisibility', 'off');

    histogram(yPred, edges, 'Normalization', 'pdf', 'DisplayStyle', 'bar', ...
        'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.2, 'LineWidth', 2, 'DisplayName', 'Raw');
    histogram(yPred, edges, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', ...
        'EdgeColor', 'b', 'LineWidth', 1, 'HandleVisibility', 'off');

    xlim([edgeMin, edgeMax])
    xlabel(varResp);
    ylabel('PDF');
    title('Distribution Reconstruction');
    legend('Location', 'best');
    grid on; axis square;
    sgtitle(['Tobit Regression: ' frml], 'Interpreter', 'none');
end

end         % EOF



%% ========================================================================
%  NOTE: BACKGROUND
%  ========================================================================
% PROBLEM DEFINITION: CENSORING ARTIFACTS IN PERTURBATION ANALYSIS
%
%   In MEA experiments involving strong inhibitors (e.g., Baclofen), a
%   subset of neurons often falls completely silent. Because the logarithm
%   of zero is undefined, standard analysis pipelines clamp the firing rate
%   of these silenced units to a small bias constant (c, typically 1
%   spike/hour).
%
%   This clamping introduces a systematic artifact visible in scatter plots
%   (eg, Perturbation Depth vs. Baseline Firing Rate). Silenced units form
%   a sharp diagonal line because their calculated depth becomes a linear
%   function of their baseline rate (Depth = log2(Baseline) - log2(c)).
%
%   Consequently, standard statistical approaches (t-tests, standard linear
%   regression) are invalid for two reasons:
%   1.  Bias against high-firing units: The artificial floor underestimates
%       the true suppression of high-baseline neurons, effectively penalizing
%       them for having a wider dynamic range.
%   2.  Violation of normality: The data distribution is not continuous; it
%       contains a spike of values at the clamping limit, violating the
%       assumptions of OLS regression.
%
% THEORETICAL FRAMEWORK: TOBIT MODEL AND LATENT VARIABLES
%
%   The correct statistical framework for this data is the Tobit (Censored)
%   Model. This model posits that there exists a latent variable, the
%   "true" biological firing rate, which continues continuously below the
%   detection threshold of the recording equipment.
%
%   For active units, the observed rate equals the true rate. For silenced
%   units, the observed rate is merely an upper bound; the true biological
%   propensity for firing is strictly lower than the bias constant c. This
%   is mathematically defined as Left Censoring.
%
%   By fitting a censored regression model, we estimate the parameters of
%   the entire population (both active and silenced) simultaneously. The
%   model reconstructs the linear relationship between Baseline and Trough
%   rates that would have been observed had the equipment possessed
%   infinite sensitivity.
%
% METHODOLOGY: CENSORING-CORRECTED IMPUTATION
%
%   To resolve the artifact and enable standard downstream analysis (like
%   averaging or hypothesis testing), we employ Tobit Imputation. This
%   process replaces the clamped artificial values with statistically
%   estimated values derived from the conditional expectation of the
%   censored tail.
%
%   The imputation process follows these sequential steps:
%
%   1.  Log-Transformation
%       Convert both Baseline and Trough firing rates to log-space (e.g.,
%       log10). This linearizes the multiplicative biological inhibition
%       and stabilizes the variance (homoscedasticity).
%
%   2.  Identification of Censored Units
%       Create a logical vector identifying units where the Trough rate is
%       equal to or lower than the bias constant c. These are the Left
%       Censored observations.
%
%   3.  Fitting the Censored Model
%       Fit a Left-Censored Linear Regression model to the data using the
%       equation: log(Trough) ~ log(Baseline).
%       This model utilizes the exact values of active units and the count of
%       silenced units to estimate:
%       -   The Slope and Intercept (The population trend).
%       -   The Sigma (The standard deviation of the residuals).
%
%   4.  Calculation of Conditional Expectation
%       For every silenced unit, calculate its imputed value. We do not
%       simply use the mean prediction from the regression line, as we know
%       the unit fell below the threshold.
%
%       The imputed value is the Predicted Mean minus a correction term.
%       The correction term is the model Sigma multiplied by the Inverse
%       Mills Ratio (probability density divided by cumulative probability)
%       evaluated at the censoring limit.
%
%       This places the imputed value correctly in the "tail" of the
%       distribution, statistically distancing it from the detection limit
%       based on how likely that specific unit was to be silenced.
%
%   5.  Data Replacement
%       Replace the raw log(Trough) values of the silenced units with these
%       new imputed values. Active units retain their original data.
%
% PRACTICAL OUTCOME AND INTERPRETATION
%
%   After imputation, the diagonal artifact in the scatter plot will
%   disperse. Silenced units will drift upward (indicating deeper
%   suppression), with high-baseline units shifting significantly more than
%   low-baseline units.
%  ========================================================================

%% ========================================================================
%  NOTE: PREDICT VS. CONDITIONAL EXPECTATION
%  ========================================================================
%   Why do we manually calculate imputation values instead of using the
%   standard 'predict' method of the 'fitlmcens' object?
%
%   1.  THE LIMITATION OF STANDARD PREDICTION
%       The standard `predict` function returns the unconditional mean of
%       the latent variable (the linear predictor $X\beta + Zb$). This
%       value represents where a neuron "would be" on average based on the
%       population parameters, ignoring the specific knowledge that this
%       particular observation was censored.
%
%       For a silenced neuron, the unconditional mean might actually lie
%       above the detection threshold (especially if the model has high
%       variance). Using this value would contradict the observed reality
%       that the neuron did not fire.
%
%   2.  THE NECESSITY OF CONDITIONAL EXPECTATION
%       Our goal is "Data Rehabilitation"—reconstructing the most likely
%       value for a specific silenced unit given that we know it fell below
%       the threshold $L$. This requires the Conditional Expectation:
%
%       $$ E[y \mid y < L] = \mu - \sigma \cdot \lambda(\alpha) $$
%
%       Where:
%       -   $\mu$: The unconditional mean ($X\beta$).
%       -   $\sigma$: The scale parameter of the model.
%       -   $\lambda(\alpha)$: The Inverse Mills Ratio, defined as
%           $\frac{\phi(\alpha)}{\Phi(\alpha)}$.
%       -   $\alpha$: The standardized threshold distance, $\frac{L - \mu}{\sigma}$.
%
%   3.  IMPLEMENTATION DETAIL
%       The manual calculation in this script (computing the Inverse Mills
%       Ratio explicitly) ensures that the imputed value is statistically
%       pushed into the tail of the distribution. This correction:
%       -   Respects the upper bound $L$ (the detection limit).
%       -   Accounts for the probability mass of the "missing" tail.
%       -   Provides a more accurate reconstruction of the silenced phenotype
%           than the population average.
%  ========================================================================


%% ========================================================================
%  NOTE: POOLING VS. PER-FILE
%  ========================================================================
%   1.  ABSENCE OF MIXED EFFECTS IN CENSORED MODELS
%       MATLAB's `fitlmcens` function fits a Censored Linear Regression
%       model using Maximum Likelihood Estimation. Unlike standard LME
%       tools (`fitlme`), `fitlmcens` does not support Random Effects.
%       Consequently, it is not possible to include a random intercept term
%       (e.g., `(1|FileID)`) to explicitly model batch-to-batch variability
%       within a single global model.
%
%   2.  IMPLICATION FOR DATA REHABILITATION
%       Because the model cannot mathematically account for the correlation
%       of units within the same recording file via random effects, there
%       are two processing strategies:
%
%       A.  Pooled Fitting (All Files)
%           Fitting the model to the aggregated table of all files.
%           This assumes a single homogeneous population. If baseline
%           firing rates or detection thresholds vary systematically
%           between recording sessions, pooling may yield an averaged slope
%           that fits no single file well (Simpson's Paradox), leading to
%           inaccurate imputation for specific batches.
%
%       B.  Per-File Fitting (Sequential)
%           Applying `mea_lmcens` independently to each recording file.
%           This implicitly treats 'FileID' as a fixed effect. It
%           allows the regression parameters (Intercept and Slope) and the
%           noise parameter (Sigma) to adapt to the specific signal-to-noise
%           ratio and baseline distribution of that specific experiment.
%
%   3.  RECOMMENDATION
%       For the purpose of calculating Perturbation Depth (a relative
%       metric), preserving the within-file relationship between Baseline
%       and Trough rates is critical. Therefore, unless the number of units
%       per file is too low for convergence (<10-15), the Per-File fitting
%       strategy is statistically superior. It ensures that the
%       rehabilitation of silenced units is calibrated to the local
%       experimental conditions rather than a global average.
%  ========================================================================


%% ========================================================================
%  NOTE: ANALYSIS EFFICACY
%  ========================================================================
%  OBSERVATION: MARGINAL UTILITY OF SOPHISTICATED IMPUTATION
%
%   Comparison of the standard "Clamping" approach against the sophisticated
%   "Tobit Censored Imputation" reveals negligible differences in the final
%   data distribution. As seen in the diagnostic plots, the imputed values
%   (Blue) deviate minimally from the censored/clamped threshold (Red). The
%   reconstructed distribution remains a "spike" near the detection limit rather
%   than expanding into a Gaussian tail.
%
%  DIAGNOSIS: PREDICTION ALIGNMENT
%
%   The magnitude of the Tobit correction depends on the distance between the
%   detection threshold (L) and the model's predicted latent mean (Mu).
%   Imputation is most aggressive when the model predicts a unit "should" be
%   firing well below the threshold based on population trends.
%
%   In this dataset, the censored units are exclusively low-baseline neurons.
%   The regression line (Population Slope) correctly predicts that low-baseline
%   units should have trough rates near the detection limit. Because the
%   Predicted Mean (Mu) is effectively equal to the Threshold (L), the
%   Conditional Expectation term approaches zero.
%
%   The model has statistically determined that these units are not "deeply
%   suppressed" into a theoretical abyss, but are simply "just barely" silent.
%   Therefore, the simple clamp is a statistically accurate approximation of
%   the latent reality.
%
%  CONCLUSION: YAK SHAVING
%
%   The implementation of `mea_lmCens` is mathematically correct but empirically
%   unnecessary for this specific biological phenotype. The complexity of
%   maintaining the Tobit pipeline outweighs the <0.1 log-unit precision gain.
%
%   Recommended Actions:
%   1.  Revert to standard bias-correction (Clamping to constant 'c').
%   2.  To improve resolution without statistical inference, increase the
%       integration window width for the Trough period in `mea_frRcv`.
%  ========================================================================