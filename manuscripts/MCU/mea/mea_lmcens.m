function [yPred, info] = mea_lmcens(tbl, frml, varargin)
% MEA_LMCENS Censored Regression (Tobit) Imputation for MEA Data.
%
%   [YPRED, INFO] = MEA_LMCENS(TBL, FRML, ...) performs Tobit censored
%   regression to impute values for units that fall below a detection
%   threshold (Left Censoring).
%
%   Process:
%       1. Prepares table (Log/Z-Score) using TBL_TRANS.
%       2. Fits Tobit Model (fitlmcens).
%       3. Imputes censored values using conditional expectation.
%
%   INPUTS:
%       tbl         - (table) Data table.
%       frml        - (char/string) Model formula.
%       varargin    - (param/value)
%                     'dist'               : (char) 'log-normal' (default) or 'Normal'.
%                     'flgPlot'            : (logical) Plot diagnostics.
%
%   OUTPUTS:
%       yPred       - (double) Rehabilitated response vector (Original scale if dist='Normal').
%                     If dist='log-normal', yPred is in Log-Space.
%       info        - (struct) Model object and stats.
%
%   See also: FITLMCENS, TBL_TRANS

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addRequired(p, 'frml', @(x) ischar(x) || isstring(x));
addParameter(p, 'thrCens', [], @isnumeric);
addParameter(p, 'dist', 'log-normal', @ischar);
addParameter(p, 'flgPlot', true, @islogical);

parse(p, tbl, frml, varargin{:});

thrCens = p.Results.thrCens;
dist = p.Results.dist;
flgPlot = p.Results.flgPlot;
frml = char(frml);

%% ========================================================================
%  FITTING
%  ========================================================================

% FORMULA
[~, varResp, ~, ~] = lme_frml2vars(frml);

% TRANSFORM
[~, ~, tblInfo, tblMdl] = lme_analyse(tbl, frml, 'dist', dist);

y = tblMdl.(varResp);
censMask = tbl.censMask;
fprintf('[MEA_LMCENS] Fitting Tobit (Censored at %.4f, N=%d, %.1f%% Censored)\n', ...
    thrCens, length(y), 100*mean(censMask));

% FIT
mdl = fitlmcens(tblMdl, varResp, frml, "Censoring", censMask); 


%% ========================================================================
%  IMPUTATION
%  ========================================================================

% Predict Latent Mean (Mu)
mu = predict(mdl, tblMdl);

sigma = mdl.Sigma;

thrCens = min(tblMdl.(varResp));

% Impute Censored Values
% E[y | y < L] = mu - sigma * (phi(alpha) / Phi(alpha))
% alpha = (L - mu) / sigma
yPred = y;

idxCen = censMask;
if any(idxCen)
    muCen = mu(idxCen);
    alpha = (thrCens - muCen) / sigma;

    pdfA = normpdf(alpha);
    cdfA = normcdf(alpha);
    cdfA(cdfA < realmin) = realmin;

    lambda = pdfA ./ cdfA;
    yPred(idxCen) = muCen - sigma * lambda;
end

% Construct Info
info.mdl = mdl;
info.transParams = tblInfo.transParams;
info.thrCens = thrCens;

%% ========================================================================
%  PLOTTING
%  ========================================================================

if flgPlot
    figure('Color', 'w', 'Name', 'mea_lmeCens');

    subplot(1,2,1); hold on;
    scatter(mu(~idxCen), y(~idxCen), 15, 'k', 'filled', 'DisplayName', 'Observed');
    scatter(mu(idxCen), y(idxCen), 15, 'r', 'DisplayName', 'Censored');
    scatter(mu(idxCen), yPred(idxCen), 15, 'b', 'filled', 'DisplayName', 'Imputed');
    yline(thrCens, '--k', 'DisplayName', 'Threshold');
    refline(1,0);
    xlabel('Predicted Latent (X\beta)'); ylabel('Response');
    title('Tobit Imputation'); legend('Location','best'); grid on;

    subplot(1,2,2); hold on;
    histogram(y, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2, 'DisplayName', 'Original');
    histogram(yPred, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 2, 'DisplayName', 'Rehabilitated');
    title('Distribution'); legend('Location','best'); grid on;
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
