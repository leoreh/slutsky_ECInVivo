function stats = lme_parkTest(tbl, frml)
% LME_PARKTEST Performs the Park Test to determine the Mean-Variance relationship.
%
%   STATS = LME_PARKTEST(TBL, FRML) estimates a provisional Mixed-Effects Model
%   and performs the Modified Park Test (Manning & Mullahy, 2001) to
%   characterize the variance structure of the data. This helps facilitate
%   the choice between Linear (Gaussian) vs. Generalized (Gamma/log-normal)
%   models.
%
%   The function automatically inspects the response variable's skewness and
%   sparseness to select an appropriate provisional distribution (Normal vs.
%   Gamma) and handles zero-inflation via `tbl_transform` if required.
%
%   INPUTS:
%       tbl         - (table) Data table containing variables.
%       frml        - (char/string) Model formula (e.g., 'y ~ x + (1|g)').
%
%   OUTPUTS:
%       stats       - (struct) Results structure containing:
%                     .Lambda         : (double) Park parameter estimate.
%                     .Recommendation : (char) Suggested distribution family.
%                     .MdlProvisional : (object) The fitted provisional model.
%                     .MdlPark        : (object) The auxiliary variance model.
%
%   See also: LME_FIT, LME_FRML2VARS, TBL_TRANSFORM, LME_COMPAREDISTS

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addRequired(p, 'frml', @(x) ischar(x) || isstring(x));
parse(p, tbl, frml);


%% ========================================================================
%  PRE-PROCESSING
%  ========================================================================

% Extract response variable
[~, varResp] = lme_frml2vars(frml);

if ~ismember(varResp, tbl.Properties.VariableNames)
    error('Response variable "%s" not found in table.', varResp);
end

yRaw = tbl.(varResp);

% Check properties
nZeros = sum(yRaw == 0);
pctZeros = nZeros / length(yRaw);
skewVal = skewness(yRaw);

% Logic: If skewed and not zero-dominated, try Gamma
isSkewed = abs(skewVal) > 1;
isSparse = pctZeros >= 0.05;

if isSkewed && ~isSparse
    % --- Gamma Candidate ---
    dist = 'Gamma';
    fitMet = 'REMPL';

    % Handle Zeros (if < 5%)
    if nZeros > 0
        warning('Response "%s" has %.1f%% zeros. Applying offset.', ...
            varResp, pctZeros*100);
        tbl = tbl_transform(tbl, 'flg0', true, 'varsInc', {varResp}, ...
            'flgZ', false, 'verbose', false);

        % Update yRaw for consistency
        yRaw = tbl.(varResp);
    end

else
    % --- Normal Candidate ---
    dist = 'Normal';
    fitMet = 'REML';
end

%  PROVISIONAL FIT
mdlProv = lme_fit(tbl, frml, 'dist', dist, 'fitMethod', fitMet);


%% ========================================================================
%  PARK TEST (VARIANCE REGRESSION)
%  ========================================================================
%  Model: ln(residuals^2) = alpha + lambda * ln(y_hat)
%  We use a Gamma GLM with Log link on squared residuals to avoid retransformation bias.

% 1. Get Predictions (Conditional on Random Effects)
yHat = predict(mdlProv);

% 2. Calculate Squared Residuals
resRaw = yRaw - yHat;
resSq = resRaw.^2;

% 3. Filter for valid log inputs
isValid = (yHat > 0) & (resSq > 0);

if sum(isValid) < length(yRaw) * 0.9
    warning('%d%% of observations excluded from Park Test (non-positive preds/resid).', ...
        round((1 - mean(isValid))*100));
end

% 4. Auxiliary Regression
% Regress Squared Residuals on Log-Predictions
% E[ResSq] = exp(alpha + lambda * ln(yHat))
parkX = log(yHat(isValid));
tblPark = table(resSq(isValid), parkX, 'VariableNames', {'ResSq', 'LogPred'});

% Fit Gamma GLM (Link=Log)
mdlPark = fitglm(tblPark, 'ResSq ~ LogPred', ...
    'Distribution', 'Gamma', 'Link', 'log');

% Extract Lambda (Slope)
lambdaEst = mdlPark.Coefficients.Estimate(2);
lambdaSe  = mdlPark.Coefficients.SE(2);
lambdaCi  = coefCI(mdlPark);
lambdaCi  = lambdaCi(2, :);



%% ========================================================================
%  INTERPRETATION
%  ========================================================================

recStr = 'Inconclusive';

% Heuristic check
% Lambda ~ 0: Gaussian (Constant Variance)
% Lambda ~ 1: Poisson (Variance ~ Mean)
% Lambda ~ 2: Gamma (Variance ~ Mean^2)
% Lambda ~ 3: Inverse Gaussian (Variance ~ Mean^3)

if (lambdaCi(1) <= 0) && (lambdaCi(2) >= 0)
    recStr = 'Gaussian';
elseif (lambdaCi(1) <= 1) && (lambdaCi(2) >= 1)
    recStr = 'Poisson';
elseif (lambdaCi(1) <= 2) && (lambdaCi(2) >= 2)
    recStr = 'Gamma';
elseif (lambdaCi(1) <= 3) && (lambdaCi(2) >= 3)
    recStr = 'InverseGaussian';
end


%% ========================================================================
%  LOG-NORMAL CHECK (If Lambda ~ 2)
%  ========================================================================
%  If variable fits Gamma, verify against Log-Normal by checking heavy tails.

stats.KurtosisLog = NaN;
stats.minDF = NaN;

if strcmp(dist, 'Gamma')
    % 1. Fit Log-Normal Model (LME on log(y))
    % Use tbl_transform to ensure log applied correctly
    % Note: We use 'e' base for mathematical consistency with Log-Normal definition
    tblLog = tbl_transform(tbl, 'logBase', 'e', 'varsInc', {varResp}, ...
        'flgZ', false, 'verbose', false);

    mdlLog = lme_fit(tblLog, frml, 'dist', 'Normal');

    % 2. Check Residual Kurtosis
    resLog = residuals(mdlLog);
    kurtLog = kurtosis(resLog);
    stats.KurtosisLog = kurtLog;

    % 3. Check DF (Small Sample Size)
    [~, ~, statsFE] = fixedEffects(mdlLog, 'DFMethod', 'Satterthwaite');
    minDF = min(statsFE.DF);
    stats.minDF = minDF;

    % Decision Rule
    % Kurtosis > 3.5 (Heavy Tails) or DF < 30 (Small Sample) -> Log-Normal
    if kurtLog > 3.5 || minDF < 30
        recStr = 'Log-Normal';
    end

    stats.MdlLog = mdlLog;
end


%% ========================================================================
%  OUTPUT
%  ========================================================================

stats.Lambda = lambdaEst;
stats.SE = lambdaSe;
stats.CI = lambdaCi;
stats.Recommendation = recStr;
stats.MdlProvisional = mdlProv;
stats.MdlPark = mdlPark;

fprintf('Park Test: Lambda = %.2f [%.2f, %.2f] -> Recommend: %s\n', ...
    lambdaEst, lambdaCi(1), lambdaCi(2), recStr);

end


%% ========================================================================
%  NOTE: THEORETICAL FRAMEWORK
%  ========================================================================
%  The Park Test is a diagnostic procedure used to characterize the family
%  of the underlying distribution by empirically estimating the relationship
%  between the mean and the variance. In Generalized Linear Models (GLMs),
%  the variance is assumed to be a power function of the mean:
%  Var(y|x) = alpha * [E(y|x)]^lambda.
%
%  The value of the exponent 'lambda' uniquely identifies the distribution
%  family.
%  - lambda of 0 implies constant variance (Gaussian/Homoscedastic),
%    indicating that a standard Linear Model is appropriate.
%  - lambda of 1 implies the variance grows linearly with the mean
%    (Poisson).
%  - lambda of 2 implies the variance grows with the square of the mean
%    (Gamma)
%  - lambda of 3 implies a cubic relationship (Inverse Gaussian).
%
%  Determining this parameter is the objective arbiter in the "Log vs.
%  Gamma" debate; if lambda is approximately 2, the variance structure
%  supports a Gamma GLM (or Log-Normal model), whereas a lambda near 0
%  supports a Gaussian model on raw data.
%
%  This function implements the "Modified Park Test" recommended by Manning
%  and Mullahy (2001). The procedure involves three steps. First, a
%  provisional model (typically Gamma with a Log link) is fitted to generate
%  predicted means. Second, raw residuals are calculated. Third, the squared
%  residuals are regressed on the predicted means in log-space. Crucially,
%  this step utilizes a Gamma GLM for the variance regression itself, rather
%  than a simple OLS on log-squared residuals. This avoids the severe bias
%  that can arise from retransforming heavy-tailed residuals, providing a
%  robust estimate of lambda.
%
%  In the specific context of Linear Mixed Models (LME/GLME), this function
%  uses conditional predictions (incorporating the estimated random effects)
%  rather than marginal predictions. This ensures that the test evaluates
%  the structure of the observation-level residual variance (epsilon) after
%  the subject-level variance (b) has been accounted for, properly isolating
%  the noise properties of the measurement itself.
%
%  REFERENCES:
%  1. Manning, W. G., & Mullahy, J. (2001). Estimating log models: to
%     transform or not to transform? Journal of Health Economics, 20(4).
%  ========================================================================


%% ========================================================================
%  NOTE: GAMMA VS. LOG-NORMAL
%  ========================================================================
%  When the Park Test indicates a variance structure where variability
%  grows with the square of the mean (Lambda approx 2), the analyst faces a
%  choice between two competing models: the Log-Normal model (LMM on
%  log-transformed data) and the Gamma GLMM (on raw data with Log link).
%  Both models mathematically imply the same Mean-Variance relationship,
%  but they optimize for different statistical properties.
%
%  1. The Log-Normal Model (LMM on log(y)):
%     This approach transforms the data to stabilize variance. Its primary
%     strength is **Precision** (Efficiency). As noted by Manning & Mullahy
%     (2001), OLS-based methods are remarkably resilient to heavy-tailed
%     residuals (high kurtosis). If the log-scale residuals are
%     leptokurtotic (Kurtosis > 3), the Log-Normal model will often yield
%     significantly smaller standard errors than the Gamma model.
%     Furthermore, in the context of mixed models (LME), this approach
%     allows for finite-sample corrections (Satterthwaite/Kenward-Roger
%     degrees of freedom), providing more accurate Type I error control for
%     small sample sizes. However, its major weakness is **Bias**. If the
%     residuals on the log-scale remain heteroscedastic (i.e., variance
%     changes across conditions even after transformation), retransforming
%     predictions back to the raw scale introduces systemic bias unless
%     complex corrections are applied.
%
%  2. The Gamma Model (GLMM on y):
%     This approach models the raw mean directly using a specific variance
%     function. Its primary strength is **Consistency**. The Gamma estimator
%     provides unbiased estimates of the arithmetic mean regardless of the
%     underlying error structure or heteroscedasticity on the log-scale. It
%     avoids the "retransformation problem" entirely. However, its weakness
%     is **inefficiency** in the face of heavy tails; standard errors can
%     explode if outliers are frequent. Additionally, most GLMM software
%     relies on asymptotic Z-tests, which can be overly optimistic (inflated
%     significance) when the number of subjects is small (< 50).
%
%  DECISION RULE: If the Park test yields Lambda approx 2, check the
%  residuals of the log-transformed data. If the log-residuals are
%  heavy-tailed (Kurtosis > 3) or the sample size is small, favor the
%  **Log-Normal** model to maximize power and inference validity. If the
%  log-residuals are heteroscedastic (variance differs significantly
%  between groups), favor the **Gamma** %  model to ensure unbiased
%  estimation of the mean effects.
%  ========================================================================


%% ========================================================================
%  NOTE: PARK TEST VS. AIC
%  ========================================================================
%  When diagnostics (Park Test) and performance metrics (AIC) disagree,
%  the AIC should take precedence as it measures the overall quality of the
%  prediction for the entire dataset.
%
%  Example from firing rate data: Following the heuristic of Manning &
%  Mullahy (2001), high kurtosis (~11) typically favors the Log-Normal
%  model to improve estimator efficiency. However, the AIC comparison
%  overwhelmingly favored the Gamma model (Delta AIC > 1000).
%
%  This discrepancy is explained by the shape of the data distribution. The
%  firing rate histogram exhibits a monotonic decay from zero (Exponential
%  shape), which is a native shape for the Gamma family (when shape
%  parameter <= 1). In contrast, the Log-Normal distribution forces the
%  probability density to approach zero as the value approaches zero (a
%  "hump" shape). Transforming this specific monotonic data to the log
%  scale likely introduced severe distortion, fitting the "body" of the data
%  poorly despite stabilizing the variance of the "tails."
%  ========================================================================


%% ========================================================================
%  NOTE: PROVISIONAL VS. AUXILIARY MODEL DISTINCTIONS
%  ========================================================================
%  A common point of confusion is why the distribution of the auxiliary
%  Park regression (`mdlPark`) often differs from the provisional model
%  (`mdlProv`).
%
%  1. The Provisional Model (`mdlProv`):
%     This models the Mean function: E[y|X].
%     Its distribution (Normal vs. Gamma) effectively represents our "best
%     guess" about the shape of the raw data `y`. We select this based on
%     the skewness and zero-inflation of `y` to ensure valid predictions.
%
%  2. The Auxiliary Park Model (`mdlPark`):
%     This models the Variance function: E[epsilon^2 | X].
%     The dependent variable here is the Squared Residuals (`resSq`). By
%     mathematical definition, squared residuals are (a) strictly non-negative
%     and (b) typically highly right-skewed.
%     The Gamma distribution is the natural choice for modeling non-negative,
%     skewed continuous data. Manning & Mullahy (2001) demonstrate that
%     using a Gamma GLM (with Log link) to model squared residuals provides
%     a robust estimate of the variance scaling parameter (lambda).
%
%  Therefore, it is entirely correct and expected for `mdlProv` to be
%  Gaussian (modeling a symmetric mean) while `mdlPark` is Gamma (modeling
%  skewed squared error terms). Validating the variance structure does not
%  require the mean and variance models to share the same distribution family.
%  ========================================================================