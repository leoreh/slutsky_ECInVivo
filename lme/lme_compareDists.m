function stats = lme_compareDists(tbl, frml)
% LME_COMPAREDISTS Compares GLME fits across multiple distributions.
%
%   stats = LME_COMPAREDISTS(TBL, FRML) fits Generalized Linear
%   mixed models using five common distributions (Normal, Log-Normal,
%   Logit-Normal, Poisson, Gamma, and InverseGaussian) and returns a table of fit statistics.
%
%   INPUTS:
%       tbl         - (table) Table containing the variables.
%       frml        - (char/formula) Model definition.
%
%   OUTPUTS:
%       stats       - (table) Summary of model performance (AIC, BIC, etc.).
%

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addRequired(p, 'frml', @(x) ischar(x) || isa(x, 'classreg.regr.LinearFormula'));
parse(p, tbl, frml);


%% ========================================================================
%  INITIALIZE
%  ========================================================================

% List of distributions to test
mdls = {'Normal', 'Log-Normal', 'Logit-Normal', 'Poisson', 'Gamma', 'InverseGaussian'};
nMdl = length(mdls);

% Common canonical/standard links
% Normal -> Identity
% Log/Logit-Normal -> Identity (on transformed data)
% Others -> Log (Enforced as per best practice for positive biological data)
links = {'Identity', 'Identity', 'Identity', 'Log', 'Log', 'Log'};

% Output table
stats = table();
stats.Model = mdls(:);
stats.Link = links(:);
stats.AIC = nan(nMdl, 1);
stats.BIC = nan(nMdl, 1);
stats.LogLikelihood = nan(nMdl, 1);
stats.R2_Ordinary = nan(nMdl, 1);
stats.R2_Adjusted = nan(nMdl, 1);
stats.Converged = false(nMdl, 1);
stats.ErrorMsg = repmat({''}, nMdl, 1);

% Hardcoded for Model Comparison
fitMethod = 'Laplace';

% Parse response variable
varResp = strtrim(strsplit(frml, '~'));
varResp = varResp{1};
yRaw = tbl.(varResp);           % Store raw values before transformations

%% ========================================================================
%  VALIDATION
%  ========================================================================

% Check for non-positive values
nZeros = sum(tbl.(varResp) == 0);
pctZeros = nZeros / length(tbl.(varResp));
isNeg = any(tbl.(varResp) < 0);

if isNeg
    error(['Response variable "%s" contains negative values.', ...
        'Cannot fit Log-Normal, Gamma, or Poisson.'], varResp);
end

% Assert non-zeros
if nZeros > 0
    warning(['Response variable "%s" has %.1f%% zeros.', ...
        'Results may be unstable.'], varResp, pctZeros*100);

    tbl = tbl_transform(tbl, 'flg0', true, 'verbose', false, ...
        'varsInc', {varResp});
end



%% ========================================================================
%  RUN FITS
%  ========================================================================

for iMdl = 1:nMdl

    % loop variables
    mdl = stats.Model{iMdl};
    dist = mdl;
    link = stats.Link{iMdl};
    mdlTbl = tbl;
    mdlFrml = frml;
    jcb = 0;

    try

        % TRANSFORMATIONS
        if strcmp(dist, 'Log-Normal')

            % Ensure Natural Log for comparison
            mdlTbl = tbl_transform(tbl, 'flgLog', true, 'logBase', 'e', ...
                'verbose', false, 'skewThr', -Inf);
            dist = 'Normal';

            % Jacobian Correction: LL_raw = LL_ln - sum(ln(y))
            % Use original y values (stored in yRaw, or retrieve from tbl)
            jcb = -sum(log(yRaw), 'omitnan');

        elseif strcmp(dist, 'Logit-Normal')

            mdlTbl = tbl_transform(tbl, 'flgLogit', true, 'verbose', false);
            dist = 'Normal';

            % Jacobian Correction: LL_raw = LL_logit - sum(ln(y*(1-y)))
            % Clip for stability (consistent with tbl_transform)
            yVal = max(eps, min(1-eps, yRaw));
            jcb = -sum(log(yVal .* (1 - yVal)), 'omitnan');
        end


        % RUN FIT
        mdl = fitglme(mdlTbl, mdlFrml, ...
            'Distribution', dist, ...
            'Link', link, ...
            'FitMethod', fitMethod);


        % STATISTICS

        % Raw stats from model
        ll = mdl.LogLikelihood;
        aic = mdl.ModelCriterion.AIC;
        bic = mdl.ModelCriterion.BIC;

        % Apply Jacobian Correction (converts back to raw scale)
        if jcb ~= 0

            % We need k (num params) to recalculate AIC/BIC
            % AIC_mdl = -2*LL_mdl + 2*k  =>  2*k = AIC_mdl + 2*LL_mdl
            k2 = aic + 2*ll;

            ll = ll + jcb; % Corrected LogLikelihood

            nObs = mdl.NumObservations;

            aic = -2*ll + k2;
            bic = -2*ll + (k2/2)*log(nObs);
        end


        % Store
        stats.AIC(iMdl) = aic;
        stats.BIC(iMdl) = bic;
        stats.LogLikelihood(iMdl) = ll;
        stats.R2_Ordinary(iMdl) = mdl.Rsquared.Ordinary;
        stats.R2_Adjusted(iMdl) = mdl.Rsquared.Adjusted;
        stats.Converged(iMdl) = true;

    catch ME
        stats.ErrorMsg{iMdl} = ME.message;
    end
end

% Sort by AIC (Low to High) - Best model first
[~, sortIdx] = sort(stats.AIC);
stats = stats(sortIdx, :);

end     % EOF


%% ========================================================================
%  NOTE: FIT STATS
%  ========================================================================
%
%  LogLikelihood (LogLik):
%     The log-probability of the data given the model. Higher (less negative)
%     is better. However, it always increases with more parameters, so it
%     cannot be used alone to compare models with different complexities.
%
%  Information Criteria (AIC & BIC):
%     These balance model fit (LogLik) against complexity (number of
%     parameters). LOWER is better.
%     - AIC (Akaike): Optimizes for predictive accuracy.
%     - BIC (Bayesian): Penalizes parameters more strongly; optimizes for
%       identifying the "true" model.
%     For model selection, a difference > 2 is positive evidence, and > 10
%     is very strong evidence.
%
%  R-Squared (R2):
%     Represents the proportion of variance explained by the model.
%     - Ordinary: Variance explained by Fixed + Random effects.
%     - Adjusted: Adjusted for number of coefficients (penalizes overfitting).
%     NOTE: In GLMEs (non-Normal), R2 is a pseudo-statistic and should be
%     interpreted with caution compared to OLS R2.
%  ========================================================================


%  ========================================================================
%  NOTE: MATLAB'S 'COMPARE'
%  ========================================================================
%  MATLAB's built-in `compare(mdl1, mdl2)` function performs a Likelihood
%  Ratio Test (LRT).
%
%  - WHEN TO USE 'COMPARE':
%    Only when models are NESTED. That is, Mdl2 is Mdl1 plus extra terms
%    (e.g., adding an interaction). They must share the same Distribution
%    and Link.
%
%  - WHEN TO USE THIS FUNCTION ('LME_COMPAREDISTS'):
%    When comparing DIFFERENT distributions (e.g., Gamma vs InverseGaussian).
%    These models are NOT nested. Therefore, you cannot use an LRT/p-value.
%    You must use Information Criteria (AIC/BIC) to select the distribution
%    that best balances fit and parsimony.
%  ========================================================================


%% ========================================================================
%  NOTE: LINK FUNCTIONS
%  ========================================================================
%  A GLM specifies the relationship between the Mean Response (mu) and the
%  linear predictor (X*beta) via a Link Function g(mu) = X*beta.
%
%  While statistical software defaults to "canonical" links for
%  computational simplicity—such as the Reciprocal link for Gamma
%  distributions - these defaults often contradict the physical reality of
%  biological data.

%  The Canonical links for Gamma/InvGauss imply that as X increases, Y
%  decreases hyperbolically (1/X). This is rarely the true biological
%  mechanism. Most positive data (like RT) follows a multiplicative process
%  (e.g., "a 10% increase"), which implies an Exponential trend. This
%  corresponds to the LOG LINK: log(mu) = X*beta.

%  Therefore, the Log link is the scientifically appropriate choice for
%  Gamma, Inverse Gaussian, and Poisson distributions as it models these
%  exponential trends. In contrast, the Normal distribution uses the
%  Identity link by default, representing the standard assumption of
%  additive, linear effects. By enforcing the Log link across the
%  heavy-tailed distributions while keeping Identity for the Normal
%  distribution, we ensure that the AIC comparison properly tests the
%  variance structure (distribution) rather than inadvertently penalizing a
%  plausible distribution for using an implausible hyperbolic trend.
%  ========================================================================


%% ========================================================================
%  NOTE: JACOBIAN CORRECTION
%  ========================================================================
%  Comparing a Log-Normal model fitlme(tbl, 'log(y) ~ ...') directly to a
%  Gamma or Inverse Gaussian model using raw Information Criteria (AIC) is
%  statistically invalid because the models operate on different scales.
%  The Log-Normal model calculates likelihoods based on the probability
%  density of the log-transformed data, which is a compressed scale with
%  naturally smaller variance. This artificial compression biases the AIC,
%  making the Log-Normal model appear to fit significantly better than it
%  actually does relative to models fitted on raw data. To rectify this, a
%  Jacobian correction is applied to the Log-Normal AIC. This mathematical
%  adjustment accounts for the derivative of the transformation,
%  effectively converting the likelihood of the log-values back into the
%  density of the original raw units. This ensures that all models in the
%  comparison are evaluated against the same underlying data variability,
%  allowing for a fair and rigorous selection between transforming the data
%  and using a generalized distribution.
%  ========================================================================


%% ========================================================================
%  NOTE: MODULATION INDEX (MI) VS. LOG-RATIO
%  ========================================================================
%  When quantifying firing rate changes (e.g., Ripple vs. Random periods), 
%  the choice of metric significantly impacts the statistical validity 
%  of the subsequent Linear Mixed-Effects (LME) analysis.
%
%  1. The Modulation Index (MI):
%     Defined as (FR_a - FR_b) / (FR_a + FR_b). This metric is a standard
%     "normalized difference" in neuroscience because it is intuitive and
%     controls for the absolute firing rate of a neuron. However, it
%     presents challenges for LME models. Because it is strictly bounded
%     between [-1, 1], it violates the assumption of normality, especially
%     near the boundaries. Analyzing MI typically requires rescaling the
%     data to (0, 1) and employing a Logit-Normal model.
%
%  2. The Log-Ratio (LR):
%     Defined as ln(FR_a / FR_b). This is often the statistically preferred
%     alternative to the MI. Unlike the MI, the Log-Ratio is not bounded
%     and is naturally symmetric around zero, which often results in
%     residuals that better approximate a Normal distribution. This
%     symmetry aligns with the "Log Link" philosophy, as it models
%     multiplicative processes on a linear scale.
%
%  3. The GLMM Approach (Raw Rates):
%     The most rigorous method is to avoid index-based metrics entirely. By
%     fitting a Generalized Linear Mixed Model (GLMM) directly to the raw
%     firing rates (e.g., using Gamma or Inverse Gaussian distributions),
%     one preserves the inherent mean-variance relationship of the data
%    . This approach accounts for the "multiplicative" nature 
%     of neural spiking while avoiding the information loss and potential 
%     biases (like the "retransformation problem") associated with 
%     collapsing observations into a single ratio.
%
%  DECISION RULE: Use the Modulation Index (MI) primarily for visualization
%  and consistency with literature. For LME modeling, favor the
%  **Log-Ratio** for simplicity of insference, or a **GLMM on raw rates**
%  (Gamma/InvGauss) to ensure the most accurate and unbiased population
%  estimates.
%  ========================================================================