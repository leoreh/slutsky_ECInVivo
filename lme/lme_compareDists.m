function tblStats = lme_compareDists(tbl, frml)
% LME_COMPAREDISTS Compares GLME fits across multiple distributions.
%
%   TBLSTATS = LME_COMPAREDISTS(TBL, FRML) fits Generalized Linear
%   Mixed Models using four common distributions (Normal, Poisson, Gamma,
%   InverseGaussian) and returns a table of fit statistics.
%
%   INPUTS:
%       tbl         - (table) Table containing the variables.
%       frml        - (char/formula) Model definition.
%       fitMethod   - (char) Fit method for fitglme. {'REMPL'}
%
%   OUTPUTS:
%       tblStats    - (table) Summary of model performance (AIC, BIC, etc.).
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
dists = {'Normal', 'Log-Normal', 'Poisson', 'Gamma', 'InverseGaussian'};
nMdl = length(dists);

% Common canonical/standard links
% Normal -> Identity
% Log-Normal -> Identity (on log data)
% Others -> Log (Enforced as per best practice for positive biological data)
links = {'Identity', 'Identity', 'Log', 'Log', 'Log'};

% Output table
tblStats = table();
tblStats.Distribution = dists(:);
tblStats.Link = links(:);
tblStats.AIC = nan(nMdl, 1);
tblStats.BIC = nan(nMdl, 1);
tblStats.LogLikelihood = nan(nMdl, 1);
tblStats.R2_Ordinary = nan(nMdl, 1);
tblStats.R2_Adjusted = nan(nMdl, 1);
tblStats.Converged = false(nMdl, 1);
tblStats.ErrorMsg = repmat({''}, nMdl, 1);


%% ========================================================================
%  RUN FITS
%  ========================================================================

% Parse response variable for Log-Normal / Jacobian Check
respVar = strtrim(strsplit(frml, '~'));
respVar = respVar{1};
yRaw = tbl.(respVar);

% Hardcoded for Model Comparison
fitMethod = 'Laplace';

for iMdl = 1:nMdl
    curDist = tblStats.Distribution{iMdl};
    curLink = tblStats.Link{iMdl};

    try
        % -----------------------------------------------------------------
        % LOG-NORMAL (Special Case)
        % -----------------------------------------------------------------
        if strcmp(curDist, 'Log-Normal')

            % Fit Normal GLME on log(y)
            logFrml = ['log(' respVar ') ~ ' extractAfter(frml, '~')];

            mdl = fitglme(tbl, logFrml, ...
                'Distribution', 'Normal', ...
                'Link', 'Identity', ...
                'FitMethod', fitMethod);

            % Extract Raw LogLikelihood (of log-data)
            llLog = mdl.LogLikelihood;

            % Jacobian Correction
            % LL_raw = LL_log - sum(log(y))
            % This converts the density from log-scale back to raw-scale
            sumLogY = nansum(log(yRaw));
            llRaw = llLog - sumLogY;

            % Recalculate ICs on raw scale
            % Back-calculate k (num params) from AIC and LL
            % AIC = -2*LL + 2*k  =>  2*k = AIC + 2*LL
            aicLog = mdl.ModelCriterion.AIC;
            k2 = aicLog + 2*llLog; % 2*k

            nObs = mdl.NumObservations;

            aicRaw = -2*llRaw + k2;
            bicRaw = -2*llRaw + (k2/2)*log(nObs);

            % Store
            tblStats.AIC(iMdl) = aicRaw;
            tblStats.BIC(iMdl) = bicRaw;
            tblStats.LogLikelihood(iMdl) = llRaw;
            tblStats.R2_Ordinary(iMdl) = mdl.Rsquared.Ordinary;
            tblStats.R2_Adjusted(iMdl) = mdl.Rsquared.Adjusted;
            tblStats.Converged(iMdl) = true;

            % -----------------------------------------------------------------
            % STANDARD GLME (Normal, Poisson, Gamma, IG)
            % -----------------------------------------------------------------
        else
            mdl = fitglme(tbl, frml, ...
                'Distribution', curDist, ...
                'Link', curLink, ...
                'FitMethod', fitMethod);

            tblStats.AIC(iMdl) = mdl.ModelCriterion.AIC;
            tblStats.BIC(iMdl) = mdl.ModelCriterion.BIC;
            tblStats.LogLikelihood(iMdl) = mdl.LogLikelihood;
            tblStats.R2_Ordinary(iMdl) = mdl.Rsquared.Ordinary;
            tblStats.R2_Adjusted(iMdl) = mdl.Rsquared.Adjusted;
            tblStats.Converged(iMdl) = true;
        end

    catch ME
        % Handle Errors
        tblStats.ErrorMsg{iMdl} = ME.message;
    end
end

% Sort by AIC (Low to High) - Best model first
[~, sortIdx] = sort(tblStats.AIC);
tblStats = tblStats(sortIdx, :);

end


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