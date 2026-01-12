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
    warning(['Response variable "%s" contains negative values.', ...
        'Cannot fit Log-Normal, Gamma, or Poisson.'], varResp);
end

% Assert non-zeros
if nZeros > 0
    warning(['Response variable "%s" has %.1f%% zeros.', ...
        'Results may be unstable.'], varResp, pctZeros*100);

    tbl = tbl_trans(tbl, 'flg0', true, 'verbose', false, ...
        'varsInc', {varResp});
end

% Update yRaw after potential transformations (offsets)
yRaw = tbl.(varResp);



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

            % Log-Transform Response
            mdlTbl = tbl_trans(tbl, 'logBase', 'e', ...
                'verbose', false, 'skewThr', -Inf, 'varsInc', {varResp});
            dist = 'Normal';

            % Jacobian Correction: LL_raw = LL_ln - sum(ln(y))
            % Use original y values (stored in yRaw, or retrieve from tbl)
            jcb = -sum(log(yRaw), 'omitnan');

        elseif strcmp(dist, 'Logit-Normal')

            % Check Domain: Logit is only valid for values in [0, 1]
            if max(yRaw) > 1 || min(yRaw) < 0
                stats.ErrorMsg{iMdl} = sprintf('Data range [%.2f, %.2f] invalid for Logit [0, 1]', min(yRaw), max(yRaw));
                continue;
            end

            mdlTbl = tbl_trans(tbl, 'logBase', 'logit', 'verbose', false, ...
                'varsInc', {varResp});
            dist = 'Normal';

            % Jacobian Correction: LL_raw = LL_logit - sum(ln(y*(1-y)))
            % Clip for stability (consistent with tbl_trans)
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
%  NOTE: LOGIT-NORMAL ARTIFACT
%  ========================================================================
%  The 'Logit-Normal' model may occasionally appear as the "best" fit
%  (lowest AIC) for some data (eg, raw firing rate), even when the data is
%  not naturally bounded between 0 and 1.
%
%  1. Why it happens:
%     The Logit transformation, ln(y / (1 - y)), is strictly designed for
%     proportions or indices bounded by [0, 1] (e.g., Modulation Index).
%     If raw firing rates (which can exceed 1 Hz) are passed to this
%     transform, the results are mathematically invalid. The extremely low
%     AIC typically arises because the Jacobian correction used to equate
%     scales becomes numerically unstable or produces nonsensical
%     likelihoods when the data range violates the (0, 1) constraint.
%
%  2. The "Scaling" Trap:
%     If the firing rate data was normalized (e.g., divided by the maximum
%     rate) to fit the 0-1 range, the Logit-Normal model might "win" because
%     the logit transform aggressively expands the tails near 0 and 1.
%     However, this creates a "humped" distribution assumption that rarely
%     matches the monotonic decay (Exponential/Gamma shape) of typical
%     neural spiking data.
%
%  TLDR: If the data is not a proportion, the Logit-Normal AIC is a
%  mathematical ghost. Move to the next distribution in the list.
%  ========================================================================