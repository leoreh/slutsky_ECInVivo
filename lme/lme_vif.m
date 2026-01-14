function vif = lme_vif(tbl, frml)
% LME_VIF Calculates Variance Inflation Factors (VIF) for numeric predictors.
%
%   VIF = LME_VIF(TBL, FRML) extracts fixed effects from the formula
%   and calculates VIF for the numeric ones.
%
%   INPUTS:
%       tbl         - (table) Data table containing predictors.
%       frml        - (char/string) LME formula string.
%
%   OUTPUTS:
%       vif         - (table) Table with columns:
%                       .Predictor (cell) Variable Name
%                       .VIF       (double) Variance Inflation Factor
%
%   NOTES:
%       - VIF = 1 / (1 - R^2_i), where R^2_i is the R-squared of the
%         regression of predictor Xi against all other predictors.
%       - Equivalent to diag(inv(corrcoef(X))).
%       - Only applicable to purely numeric predictors. Categorical variables
%         are ignored.
%
%   See also: LME_ANALYSE, CORRCOEF

%% ========================================================================
%  INPUT PARSING
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addRequired(p, 'frml', @(x) ischar(x) || isstring(x));
parse(p, tbl, frml);

frml = p.Results.frml;

%% ========================================================================
%  VARIABLE EXTRACTION
%  ========================================================================

% Extract Fixed Effects from Formula
[varsFxd, ~, ~, ~] = lme_frml2vars(char(frml));
varsTrgt = varsFxd;

% Filter for existence in table
inTable = ismember(varsTrgt, tbl.Properties.VariableNames);
if ~all(inTable)
    varsMiss = varsTrgt(~inTable);
    warning('Variables not found in table: %s', strjoin(varsMiss, ', '));
    varsTrgt = varsTrgt(inTable);
end

if isempty(varsTrgt)
    vif = table({}, [], 'VariableNames', {'Predictor', 'VIF'});
    return;
end

%% ========================================================================
%  NUMERIC FILTERING
%  ========================================================================

% Filter numeric
isNum = varfun(@isnumeric, tbl(:, varsTrgt), 'OutputFormat', 'uniform');
varsTrgt = varsTrgt(isNum);

% Need at least 2 variables for VIF
if length(varsTrgt) < 2
    vif = table({}, [], 'VariableNames', {'Predictor', 'VIF'});
    return;
end

%% ========================================================================
%  CALCULATION
%  ========================================================================

% Extract Matrix
X = table2array(tbl(:, varsTrgt));

% Correlation Matrix
% Using 'Rows','pairwise' in corrcoef handles NaNs robustly.
R = corrcoef(X, 'Rows', 'pairwise');

% Invert Correlation Matrix
Rinv = inv(R);
vifVal = diag(Rinv);


%% ========================================================================
%  OUTPUT PACKAGING
%  ========================================================================

% Check for high VIF
if any(vifVal > 5)
    warning('Some predictors have VIF > 5\n');
end

vif = table(varsTrgt(:), vifVal(:), 'VariableNames', {'Predictor', 'VIF'});

end     % EOF




%% ========================================================================
%  NOTE: COLLINEARITY
%  ========================================================================
% Collinearity, or multicollinearity, occurs when two or more independent
% predictors in your model are highly correlated with one another. In
% physiological datasets, this is a common challenge because many
% metrics—such as mean firing rate, burst fraction, and network
% dimensionality—often capture overlapping aspects of the same underlying
% neural state. Technically, collinearity means that one predictor can be
% linearly predicted from the others with a high degree of accuracy, which
% violates the assumption that each predictor provides unique information
% to the model.
%
% While collinearity does not necessarily reduce the overall predictive
% power of the model, it severely compromises the interpretability of
% individual coefficients. In a mixed-effects framework, high collinearity
% inflates the standard errors of the fixed effects, making it difficult to
% achieve statistical significance even for variables with strong
% biological effects. Furthermore, the estimated coefficients become highly
% unstable; small changes in your dataset can cause the sign of a
% coefficient to flip from positive to negative, leading to erroneous
% physiological conclusions.
%
% To assess whether you need to intervene, you should calculate the
% Variance Inflation Factor (VIF) for each continuous predictor. The VIF
% quantifies how much the variance of an estimated coefficient is increased
% due to collinearity. A general rule of thumb in the literature is that a
% VIF exceeding 5 to 10 indicates problematic redundancy. Additionally,
% generating a pairwise correlation matrix of your predictors before
% running `fitglme` provides a direct visual assessment of which variables
% are "tracking" each other too closely.
%
%  Impact:
%     - Unstable coefficient estimates (Beta).
%     - High sensitivity to outliers.
%     - Masked significance of predictors.
%
% See:
% https://www.mathworks.com/help/econ/time-series-regression-ii-collinearity-and-estimator-variance.html
%
%  If VIF analysis reveals severe redundancy between variables (e.g., fr
%  and dimensionality), consider the following technical steps:
%
%  1. Variable Removal:
%     Drop the less biologically relevant predictor.
%
%  2. Dimensionality Reduction (PCA):
%     Combine correlated metrics into one.
%     Use the first Principal Component (PC1).
%
%  3. Feature Centering:
%     Standardizing continuous variables reduces collinearity.
%     Particularly helpful for interaction terms.
%  ========================================================================

%% ========================================================================
%  NOTE: RESIDUALIZATION (FEATURE DECOUPLING)
%  ========================================================================
% Residualization is a statistical technique used to isolate the unique
% variance of a predictor by "stripping away" its association with another
% variable. In neural data, it is frequently proposed when two metrics -
% such as firing rate (FR) and burst fraction (pBspk)— - are suspected of
% being mathematically or biologically coupled. The goal is to create a
% "purified" version of a metric that represents only the variance that
% cannot be explained by a primary driver.
%
% 1. Technical Implementation:
%    To residualize Burstiness (B) against Firing Rate (FR):
%    - Step A: Fit a linear regression model: B ~ FR.
%    - Step B: Extract the residuals (observed B - predicted B).
%    - Step C: Use these residuals as the new predictor in your GLME.
%    The resulting "Residualized Burstiness" is, by definition, perfectly
%    uncorrelated () with Firing Rate.
%
% 2. Theoretical Utility:
%    Residualization is technically valid in two specific scenarios:
%    - Severe Collinearity: When VIF > 5-10, making it impossible for the
%      model to distinguish between two predictors.
%    - Hierarchical Hypothesis Testing: When you want to force the model to
%      assign all shared variance to a "nuisance" variable first, testing
%      if the secondary variable has any remaining incremental value.
%
% 3. Why Residualization is Discouraged in This Context:
%    - Loss of Physiological Meaning: A residualized metric is an abstract
%      statistical construct. "Burstiness beyond what is expected for a
%      given FR" is difficult to interpret biologically and cannot be
%      compared across different experimental datasets.
%    - Low Redundancy (VIF < 1.5): If the raw correlation between FR and
%      pBspk is already low and the VIF is near 1, residualization is
%      mathematically redundant. The model is already capable of
%      partitioning the variance accurately.
%    - Altering the Null Hypothesis: Residualization changes the coefficient
%      of the primary variable (FR). It effectively "over-adjusts" the model,
%      potentially leading to Type I errors if the shared variance was
%      biologically meaningful.
%    - Masking the Interaction: If your hypothesis is that bursts are
%      "ineffective" in a specific Group (MCU-KO), residualization can
%      suppress the very interaction you seek to find. The "ineffectiveness"
%      is best captured by a `pBspk * Group` interaction term, which
%      directly models how the relationship between burstiness and
%      recovery changes across conditions.
%% ========================================================================