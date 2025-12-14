function stats = lme_parkTest(tbl, frml, varargin)
% LME_PARKTEST Performs the Park Test to determine the Mean-Variance relationship.
%
%   STATS = LME_PARKTEST(TBL, FRML, VARARGIN) estimates a provisional
%   Generalized Linear Mixed Model (GLME) and performs the Modified Park Test
%   to determine the appropriate distribution family.
%
%   INPUTS:
%       tbl         - (table) Table containing the variables.
%       frml        - (char/formula) Model definition (e.g., 'y ~ x1 + (1|g)').
%                     NOTE: Response variable must be RAW data.
%       dist        - (char) Provisional distribution. {'Gamma'}
%       link        - (char) Provisional link function. {'Log'}
%       fitMethod   - (char) Fit method for fitglme. {'REMPL'}
%
%   OUTPUTS:
%       stats       - (struct) Structure containing Park Test results.
%

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addRequired(p, 'frml', @(x) ischar(x) || isa(x, 'classreg.regr.LinearFormula'));
addOptional(p, 'dist', 'Gamma', @ischar);
addOptional(p, 'link', 'Log', @ischar);
addOptional(p, 'fitMethod', 'REMPL', @ischar);

parse(p, tbl, frml, varargin{:});
dist = p.Results.dist;
link = p.Results.link;
fitMethod = p.Results.fitMethod;


%% ========================================================================
%  FIT PROVISIONAL MODEL (GLME)
%  ========================================================================

% Fit the Generalized Linear Mixed Model using the provisional settings
% default is usually Gamma/Log as a robust starting point.
glmeProv = fitglme(tbl, frml, 'Distribution', dist, 'Link', link, ...
    'FitMethod', fitMethod);


%% ========================================================================
%  CALCULATE RESIDUALS & PREDICTIONS
%  ========================================================================

respName = glmeProv.ResponseName;
yRaw = tbl.(respName);

% Conditional predictions (yHat)
% We use conditional predictions (including random effects) to isolate
% observation-level variance.
yHat = predict(glmeProv);

% Raw residuals
resRaw = yRaw - yHat;

% Squared residuals (Proxy for variance)
resSq = resRaw.^2;


%% ========================================================================
%  PARK REGRESSION
%  ========================================================================
%  Model: ln(residuals^2) = alpha + lambda * ln(y_hat)

% Filter for valid log inputs
isValid = (yHat > 0) & (resSq > 0);

if sum(isValid) < length(yRaw)
    warning('lme_parkTest:exclusions', ...
        '%d observations excluded due to non-positive predictions/residuals.', ...
        length(yRaw) - sum(isValid));
end

parkY = log(resSq(isValid));
parkX = log(yHat(isValid));

% Fit auxiliary OLS
lmPark = fitlm(parkX, parkY, 'VarNames', {'LogPred', 'LogResSq'});


%% ========================================================================
%  INTERPRET RESULTS
%  ========================================================================

lambdaEst = lmPark.Coefficients.Estimate(2);
lambdaSe  = lmPark.Coefficients.SE(2);
lambdaCi  = coefCI(lmPark);
lambdaCi  = lambdaCi(2, :);

recStr = 'Inconclusive';

% Heuristic check
is0 = (lambdaCi(1) <= 0) && (lambdaCi(2) >= 0);
is1 = (lambdaCi(1) <= 1) && (lambdaCi(2) >= 1);
is2 = (lambdaCi(1) <= 2) && (lambdaCi(2) >= 2);
is3 = (lambdaCi(1) <= 3) && (lambdaCi(2) >= 3);

if is2
    recStr = 'Gamma';
elseif is1
    recStr = 'Poisson';
elseif is0
    recStr = 'Gaussian';
elseif is3
    recStr = 'InverseGaussian';
end


%% ========================================================================
%  OUTPUT
%  ========================================================================

stats = struct();
stats.Lambda = lambdaEst;
stats.SE = lambdaSe;
stats.CI = lambdaCi;
stats.Recommendation = recStr;
stats.ProvisionalModel = glmeProv;
stats.ParkModel = lmPark;



end


%% ========================================================================
%  NOTE: THEORY
%  ========================================================================
%   The choice between Log-Linear (OLS on log(y)) and GLM (Gamma/Poisson)
%   depends on the structure of the variance. In GLMs, the variance is
%   assumed to be a function of the mean:
%
%       Var(y|x) = alpha * [E(y|x)]^lambda
%
%   The Park Test empirically estimates 'lambda' to identify the underlying
%   distribution family:
%
%       1. Fit a provisional model (usually Gamma GLM/GLME) to get predictions yHat.
%       2. Calculate raw residuals: r = (y - yHat).
%       3. Regress squared residuals on predictions in log-log space:
%          ln((y - yHat)^2) = alpha + lambda * ln(yHat) + error
%
%   INTERPRETATION OF LAMBDA:
%       lambda = 0 : Gaussian (Homoscedastic) -> Linear Model appropriate.
%       lambda = 1 : Poisson -> Poisson GLM appropriate.
%       lambda = 2 : Gamma   -> Gamma GLM appropriate.
%       lambda = 3 : Inverse Gaussian / Wald.
%
%   NOTE ON MIXED MODELS:
%   When dealing with LME/GLME, we use the conditional predictions (incorporating
%   random effects) to isolate the observation-level variance structure.
%
%   REFERENCES:
%   1. Manning, W. G., & Mullahy, J. (2001). Estimating log models: to transform
%      or not to transform? Journal of Health Economics, 20(4), 461-494.
%  ========================================================================