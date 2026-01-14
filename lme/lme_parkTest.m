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
%   Gamma) and handles zero-inflation via `tbl_trans` if required.
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
    fitMet = 'Laplace';

    % Handle Zeros (if < 5%)
    if nZeros > 0
        warning('Response "%s" has %.1f%% zeros. Applying offset.', ...
            varResp, pctZeros*100);
        tbl = tbl_trans(tbl, 'flg0', true, 'varsInc', {varResp}, ...
            'flgZ', false, 'verbose', false);

        % Update yRaw for consistency
        yRaw = tbl.(varResp);
    end

else
    % --- Normal Candidate ---
    dist = 'Normal';
    fitMet = 'ML';
end

%  PROVISIONAL FIT
mdlProv = lme_fit(tbl, frml, 'dist', dist, 'fitMethod', fitMet);


%% ========================================================================
%  PARK TEST (VARIANCE REGRESSION)
%  ========================================================================
%  Model: ln(residuals^2) = alpha + lambda * ln(y_hat)
%  We use a Gamma GLM with Log link on squared residuals to avoid retransformation bias.

% Get Predictions (Conditional on Random Effects)
yHat = predict(mdlProv);

% Calculate Squared Residuals
resRaw = yRaw - yHat;
resSq = resRaw.^2;

% Filter for valid log inputs
isValid = (yHat > 0) & (resSq > 0);

if sum(isValid) < length(yRaw) * 0.9
    warning('%d%% of observations excluded from Park Test (non-positive preds/resid).', ...
        round((1 - mean(isValid))*100));
end

% Auxiliary Regression
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
    % Fit Log-Normal Model (LME on log(y))
    % Use tbl_trans to ensure log applied correctly
    % Note: We use 'e' base for mathematical consistency with Log-Normal definition
    tblLog = tbl_trans(tbl, 'logBase', 'e', 'varsInc', {varResp}, ...
        'flgZ', false, 'verbose', false);

    mdlLog = lme_fit(tblLog, frml, 'dist', 'Normal');

    % Check Residual Kurtosis
    resLog = residuals(mdlLog);
    kurtLog = kurtosis(resLog);
    stats.KurtosisLog = kurtLog;

    % Check DF (Small Sample Size)
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

end     % EOF

