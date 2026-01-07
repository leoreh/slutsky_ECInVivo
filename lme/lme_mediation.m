function res = lme_mediation(tbl, frml, xVar, mVar, varargin)
% LME_MEDIATION Performs a causal legacy Mediation Analysis using Mixed-Effects Models.
%
%   RES = LME_MEDIATION(TBL, FRML, XVAR, MVAR, ...) conducts a 4-step
%   mediation analysis.
%
%   STEPS (Baron & Kenny, 1986 adapted):
%       1. Path C  (Total):  Y ~ X + Cov + RE  (Input Formula)
%       2. Path A  (Med):    M ~ X + Cov + RE
%       3. Path B/C' (Dir):  Y ~ X + M + Cov + RE
%
%   INPUTS:
%       tbl         - (table) Raw Data table.
%       frml        - (char) Formula for Total Effect (Y ~ X + ...).
%       xVar        - (char) Independent Variable (Treatment).
%       mVar        - (char) Mediator Variable (Mechanism).
%       varargin    - (param/value) Optional parameters:
%                     'distM' : (char) Distribution for Mediator.
%                     'distY' : (char) Distribution for Outcome.
%                     'verbose': (logical, default true).
%
%   OUTPUTS:
%       res         - (struct) Results with models and path table.
%
%   See also: LME_ANALYSE, LME_FIT

%% ========================================================================
%  INPUT PARSING
%  ========================================================================
p = inputParser;
addRequired(p, 'tbl', @istable);
addRequired(p, 'frml', @ischar);
addRequired(p, 'xVar', @ischar);
addRequired(p, 'mVar', @ischar);
addParameter(p, 'distM', '', @ischar);
addParameter(p, 'distY', '', @ischar);
addParameter(p, 'verbose', true, @islogical);
parse(p, tbl, frml, xVar, mVar, varargin{:});

distM = p.Results.distM;
distY = p.Results.distY;
flgVerbose = p.Results.verbose;

% Parse Formula
[~, yVar] = lme_frml2vars(frml);
rhs = strtrim(extractAfter(frml, '~'));


%% ========================================================================
%  PATH C: X -> Y (TOTAL EFFECT)
%  ========================================================================
%  Does X predict Y?
if flgVerbose, fprintf('\n[LME_MEDIATION] Step 1: Total Effect (X->Y)\n'); end

% Use input formula directly
frmlC = frml;
[mdlC, ~, infoC] = lme_analyse(tbl, frmlC, 'dist', distY);
[betaC, pC, seC] = get_coeff(mdlC, xVar);


%% ========================================================================
%  PATH A: X -> M (MEDIATOR MODEL)
%  ========================================================================
%  Does X predict M? (Using same covariates/RE as Y model)
if flgVerbose, fprintf('\n[LME_MEDIATION] Step 2: Mediator Model (X->M)\n'); end

frmlA = sprintf('%s ~ %s', mVar, rhs);
[mdlA, ~, infoA] = lme_analyse(tbl, frmlA, 'dist', distM);
[betaA, pA, seA] = get_coeff(mdlA, xVar);


%% ========================================================================
%  PATH B & C': X + M -> Y (OUTCOME MODEL)
%  ========================================================================
%  Does M predict Y (Path B) controlling for X?
%  Does X still predict Y (Path C') controlling for M?
if flgVerbose, fprintf('\n[LME_MEDIATION] Step 3: Outcome Model (X+M->Y)\n'); end

% Add M to the predictors (RHS)
frmlBC = sprintf('%s ~ %s + %s', yVar, mVar, rhs);

% Reuse distY logic from Path C if it was auto-selected?
[mdlBC, ~, infoBC] = lme_analyse(tbl, frmlBC, 'dist', distY);

[betaB, pB, seB] = get_coeff(mdlBC, mVar);
[betaC_prime, pC_prime, seC_prime] = get_coeff(mdlBC, xVar);


%% ========================================================================
%  RESULTS & STATISTICS
%  ========================================================================

% Sobel Test for Indirect Effect (A * B)
% Z = (a*b) / sqrt(b^2*sa^2 + a^2*sb^2)
% Valid even if variables are scaled differently (Z-invariance).
indirectEffect = betaA * betaB;
seIndirect = sqrt(betaB^2 * seA^2 + betaA^2 * seB^2);
zSobel = indirectEffect / seIndirect;
pSobel = 2 * (1 - normcdf(abs(zSobel))); % Two-tailed

res.mdlA  = mdlA;
res.mdlC  = mdlC;
res.mdlBC = mdlBC;
res.infoA = infoA;
res.infoC = infoC;

% Summary Table
RowNames = {'Path A (X->M)'; 'Path B (M->Y|X)'; 'Path C (Total X->Y)'; 'Path C'' (Direct X->Y|M)'; 'Mediation (Sobel)'};
Variable  = {xVar; mVar; xVar; xVar; 'Indirect'};
Estimate  = [betaA; betaB; betaC; betaC_prime; indirectEffect];
SE        = [seA; seB; seC; seC_prime; seIndirect];
PValue    = [pA; pB; pC; pC_prime; pSobel];

res.paths = table(Variable, Estimate, SE, PValue, 'RowNames', RowNames);


%% ========================================================================
%  DISPLAY
%  ========================================================================
if flgVerbose
    fprintf('\n=======================================================\n');
    fprintf(' MEDIATION ANALYSIS: %s -> %s -> %s\n', xVar, mVar, yVar);
    fprintf('=======================================================\n');
    fprintf('%-25s | Est=%8.3f | SE=%7.3f | p=%8.4f\n', 'Path A (X->M)', betaA, seA, pA);
    fprintf('%-25s | Est=%8.3f | SE=%7.3f | p=%8.4f\n', 'Path B (M->Y)', betaB, seB, pB);
    fprintf('%-25s | Est=%8.3f | SE=%7.3f | p=%8.4f\n', 'Path C (Total)', betaC, seC, pC);
    fprintf('%-25s | Est=%8.3f | SE=%7.3f | p=%8.4f\n', 'Path C'' (Direct)', betaC_prime, seC_prime, pC_prime);
    fprintf('-------------------------------------------------------\n');
    fprintf('%-25s | Est=%8.3f | Z =%7.3f | p=%8.4f\n', 'Sobel (Indirect)', indirectEffect, zSobel, pSobel);

    % Interpretation
    if pSobel < 0.05
        if pC_prime > 0.05
            fprintf('RESULT: Full Mediation (Significant Indirect, Non-sig Direct)\n');
        else
            fprintf('RESULT: Partial Mediation (Significant Indirect & Direct)\n');
        end
    else
        fprintf('RESULT: No Significant Mediation (pSobel > 0.05)\n');
    end
    fprintf('=======================================================\n');
end

end


%% ========================================================================
%  HELPER: GET COEFFICIENT
%  ========================================================================
function [est, pval, se] = get_coeff(mdl, varName)

allNames = mdl.Coefficients.Name;
idx = find(strcmp(allNames, varName));

if isempty(idx)
    % Try partial match (e.g. 'Group_TG')
    idx = find(contains(allNames, varName) & ~contains(allNames, ':'));
end

if isempty(idx)
    est = NaN; pval = NaN; se = NaN;
    warning('Variable %s not found in coefficients.', varName);
else
    % Take first match (Reference)
    est  = mdl.Coefficients.Estimate(idx(1));
    pval = mdl.Coefficients.pValue(idx(1));
    se   = mdl.Coefficients.SE(idx(1));
end

end


%% ========================================================================
%  NOTE: MEDIATION ANALYSIS 
%  ========================================================================
% Mediation analysis is a statistical method used to elucidate the mechanism
% or "pathway" through which an independent variable (X) influences a
% dependent variable (Y). It posits that X influences a third variable, the
% mediator (M), which in turn influences Y.
%
% 1. The Four Conditions (Baron & Kenny, 1986):
%    To establish mediation, four conditions typically need to be met:
%
%    Path A (X -> M): Use LME/GLME to show X significantly predicts M.
%        Formula: M ~ X + (1|Covariates)
%        Interpretation: The treatment must affect the proposed mechanism.
%
%    Path C (Total Effect, X -> Y): Use LME/GLME to show X predicts Y.
%        Formula: Y ~ X + (1|Covariates)
%        Interpretation: There is an effect to be mediated.
%
%    Path B (M -> Y | X): Use LME/GLME to show M predicts Y when controlling for X.
%        Formula: Y ~ X + M + (1|Covariates)
%        Interpretation: The mechanism affects the outcome independent of the treatment.
%
%    Path C' (Direct Effect, X -> Y | M): In the same model as Path B, the
%    effect of X on Y should act as follows:
%        - Full Mediation: Path C' is no longer significant.
%        - Partial Mediation: Path C' is smaller than Path C but still significant.
%
% 2. Mixed-Effects Context:
%    Standard mediation relies on OLS regression (General Linear Model).
%    However, in physiological experiments with hierarchical data (cells
%    nested within animals), we MUST use Mixed-Effects Models. Ignoring
%    clustering typically leads to Type I errors (false positives) for Path A
%    and Path C. This function wraps `fitglme` to perform these steps
%    correctly while respecting the random effects structure.
%
% 3. Causality Warning:
%    Mediation is a statistical test of correlations, not a proof of
%    causality. Even if all paths are significant, M could be a correlate
%    of the true cause, or Y could cause M (reverse causality). Strong
%    causal claims require experimental manipulation of the mediator (e.g.,
%    blocking bFrac directly) rather than just statistical adjustment.
%
% ========================================================================

%% ========================================================================
%  NOTE: INTERPRETATION OF COEFFICIENTS
%  ========================================================================
% - Total Effect (Path C): The overall impact of X on Y.
% - Direct Effect (Path C'): The impact of X on Y that is NOT seemingly
%   due to M.
% - Indirect Effect (A * B): The portion of the effect passing through M.
%
% Significance Testing:
% The Sobel test is a common method to test the significance of the
% indirect efffect (A*B). However, it assumes normal sampling distributions
% which often doesn't hold for the product of coefficients. Bootstrapping
% is the modern gold standard but is computationally expensive for GLMEs.
% This function relies on the joint significance logic of Paths A and B.
%
% ========================================================================

%% ========================================================================
%  NOTE: COMPETITIVE MEDIATION (SUPPRESSION)
%  ========================================================================
% The analysis may reveal a phenomenon known as Competitive Mediation (or
% Suppression) where the Direct Effect (C') is larger in magnitude than the
% Total Effect (C).
%
% 1. Mechanism:
%    This occurs when the two pathways work in opposite directions:
%    - Direct Path: The Treatment (X) has a negative impact on the Outcome (Y).
%    - Indirect Path: The Treatment (X) increases the Mediator (M), and the
%      Mediator (M) has a positive impact on the Outcome (Y).
%
% 2. Interpretation:
%    In this scenario, the Mediator acts as a "suppressor" variable. It
%    "hides" a portion of the Treatment's negative effect by providing a
%    compensatory boost. When you control for the Mediator in the model
%    (Path C'), the "pure" negative impact of the Treatment becomes more
%    pronounced (larger beta) because the masking effect is removed.
%
% 3. Impact on Feature Importance:
%    Competitive mediation can make Importance Analysis (Ablation) confusing.
%    - Information Overlap: The model sees X and M pulling predictions in
%      opposite directions.
%    - Accuracy Sensitivity: Removing the "negative" predictor (X) might be
%      partially offset by the "positive" predictor (M) remaining, leading
%      to a deceptively small drop in accuracy.
%    - Conclusion: In these cases, Likelihood Ratio Tests (LRT) or p-values
%      are often more reliable metrics of feature relevance than simple
%      prediction accuracy.
% ========================================================================