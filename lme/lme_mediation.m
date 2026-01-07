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