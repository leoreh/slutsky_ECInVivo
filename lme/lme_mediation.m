function res = lme_mediation(mdlA, mdlC, mdlBC, xVar, mVar, varargin)
% LME_MEDIATION Sobel test and summary for a pre-fitted mediation triplet.
%
%   RES = LME_MEDIATION(MDLA, MDLC, MDLBC, XVAR, MVAR, ...) computes the
%   Sobel test for the indirect effect (A*B) from three pre-fitted
%   mixed-effects models representing the Baron & Kenny (1986) mediation
%   framework:
%
%       Model A  (X->M):     M ~ X + Cov + RE
%       Model C  (Total):    Y ~ X + Cov + RE
%       Model BC (Direct):   Y ~ X + M + Cov + RE
%
%   The function does NOT fit any models. Models should be fitted externally
%   using LME_ANALYSE.
%
%   When 'grpVar' is specified and an interaction between a continuous
%   predictor (xVar or mVar) and the categorical grpVar is present,
%   group-specific coefficients are extracted via contrast vectors using the
%   same coefTest pattern as LME_POSTHOC. A separate Sobel test is computed
%   for each level of grpVar. This only handles interactions between a
%   continuous predictor and a single categorical variable; multi-way
%   categorical interactions are not supported.
%
%   LIMITATION: Predictors should NOT be Z-scored (flgStnd = false in
%   lme_analyse). If M is Z-scored in Model BC, the units of betaB become
%   Y/Z_M rather than Y/M, making the A*B product inconsistent. When using
%   unstandardized predictors (the typical use case), this is a non-issue.
%
%   INPUTS:
%       mdlA        - (LME/GLME) Pre-fitted Model A.
%       mdlC        - (LME/GLME) Pre-fitted Model C.
%       mdlBC       - (LME/GLME) Pre-fitted Model BC.
%       xVar        - (char) Independent variable name.
%       mVar        - (char) Mediator variable name.
%       varargin    - (param/value) Optional parameters:
%                     'grpVar' : (char) Categorical variable in the model
%                                interaction. When provided, a separate
%                                Sobel test is computed for each level {''}.
%                     'verbose': (logical) Print summary {true}.
%
%   OUTPUTS:
%       res         - (struct array) One element per group level (or scalar
%                     if no grpVar). Each element contains:
%                     .grpLevel : (char) Group level name ('' if ungrouped).
%                     .paths    : (table) 5-row summary (A, B, C, C', Sobel).
%
%   See also: LME_ANALYSE, LME_POSTHOC, LME_MEDIATIONPLOT

%% ========================================================================
%  INPUT PARSING
%  ========================================================================

flgMdl = @(x) isa(x, 'LinearMixedModel') || isa(x, 'GeneralizedLinearMixedModel');
p = inputParser;
addRequired(p, 'mdlA', flgMdl);
addRequired(p, 'mdlC', flgMdl);
addRequired(p, 'mdlBC', flgMdl);
addRequired(p, 'xVar', @ischar);
addRequired(p, 'mVar', @ischar);
addParameter(p, 'grpVar', '', @ischar);
addParameter(p, 'verbose', true, @islogical);
parse(p, mdlA, mdlC, mdlBC, xVar, mVar, varargin{:});

grpVar     = p.Results.grpVar;
flgVerbose = p.Results.verbose;


%% ========================================================================
%  DETERMINE GROUP LEVELS
%  ========================================================================

if isempty(grpVar)
    grpLevels = {''};
else
    grpLevels = categories(mdlA.Variables.(grpVar));
end
nGrp = numel(grpLevels);


%% ========================================================================
%  COMPUTE SOBEL TEST PER GROUP
%  ========================================================================

for iGrp = nGrp : -1 : 1

    grpLevel = grpLevels{iGrp};

    % --- Extract path coefficients ---
    if isempty(grpVar)
        [betaA, pA, seA, statA, dfA, ciA]      = get_coeff(mdlA,  xVar);
        [betaC, pC, seC, statC, dfC, ciC]      = get_coeff(mdlC,  xVar);
        [betaB, pB, seB, statB, dfB, ciB]      = get_coeff(mdlBC, mVar);
        [betaCp, pCp, seCp, statCp, dfCp, ciCp] = get_coeff(mdlBC, xVar);
    else
        [betaA, pA, seA, statA, dfA, ciA]      = get_coeff_grp(mdlA,  xVar, grpVar, grpLevel);
        [betaC, pC, seC, statC, dfC, ciC]      = get_coeff_grp(mdlC,  xVar, grpVar, grpLevel);
        [betaB, pB, seB, statB, dfB, ciB]      = get_coeff_grp(mdlBC, mVar, grpVar, grpLevel);
        [betaCp, pCp, seCp, statCp, dfCp, ciCp] = get_coeff_grp(mdlBC, xVar, grpVar, grpLevel);
    end

    % --- Sobel test ---
    % Z = (a*b) / sqrt(b^2*sa^2 + a^2*sb^2)
    ab     = betaA * betaB;
    seAB   = sqrt(betaB^2 * seA^2 + betaA^2 * seB^2);
    zSobel = ab / seAB;
    pSobel = 2 * (1 - normcdf(abs(zSobel)));
    ciAB   = [ab - 1.96 * seAB, ab + 1.96 * seAB];

    % --- Results table ---
    res(iGrp).grpLevel = grpLevel;
    res(iGrp).paths = table(...
        string({'Path A (X->M)'; 'Path B (M->Y|X)'; 'Path C (Total X->Y)'; ...
                'Path C'' (Direct X->Y|M)'; 'Sobel (Indirect X->M->Y)'}), ...
        round([betaA; betaB; betaC; betaCp; ab], 3), ...
        string({mat2str(round(ciA, 2)); mat2str(round(ciB, 2)); ...
                mat2str(round(ciC, 2)); mat2str(round(ciCp, 2)); ...
                mat2str(round(ciAB, 2))}), ...
        round([seA; seB; seC; seCp; seAB], 3), ...
        round([statA; statB; statC; statCp; zSobel], 2), ...
        [dfA; dfB; dfC; dfCp; NaN], ...
        round([pA; pB; pC; pCp; pSobel], 4), ...
        'VariableNames', {'Description', 'Estimate', 'CI95', 'SE', ...
                          'Statistic', 'DF', 'pValue'});

    % --- Display ---
    if flgVerbose
        yVar  = mdlBC.ResponseName;
        mName = mdlA.ResponseName;

        grpStr = '';
        if ~isempty(grpVar)
            grpStr = sprintf(' [%s = %s]', grpVar, grpLevel);
        end

        fprintf('\n=======================================================\n');
        fprintf(' MEDIATION: %s -> %s -> %s%s\n', xVar, mName, yVar, grpStr);
        fprintf('=======================================================\n');
        fprintf('%-25s | Est=%8.3f | SE=%7.3f | p=%8.4f\n', 'Path A (X->M)',    betaA, seA, pA);
        fprintf('%-25s | Est=%8.3f | SE=%7.3f | p=%8.4f\n', 'Path B (M->Y)',    betaB, seB, pB);
        fprintf('%-25s | Est=%8.3f | SE=%7.3f | p=%8.4f\n', 'Path C (Total)',   betaC, seC, pC);
        fprintf('%-25s | Est=%8.3f | SE=%7.3f | p=%8.4f\n', 'Path C'' (Direct)', betaCp, seCp, pCp);
        fprintf('-------------------------------------------------------\n');
        fprintf('%-25s | Est=%8.3f | Z =%7.3f | p=%8.4f\n', 'Sobel (Indirect)', ab, zSobel, pSobel);

        if pSobel < 0.05
            if pCp > 0.05
                fprintf('  => Full Mediation\n');
            else
                fprintf('  => Partial Mediation\n');
            end
        else
            fprintf('  => No Significant Mediation\n');
        end
        fprintf('=======================================================\n');
    end
end

end


%% ========================================================================
%  HELPER: GET COEFFICIENT (MAIN EFFECT)
%  ========================================================================

function [est, pval, se, stat, df, ci] = get_coeff(mdl, varName)
% GET_COEFF Extracts a coefficient from a fitted model by exact name match.

coefTbl  = mdl.Coefficients;
allNames = coefTbl.Name;
idx = find(strcmp(allNames, varName), 1);

if isempty(idx)
    % Partial match excluding interactions (e.g., 'Group_TG')
    idx = find(contains(allNames, varName) & ~contains(allNames, ':'), 1);
end

if isempty(idx)
    est = NaN; pval = NaN; se = NaN; stat = NaN; df = NaN; ci = [NaN, NaN];
    warning('LME_MEDIATION:CoefNotFound', ...
        'Variable "%s" not found in model coefficients.', varName);
    return;
end

est  = coefTbl.Estimate(idx);
pval = coefTbl.pValue(idx);
se   = coefTbl.SE(idx);
stat = coefTbl.tStat(idx);
df   = coefTbl.DF(idx);
ci   = [coefTbl.Lower(idx), coefTbl.Upper(idx)];

end


%% ========================================================================
%  HELPER: GET COEFFICIENT (GROUP-SPECIFIC VIA H-VECTOR)
%  ========================================================================

function [est, pval, se, stat, df, ci] = get_coeff_grp(mdl, varName, grpVar, grpLevel)
% GET_COEFF_GRP Group-specific coefficient from a continuous * categorical
% interaction model. For the reference level, returns the main effect. For
% non-reference levels, constructs an H-vector summing the main effect and
% interaction term, then uses coefTest (same pattern as LME_POSTHOC).

% --- Determine reference level ---
cats     = categories(mdl.Variables.(grpVar));
refLevel = cats{1};

if strcmp(grpLevel, refLevel)
    [est, pval, se, stat, df, ci] = get_coeff(mdl, varName);
    return;
end

% --- Build H-vector for non-reference level ---
coefNames = mdl.CoefficientNames;
nCoefs    = numel(coefNames);
coefMap   = containers.Map(coefNames, 1:nCoefs);
coefEst   = mdl.Coefficients.Estimate;
coefCov   = mdl.CoefficientCovariance;

hVec = zeros(1, nCoefs);

% Main effect index
idxMain = find(strcmp(coefNames, varName), 1);
if isempty(idxMain)
    idxMain = find(contains(coefNames, varName) & ~contains(coefNames, ':'), 1);
end
if isempty(idxMain)
    est = NaN; pval = NaN; se = NaN; stat = NaN; df = NaN; ci = [NaN, NaN];
    warning('LME_MEDIATION:CoefNotFound', ...
        'Variable "%s" not found in model coefficients.', varName);
    return;
end
hVec(idxMain) = 1;

% Interaction index (try both orderings via find_coef_name)
grpCoefStr = sprintf('%s_%s', grpVar, grpLevel);
mainName   = coefNames{idxMain};
[~, found, idxIntr] = find_coef_name({mainName, grpCoefStr}, coefMap);

if ~found
    warning('LME_MEDIATION:IntrNotFound', ...
        'Interaction %s:%s not found. Returning main effect only.', varName, grpCoefStr);
    [est, pval, se, stat, df, ci] = get_coeff(mdl, varName);
    return;
end
hVec(idxIntr) = 1;

% --- Compute statistics via coefTest ---
est = hVec * coefEst(:);
varH = hVec * coefCov * hVec';
if varH < 0 && abs(varH) < 1e-10
    varH = 0;
end
se = sqrt(varH);

% DF method: LME supports Satterthwaite, GLME does not
if isa(mdl, 'LinearMixedModel')
    dfMethod = 'Satterthwaite';
else
    dfMethod = 'Residual';
end

[pval, FVal, ~, df] = coefTest(mdl, hVec, 0, 'DFMethod', dfMethod);
stat  = sqrt(FVal) * sign(est);
tCrit = tinv(0.975, df);
ci    = [est - tCrit * se, est + tCrit * se];

end


%% ========================================================================
%  HELPER: FIND COEFFICIENT NAME (PERMUTATION MATCHING)
%  ========================================================================

function [coefName, found, idx] = find_coef_name(partStrs, coefMap)
% FIND_COEF_NAME Finds a coefficient in coefMap matching component strings
% joined by ':', regardless of ordering. Same pattern as LME_POSTHOC.

n    = numel(partStrs);
perm = perms(1:n);
for iP = 1:size(perm, 1)
    candidate = strjoin(partStrs(perm(iP, :)), ':');
    if isKey(coefMap, candidate)
        coefName = candidate;
        found    = true;
        idx      = coefMap(candidate);
        return;
    end
end
coefName = '';
found    = false;
idx      = [];

end


%% ========================================================================
%  NOTE: MEDIATION ANALYSIS
%  ========================================================================
% Mediation analysis is a statistical method to elucidate the mechanism
% through which an independent variable (X) influences a dependent variable
% (Y). It posits that X influences a mediator (M), which in turn
% influences Y.
%
% Steps (Baron & Kenny, 1986):
%   Path A  (X->M):   X must significantly predict M.
%   Path C  (Total):  X must significantly predict Y.
%   Path B  (M->Y|X): M must predict Y controlling for X.
%   Path C' (Direct): Effect of X on Y controlling for M.
%     - Full Mediation:    C' non-significant.
%     - Partial Mediation: C' significant but smaller than C.
%
% The Sobel test assesses the significance of the indirect effect (A*B)
% using a first-order delta method approximation for the SE. Bootstrapping
% is the modern gold standard but computationally expensive for GLMEs.
%
% Causality Warning: Mediation is a statistical test of correlations, not
% a proof of causality. Strong causal claims require experimental
% manipulation of the mediator.
% ========================================================================
