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
%                     'distM'   : (char) Distribution for Mediator.
%                     'distY'   : (char) Distribution for Outcome.
%                     'transTemplate' : (struct) Explicit template to standardize parameters for the predictors {[]}.
%                     'verbose' : (logical, default true).
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
addParameter(p, 'transTemplate', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'verbose', true, @islogical);
parse(p, tbl, frml, xVar, mVar, varargin{:});

distM = p.Results.distM;
distY = p.Results.distY;
transTemplate = p.Results.transTemplate;
flgVerbose = p.Results.verbose;

%% ========================================================================
%  PREP
%  =======================================================================

% Determine fitMethod for Gamma distributions
fitMethodM = ''; if strcmpi(distM, 'gamma'), fitMethodM = 'Laplace'; end
fitMethodY = ''; if strcmpi(distY, 'gamma'), fitMethodY = 'Laplace'; end

% Parse Formula
[varsFxd, yVar, varsRand] = lme_frml2vars(frml);
rhs = strtrim(extractAfter(frml, '~'));

varsGrp = {};
for iVar = 1:numel(varsRand)
    tok = regexp(varsRand{iVar}, '\|([^)]+)\)', 'tokens', 'once');
    if ~isempty(tok)
        varsGrp = [varsGrp, strtrim(strsplit(tok{1}, {':', '*'}))]; %#ok<AGROW>
    end
end
vars = unique([{yVar, xVar, mVar}, varsFxd, varsGrp], 'stable');

% Keep only necessary columns 
tbl = tbl(:, vars);

% Remove rows with missing values (NaN/Missing)
idxMissing = any(ismissing(tbl), 2);
if any(idxMissing)
    nBef = size(tbl, 1);
    tbl = tbl(~idxMissing, :);
    nRem = nBef - size(tbl, 1);
    warning('LME_MEDIATION: Detected %d missing rows', nRem);
end

%  =======================================================================
%  PATH C: X -> Y (TOTAL EFFECT)
%  ========================================================================
%  Does X predict Y?
if flgVerbose, fprintf('\n[LME_MEDIATION] Step 1: Total Effect (X->Y)\n'); end

% Use input formula directly
frmlC = frml;
[mdlC, statsC, infoC] = lme_analyse(tbl, frmlC, 'dist', distY, 'fitMethod', fitMethodY, 'transTemplate', transTemplate, 'verbose', false);
[betaC, pC, seC, statC, dfC, ciC] = get_coeff(mdlC, xVar);


%% ========================================================================
%  PATH A: X -> M (MEDIATOR MODEL)
%  ========================================================================
%  Does X predict M? (Using same covariates/RE as Y model)
if flgVerbose, fprintf('\n[LME_MEDIATION] Step 2: Mediator Model (X->M)\n'); end

frmlA = sprintf('%s ~ %s', mVar, rhs);
[mdlA, statsA, infoA] = lme_analyse(tbl, frmlA, 'dist', distM, 'fitMethod', fitMethodM, 'transTemplate', transTemplate, 'verbose', true);
[betaA, pA, seA, statA, dfA, ciA] = get_coeff(mdlA, xVar);


%% ========================================================================
%  PATH B & C': X + M -> Y (OUTCOME MODEL)
%  ========================================================================
%  PATH B:  Does M predict Y controlling for X?
%  PATH C': Does X still predict Y controlling for M?
if flgVerbose, fprintf('\n[LME_MEDIATION] Step 3: Outcome Model (X+M->Y)\n'); end

% --- LINK CONSISTENCY ---
% If Path A used a non-identity Link (e.g. Log for Gamma/Poisson), 
% then betaA is in units of Link(M).
% For the mediation product (A*B) to be valid, M must enter Path B 
% as a predictor on that same Link scale.

% If Log-Link was used but not captured as a template (Gamma/Poisson case), 
% inject it into the template so tbl_trans can apply it to Path B.
isLogM = isprop(mdlA, 'Link') && strcmpi(mdlA.Link.Name, 'Log');
if isLogM
    infoA.transParams.varsTrans.(mVar).logBase = 'e';
end

% Construct explicit template for Path BC
tmplBC = transTemplate;
if isempty(tmplBC)
    tmplBC = struct('varsTrans', struct(), 'varsGrp', [], 'varNorm', '', 'catRef', []);
end
% Inject M's explicit transformation block from Path A into the Path BC template
% NOTE: lme_analyse will apply this transformation to M as a predictor in Step 3.
if isfield(infoA, 'transParams') && isfield(infoA.transParams, 'varsTrans') && isfield(infoA.transParams.varsTrans, mVar)
    tmplBC.varsTrans.(mVar) = infoA.transParams.varsTrans.(mVar);
    % tmplBC.varsTrans.(mVar).flgZ = true; % Force standardization for predictor role
end

% Add M to the predictors (RHS)
frmlBC = sprintf('%s ~ %s + %s', yVar, mVar, rhs);

% Reuse distY logic from Path C
[mdlBC, statsBC, infoBC, tblBC] = lme_analyse(tbl, frmlBC, 'dist', distY, 'fitMethod', fitMethodY, 'transTemplate', tmplBC, 'verbose', true);

[betaB, pB, seB, statB, dfB, ciB] = get_coeff(mdlBC, mVar);
[betaC_prime, pC_prime, seC_prime, statC_prime, dfC_prime, ciC_prime] = get_coeff(mdlBC, xVar);


%% ========================================================================
%  RESULTS & STATISTICS
%  ========================================================================

% UNIT CONSISTENCY CHECK:
% Path A (X->M): betaA = d(TransM) / d(Z_X)
% Path B (M->Y): betaB = d(Y)      / d(TransM)
%
% If Path BC was standardized, betaB is d(Y)/d(Z_TransM).
% We must un-Z-score it to get $d(Y)/d(TransM)$.
flgZ_M = false;
if isfield(infoBC.transParams.varsTrans, mVar)
    flgZ_M = infoBC.transParams.varsTrans.(mVar).flgZ;
end

if flgZ_M
    % We retrieve the standard deviation on the transformed scale from Path A.
    if isfield(infoA.transParams.varsTrans, mVar)
        sdM = infoA.transParams.varsTrans.(mVar).stats.SD(1);
    else
        sdM = 1;
    end
    betaB = betaB / sdM;
    seB   = seB   / sdM;
    ciB   = {[ciB{1}(1) / sdM, ciB{1}(2) / sdM]};
end

% Sobel Test for Indirect Effect (A * B)
% Z = (a*b) / sqrt(b^2*sa^2 + a^2*sb^2)
indirectEffect = betaA * betaB;
seIndirect = sqrt(betaB^2 * seA^2 + betaA^2 * seB^2);
zSobel = indirectEffect / seIndirect;
pSobel = 2 * (1 - normcdf(abs(zSobel))); % Two-tailed
statSobel = zSobel;
dfSobel = NaN;
ciSobel = {[indirectEffect - 1.96*seIndirect, indirectEffect + 1.96*seIndirect]};

res.mdlA  = mdlA;
res.mdlC  = mdlC;
res.mdlBC = mdlBC;
res.statsA  = statsA;
res.statsC  = statsC;
res.statsBC = statsBC;
res.infoA = infoA;
res.infoC = infoC;
res.infoBC = infoBC;

% Summary Table
Description = string({'Path A (X->M)'; 'Path B (M->Y|X)'; 'Path C (Total X->Y)'; 'Path C'' (Direct X->Y|M)'; 'Sobel test (Indirect X->M->Y)'});
Estimate    = round([betaA; betaB; betaC; betaC_prime; indirectEffect], 3);
CI95        = {mat2str(round(ciA{1}, 2)); mat2str(round(ciB{1}, 2)); mat2str(round(ciC{1}, 2)); mat2str(round(ciC_prime{1}, 2)); mat2str(round(ciSobel{1}, 2))};
SE          = round([seA; seB; seC; seC_prime; seIndirect], 3);
tStatistic  = round([statA; statB; statC; statC_prime; statSobel], 2);
DF          = [dfA; dfB; dfC; dfC_prime; dfSobel];
PValue      = round([pA; pB; pC; pC_prime; pSobel], 4);

res.paths = table(Description, Estimate, string(CI95), SE, tStatistic, DF, PValue, ...
    'VariableNames', {'Description', 'Estimate', 'CI95', 'SE', 't-statistic', 'DF', 'P-value'});

res.xlsTbls = mediation2xls(res);

%% ========================================================================
%  PLOTTING DATA
%  ========================================================================

res.data = tbl; 

% --- CALCULATE PARTIAL DATA FOR PLOTTING ---
% Strategy: Use simple vector algebra to adjust Y.
%
% 1. Path B Update (M -> Y | X):
%    We want to visualize the effect of M on Y, controlling for X.
%    Partial Residual = Residuals + Beta_B * M
%    This is equivalent to Y_adjusted = Beta_B * M + Epsilon
%    (Component + Residual plot)
%
% 2. Path C' Update (X -> Y | M):
%    We want to visualize the effect of X on Y, controlling for M.
%    Y_adjusted = Y - Beta_B * M
%    (Removes the effect of M from Y, leaving X + Epsilon)

rawResid = residuals(mdlBC, 'ResidualType', 'Raw');

% Store Data for Plotting (Table Format)
X = tbl.(xVar);          % X is from input table (original)
M = tbl.(mVar);          % M is from input table (transformed if applicable)
Y = tblBC.(yVar);        % Y from fitted table (transformed if log-normal)

fitA = fitted(mdlA);     % Path A Fit (M ~ X)
fitC = fitted(mdlC);     % Path C Fit (Y ~ X)

% Partial Calculations
% Note: betaB is already un-standardized (in units of Y / M_transformed).
Y_part_M = rawResid + (betaB * M);             % Path B Partial (Effect of M)
Y_part_X = Y - (betaB * M);                    % Path C' Partial (Effect of X)

res.plot = table(X, M, Y, fitA, fitC, Y_part_M, Y_part_X);



%% ========================================================================
%  NOTE: PARTIAL RESIDUALS & LME_PR
%  ========================================================================
%  The Partial Residuals calculated here for Path B (Y_part_M) are conceptually
%  identical to those produced by LME_PR in 'residual' mode (Component + Residual).
%
%  Formula:
%    Y_part_M = Residuals(Y|X,M) + Beta_B * M
%
%  DIFFERENCE IN SCALING (X-AXIS):
%  - Here (LME_MEDIATION): The X-axis (M) is plotted in its Raw (or Log-Transformed)
%    units. This preserves the physical interpretation of the Mediator's scale.
%  - LME_PR: By default, LME_PR plots against the predictor values stored in
%    the model object. LME_ANALYSE automatically Z-scores continuous predictors.
%    Therefore, LME_PR's X-axis will be in Standard Deviations (Z-scores).
%
%  To replicate this plot exactly using LME_PR, you must un-standardize the X-axis
%  by providing the 'transParams' structure.
%  ========================================================================

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
    fprintf('%-30s | Est=%8.3f | Z =%7.3f | p=%8.4f\n', 'Sobel test (Indirect X->M->Y)', indirectEffect, zSobel, pSobel);

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
function [est, pval, se, stat, df, ci95] = get_coeff(mdl, varName)

allNames = mdl.Coefficients.Name;
idx = find(strcmp(allNames, varName));

if isempty(idx)
    % Try partial match (e.g. 'Group_TG')
    idx = find(contains(allNames, varName) & ~contains(allNames, ':'));
end

if isempty(idx)
    est = NaN; pval = NaN; se = NaN; stat = NaN; df = NaN; ci95 = {[NaN, NaN]};
    warning('Variable %s not found in coefficients.', varName);
else
    % Take first match (Reference)
    est  = mdl.Coefficients.Estimate(idx(1));
    pval = mdl.Coefficients.pValue(idx(1));
    se   = mdl.Coefficients.SE(idx(1));
    stat = mdl.Coefficients.tStat(idx(1));
    df   = mdl.Coefficients.DF(idx(1));
    ci95 = {[mdl.Coefficients.Lower(idx(1)), mdl.Coefficients.Upper(idx(1))]};
end

end


%% ========================================================================
%  HELPER: MEDIATION TO XLS
%  ========================================================================

function medTbls = mediation2xls(res)
% MEDIATION2XLS Formats mediation results into a structure array for lme_save.

medTbls = struct('Title', {}, 'Table', {});

% Main summary
medTbls(end+1).Title = 'MEDIATION SUMMARY';
medTbls(end).Table = res.paths;

% PATH A
medTbls(end+1).Title = 'PATH A (X->M)';
medTbls(end).Table = table();
tblA = lme_mdl2tbls(res.mdlA, res.statsA, res.infoA);
for i = 1:length(tblA)
    medTbls(end+1).Title = tblA(i).Title;
    medTbls(end).Table = tblA(i).Table;
end

% PATH B & C'
medTbls(end+1).Title = 'PATH B & C'' (X+M->Y)';
medTbls(end).Table = table();
tblBC = lme_mdl2tbls(res.mdlBC, res.statsBC, res.infoBC);
for i = 1:length(tblBC)
    medTbls(end+1).Title = tblBC(i).Title;
    medTbls(end).Table = tblBC(i).Table;
end

% PATH C
medTbls(end+1).Title = 'PATH C (Total X->Y)';
medTbls(end).Table = table();
tblC = lme_mdl2tbls(res.mdlC, res.statsC, res.infoC);
for i = 1:length(tblC)
    medTbls(end+1).Title = tblC(i).Title;
    medTbls(end).Table = tblC(i).Table;
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
% ========================================================================