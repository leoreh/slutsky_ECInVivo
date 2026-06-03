function tbl = lme_simpleAt(mdl, factorVar, covVar, atVals, varargin)
% LME_SIMPLEAT Simple effect of a categorical factor at set covariate values.
%
%   TBL = LME_SIMPLEAT(MDL, FACTORVAR, COVVAR, ATVALS, ...) computes, for a
%   fitted (G)LME containing a FACTORVAR*COVVAR interaction, the contrast of
%   each non-reference level of FACTORVAR against its reference, evaluated at
%   each value in ATVALS of the continuous predictor COVVAR. This is a
%   "spotlight" / simple-slopes analysis: it reads the interaction at chosen
%   covariate values rather than only at the (reference-coded) zero point.
%
%   The contrast vector L is built from the factor main coefficient plus
%   COVVAL * interaction coefficient, and statistics are computed from the
%   model coefficient covariance via COEFTEST, identically to LME_POSTHOC.
%
%   INPUTS:
%       mdl         - (object) Fitted LinearMixedModel / GeneralizedLinearMixedModel.
%       factorVar   - (char) Categorical predictor (reference = first category).
%       covVar      - (char) Continuous predictor interacting with factorVar.
%       atVals      - (numeric) Covariate values to evaluate at, in MODEL space
%                     (i.e. the same units the model was fit on; if the
%                     predictor was transformed, pass transformed values such
%                     as group means of mdl.Variables.(covVar)).
%       ...         - (param/value):
%                     'atLabels'    : (cellstr) one label per ATVALS for the
%                                     Description column {numeric value}.
%                     'transParams' : (struct) from LME_ANALYSE; if the response
%                                     was log-transformed, adds Fold + FoldCI95
%                                     columns on the response scale.
%                     'dfMethod'    : (char) 'Satterthwaite' (default) / 'Residual'.
%
%   OUTPUT:
%       tbl         - (table) One row per (non-ref level x atVals) with columns
%                     Description, Estimate, CI95, SE, Statistic, DF, P-value
%                     (and Fold, FoldCI95 when a log response transform is found).
%
%   See also: LME_POSTHOC, LME_LSMEANS, COEFTEST

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'mdl', @(x) isa(x,'LinearMixedModel') || isa(x,'GeneralizedLinearMixedModel'));
addRequired(p, 'factorVar', @(x) ischar(x) || isstring(x));
addRequired(p, 'covVar', @(x) ischar(x) || isstring(x));
addRequired(p, 'atVals', @isnumeric);
addParameter(p, 'atLabels', {}, @(x) iscell(x) || isempty(x));
addParameter(p, 'transParams', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'dfMethod', 'Satterthwaite', @ischar);
parse(p, mdl, factorVar, covVar, atVals, varargin{:});

factorVar = char(factorVar);
covVar    = char(covVar);
atVals    = atVals(:);
atLabels  = p.Results.atLabels;
transP    = p.Results.transParams;
dfMethod  = p.Results.dfMethod;

if isa(mdl, 'GeneralizedLinearMixedModel'), dfMethod = 'Residual'; end
if isempty(atLabels)
    atLabels = arrayfun(@(v) sprintf('%s=%.3g', covVar, v), atVals, 'Uni', false);
end

%% ========================================================================
%  MODEL DETAILS
%  ========================================================================

cn    = mdl.CoefficientNames;
b     = mdl.Coefficients.Estimate;
V     = mdl.CoefficientCovariance;

% Factor levels (reference = first category)
lvls   = categories(mdl.Variables.(factorVar));
refLvl = lvls{1};
nonRef = lvls(2:end);

% Detect a natural-log response transform for back-transformation
isLogResp = false;
if ~isempty(transP) && isfield(transP, 'varsTrans') && ...
        isfield(transP.varsTrans, mdl.ResponseName)
    lb = transP.varsTrans.(mdl.ResponseName).logBase;
    isLogResp = (ischar(lb) && strcmpi(lb,'e'));
end

%% ========================================================================
%  BUILD & TEST CONTRASTS
%  ========================================================================

rows = {};
for iLvl = 1:numel(nonRef)
    lvl     = nonRef{iLvl};
    iMain   = find(strcmp(cn, sprintf('%s_%s', factorVar, lvl)));
    iInt    = find_interaction(cn, covVar, sprintf('%s_%s', factorVar, lvl));
    if isempty(iMain) || isempty(iInt)
        error('lme_simpleAt:CoefNotFound', ...
            'Could not locate main/interaction coefficients for %s_%s x %s.', ...
            factorVar, lvl, covVar);
    end

    for iAt = 1:numel(atVals)
        L        = zeros(1, numel(cn));
        L(iMain) = 1;
        L(iInt)  = atVals(iAt);

        est = L * b;
        se  = sqrt(max(0, L * V * L'));
        [pVal, ~, ~, df] = coefTest(mdl, L, 0, 'DFMethod', dfMethod);
        tCrit = tinv(0.975, df);
        ci    = [est - tCrit*se, est + tCrit*se];

        desc = sprintf('(%s vs %s) at %s', refLvl, lvl, atLabels{iAt});
        row  = {desc, round(est,4), sprintf('[%.4g %.4g]', ci(1), ci(2)), ...
                round(se,4), round(est/se,2), round(df,1), round(pVal,4)};

        if isLogResp
            row = [row, {round(exp(est),3), ...
                sprintf('[%.3g %.3g]', exp(ci(1)), exp(ci(2)))}]; %#ok<AGROW>
        end
        rows(end+1, :) = row; %#ok<AGROW>
    end
end

varNames = {'Description','Estimate','CI95','SE','Statistic','DF','P-value'};
if isLogResp, varNames = [varNames, {'Fold','FoldCI95'}]; end
tbl = cell2table(rows, 'VariableNames', varNames);

end     % EOF

%% ========================================================================
%  HELPER
%  ========================================================================

function idx = find_interaction(cn, covVar, factorLvlStr)
% Locate the interaction coefficient covVar:factorLvlStr regardless of the
% order MATLAB used when naming it.
cand = {sprintf('%s:%s', covVar, factorLvlStr), ...
        sprintf('%s:%s', factorLvlStr, covVar)};
idx = find(ismember(cn, cand), 1);
end
