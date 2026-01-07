function [lmeMdl, lmeStats, lmeInfo, lmeTbl] = lme_analyse(tbl, frml, varargin)
% LME_ANALYSE Wrapper for LME/GLME analysis pipeline.
%
%   [STATS, MDL, INFO] = LME_ANALYSE(TBL, FRML, ...) orchestrates the full
%   analysis workflow:
%
%   Distribution Selection: If 'dist' is not provided, runs
%   LME_COMPAREDISTS to select the best distribution via AIC and
%   LME_PARKTEST to allow for Gamma/Log-Normal checks.
%   Transformation:
%       - Predictors: Automatically Log-transforms skewed continuous
%         predictors and Z-scores all continuous predictors.
%       - Response: Applies Log or Logit transforms if required by dist.
%         Adds offset for Zero-inflated data if Gamma/Log is selected.
%   Fitting: Calls LME_FIT to fit the model.
%   Analysis: Calls LME_EFFECTS to generate ANOVA and contrasts.
%
%   INPUTS:
%       tbl         - (table) Data table.
%       frml        - (char/string) Model formula.
%       varargin    - (param/value) Optional parameters:
%                     'dist'        : (char) Force distribution. If empty, auto-selects.
%                     'contrasts', 'correction', 'dfMethod': See LME_EFFECTS.
%
%   OUTPUTS:
%       lmeStats    - (table) Statistical results (ANOVA, Effects).
%       lmeMdl      - (object) Fitted model object.
%       lmeInfo     - (struct) Metadata (dist used, transforms applied).
%       lmeTbl      - (table) with data actually used for the model
%
%   See also: LME_FIT, LME_EFFECTS, LME_COMPAREDISTS, TBL_TRANSFORM

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addRequired(p, 'frml', @(x) ischar(x) || isstring(x));

% Wrapper Options
addParameter(p, 'dist', '', @ischar);

% Pass-through Options (LME_EFFECTS)
addParameter(p, 'contrasts', 'all');
addParameter(p, 'correction', 'holm', @ischar);
addParameter(p, 'dfMethod', 'Satterthwaite', @ischar);
parse(p, tbl, frml, varargin{:});

dist = p.Results.dist;
frml = char(frml);

% Run Info Storage
lmeInfo = struct();
lmeInfo.distInput = dist;


%% ========================================================================
%  PRE-PROCESSING & TRANSFORMATION
%  ========================================================================

% Vars
[varsFxd, varResp, varsRnd] = lme_frml2vars(frml);

% Extract actual grouping variables from random effect terms (e.g., '(1|Name)' -> 'Name')
varsGrp = {};
if ~isempty(varsRnd)
    for iRand = 1:length(varsRnd)
        tokens = regexp(varsRnd{iRand}, '\|([^)]+)\)', 'tokens');
        varsGrp = [varsGrp; strtrim(tokens{1}{1})];
    end
end
varsGrp = unique(varsGrp);

% Truncate Table
vars = unique([{varResp}, varsFxd, varsGrp]);
lmeTbl = tbl(:, vars);

% PREDICTOR TRANSFORMATION
% Identify numeric/continuous predictors (exclude categorical)
isNum = cellfun(@(x) isnumeric(lmeTbl.(x)) && ~iscategorical(lmeTbl.(x)), varsFxd);
varsCont = varsFxd(isNum);

if ~isempty(varsCont)
    fprintf('[LME_ANALYSE] Transforming Predictors (Log10[Skew>2] + Z-Score)\n');
    lmeTbl = tbl_transform(lmeTbl, 'varsInc', varsCont, ...
        'logBase', 10, 'skewThr', 2, ...
        'flgZ', true, ...
        'verbose', true);
end


%% ========================================================================
%  DISTRIBUTION SELECTION
%  ========================================================================

wMsg = ['stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:',...
    'Message_PosteriorModeLineSearch'];
warning('off', wMsg);

if isempty(dist)
    fprintf('[LME_ANALYSE] Auto-selecting distribution...\n');

    % Check for Binomial (Binary) Response
    uResp = unique(lmeTbl.(varResp));
    isBinomial = islogical(lmeTbl.(varResp)) || ...
        (length(uResp) == 2 && all(ismember(uResp, [0, 1])));

    if isBinomial
        fprintf('[LME_ANALYSE] Detected Binary Response -> Distribution: Binomial\n');
        dist = 'Binomial';
    else

        parkStats = lme_parkTest(lmeTbl, frml);
        compStats = lme_compareDists(lmeTbl, frml);

        % Pick best model (lowest AIC)
        bestModel = compStats.Model{1};
        fprintf('[LME_ANALYSE] Selected: %s (AIC diff to 2nd: %.1f)\n', ...
            bestModel, compStats.AIC(2) - compStats.AIC(1));

        dist = bestModel;
        lmeInfo.compStats = compStats;
        lmeInfo.parkStats = parkStats;
    end
end

lmeInfo.distSelected = dist;

warning('on', wMsg);


%% ========================================================================
%  FINAL RESPONSE TRANSFORMATION
%  ========================================================================

switch lower(dist)
    case 'log-normal'
        % Transform: Log(y), Dist: Normal
        fprintf('[LME_ANALYSE] Transforming Response: Log-Normal -> Log(Y)\n');
        lmeTbl = tbl_transform(lmeTbl, 'varsInc', {varResp}, ...
            'logBase', 'e', ...
            'verbose', true);
        dist = 'Normal';

    case 'logit-normal'
        % Transform: Logit(y), Dist: Normal
        fprintf('[LME_ANALYSE] Transforming Response: Logit-Normal -> Logit(Y)\n');
        lmeTbl = tbl_transform(lmeTbl, 'varsInc', {varResp}, ...
            'logBase', 'logit', ...
            'verbose', true);
        dist = 'Normal';

    case {'gamma', 'poisson', 'inversegaussian'}
        fprintf('[LME_ANALYSE] Zero-inflation detected for %s. Adding offset.\n', dist);
        lmeTbl = tbl_transform(lmeTbl, 'varsInc', {varResp}, ...
            'flg0', true, 'verbose', true);
end

lmeInfo.distFinal = dist;


%% ========================================================================
%  MODEL FITTING
%  ========================================================================

lmeMdl = lme_fit(lmeTbl, frml, 'dist', dist);


%% ========================================================================
%  EFFECTS ANALYSIS
%  ========================================================================

lmeStats = lme_effects(lmeMdl, ...
    'contrasts', p.Results.contrasts, ...
    'correction', p.Results.correction, ...
    'dfMethod', p.Results.dfMethod);


end     % EOF




