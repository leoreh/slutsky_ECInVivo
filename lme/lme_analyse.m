function [lmeMdl, lmeStats, lmeInfo, lmeTbl] = lme_analyse(tbl, frml, varargin)
% LME_ANALYSE Wrapper for LME/GLME analysis pipeline.
%
%   [STATS, MDL, INFO] = LME_ANALYSE(TBL, FRML, ...) orchestrates the full
%   analysis workflow:
%
%   LME_VIF to examine collinearity between numeric predictors
%   LME_COMPAREDISTS to select the best distribution via AIC
%   LME_PARKTEST to allow for Gamma/Log-Normal checks.
%   TBL_TRANS:
%       - Predictors: Automatically Log-transforms skewed continuous
%         predictors and Z-scores all continuous predictors.
%       - Response: Applies Log or Logit transforms if required by dist.
%         Adds offset for Zero-inflated data if Gamma/Log is selected.
%   LME_FIT to fit the model.
%   LME_EFFECTS to generate ANOVA and contrasts.
%   LME_PLOTRES to plot residuals
%
%   INPUTS:
%       tbl         - (table) Data table.
%       frml        - (char/string) Model formula.
%       varargin    - (param/value) Optional parameters:
%                     'dist'        : (char) Force distribution. If empty, auto-selects.
%                     'contrasts', 'correction', 'dfMethod': See LME_EFFECTS.
%                     'flgPlot'
%
%   OUTPUTS:
%       lmeStats    - (table) Statistical results (ANOVA, Effects).
%       lmeMdl      - (object) Fitted model object.
%       lmeInfo     - (struct) Metadata (dist used, transforms applied).
%       lmeTbl      - (table) with data actually used for the model
%
%   See also: LME_FIT, LME_EFFECTS, LME_COMPAREDISTS, TBL_TRANS

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addRequired(p, 'frml', @(x) ischar(x) || isstring(x));

% Wrapper Options
addParameter(p, 'dist', '', @ischar);
addParameter(p, 'flgPlot', false, @islogical);
addParameter(p, 'fitMethod', '', @ischar);

% Pass-through Options (LME_EFFECTS)
addParameter(p, 'contrasts', 'all');
addParameter(p, 'correction', 'holm', @ischar);
addParameter(p, 'dfMethod', 'Satterthwaite', @ischar);
addParameter(p, 'verbose', true, @islogical);
parse(p, tbl, frml, varargin{:});

dist = p.Results.dist;
verbose = p.Results.verbose;
flgPlot = p.Results.flgPlot;
fitMethod = p.Results.fitMethod;
frml = char(frml);

% Run Info Storage
lmeInfo = struct();
lmeInfo.distInput = dist;


%% ========================================================================
%  PREDICTOR TRANSFORMATION
%  ========================================================================

% Vars
[varsFxd, varResp, varsRnd, intrTerms] = lme_frml2vars(frml);

% Extract actual grouping variables from random effect terms
varsGrp = {};
if ~isempty(varsRnd)
    for iRand = 1:length(varsRnd)
        tokens = regexp(varsRnd{iRand}, '\|([^)]+)\)', 'tokens');
        grpTerm = strtrim(tokens{1}{1});
        
        % Split by ':' or '*' to handle interactions (e.g., Name:UnitID)
        grpVars = strtrim(split(grpTerm, {':', '*'}));
        varsGrp = [varsGrp; grpVars(:)];
    end
end
varsGrp = unique(varsGrp);
varsFxd = unique(varsFxd);

% Extract variables from interaction terms (varsIntr)
% varsIntr contains strings like 'A:B'. We need 'A' and 'B'.
varsIntr = {};
if ~isempty(intrTerms)
    for iIntr = 1:length(intrTerms)
        atoms = strtrim(strsplit(intrTerms{iIntr}, ':'));
        varsIntr = [varsIntr, atoms];
    end
end
varsIntr = unique(varsIntr(:));

% Combine all required columns (Response + Fixed + Interaction + Random Groups)
varsTbl = unique([varsFxd(:); varsIntr(:); varsGrp(:); {varResp}]);

% Truncate Table
lmeTbl = tbl(:, varsTbl);

% Identify numeric/continuous predictors (exclude categorical)
isNum = cellfun(@(x) isnumeric(lmeTbl.(x)) && ~iscategorical(lmeTbl.(x)), varsTbl);
varsNum = varsTbl(isNum);
varsNum = setdiff(varsNum, varResp); % Exclude Response from Predictor Z-Scoring

if ~isempty(varsNum)
    if verbose
        fprintf('[LME_ANALYSE] Transforming Predictors\n');
    end
    [lmeTbl, transParams] = tbl_trans(lmeTbl, 'varsInc', varsNum, ...
        'logBase', 10, 'skewThr', 2, ...
        'flgZ', true, ...
        'verbose', verbose);

    lmeInfo.transParams = transParams;
else
    lmeInfo.transParams = struct();
    lmeInfo.transParams.varsTrans = struct();
    lmeInfo.transParams.varsGrp = {}; % Initialize if empty
end

%% ========================================================================
%  VIF (Collinearity)
%  ========================================================================
if numel(varsNum) >= 2
    if verbose; fprintf('[LME_ANALYSE] Checking Collinearity (VIF)...\n'); end
    lmeInfo.vif = lme_vif(lmeTbl, frml);
else
    lmeInfo.vif = table();
end


%% ========================================================================
%  DISTRIBUTION SELECTION
%  ========================================================================

wMsg = ['stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:',...
    'Message_PosteriorModeLineSearch'];
warning('off', wMsg);

if isempty(dist)
    if verbose
        fprintf('[LME_ANALYSE] Auto-selecting distribution...\n');
    end

    % Check for Binomial (Binary) Response
    uResp = unique(lmeTbl.(varResp));
    isBinomial = islogical(lmeTbl.(varResp)) || ...
        (length(uResp) == 2 && all(ismember(uResp, [0, 1])));

    if isBinomial
        if verbose
            fprintf('[LME_ANALYSE] Detected Binary Response -> Distribution: Binomial\n');
        end
        dist = 'Binomial';
    else

        parkStats = lme_parkTest(lmeTbl, frml, 'flgPlot', flgPlot);
        compStats = lme_compareDists(lmeTbl, frml);

        % Pick best model (lowest AIC)
        bestModel = compStats.Model{1};
        if verbose
            fprintf('[LME_ANALYSE] Selected: %s (AIC diff to 2nd: %.1f)\n', ...
                bestModel, compStats.AIC(2) - compStats.AIC(1));
        end

        dist = bestModel;
        lmeInfo.compStats = compStats;
        lmeInfo.parkStats = parkStats;
    end
end

lmeInfo.distSelected = dist;

warning('on', wMsg);


%% ========================================================================
%  RESPONSE TRANSFORMATION
%  ========================================================================

switch lower(dist)
    case 'log-normal'
        % Transform: Log(y), Dist: Normal
        if verbose
            fprintf('[LME_ANALYSE] Transforming Response: Log-Normal -> Log(Y)\n');
        end
        [lmeTbl, pResp] = tbl_trans(lmeTbl, 'varsInc', {varResp}, ...
            'logBase', 'e', ...
            'verbose', verbose);
        % Merge Response Params
        lmeInfo.transParams.varsTrans.(varResp) = pResp.varsTrans.(varResp);
        dist = 'Normal';

    case 'logit-normal'
        % Transform: Logit(y), Dist: Normal
        if verbose
            fprintf('[LME_ANALYSE] Transforming Response: Logit-Normal -> Logit(Y)\n');
        end
        [lmeTbl, pResp] = tbl_trans(lmeTbl, 'varsInc', {varResp}, ...
            'logBase', 'logit', ...
            'verbose', verbose);
        % Merge Response Params
        lmeInfo.transParams.varsTrans.(varResp) = pResp.varsTrans.(varResp);
        dist = 'Normal';

    case {'gamma', 'poisson', 'inversegaussian'}
        [lmeTbl, pResp] = tbl_trans(lmeTbl, 'varsInc', {varResp}, ...
            'flg0', true, 'verbose', verbose);
        if isfield(pResp, 'varsTrans') && isfield(pResp.varsTrans, varResp)
            lmeInfo.transParams.varsTrans.(varResp) = pResp.varsTrans.(varResp);
        end
end

lmeInfo.distFinal = dist;


%% ========================================================================
%  MODEL FITTING
%  ========================================================================

lmeMdl = lme_fit(lmeTbl, frml, 'dist', dist, 'fitMethod', fitMethod);

if flgPlot
    lme_plotRes(lmeMdl, 'Name', 'LME Analyse Diagnostics');
end


%% ========================================================================
%  EFFECTS ANALYSIS
%  ========================================================================

lmeStats = lme_effects(lmeMdl, ...
    'contrasts', p.Results.contrasts, ...
    'correction', p.Results.correction, ...
    'dfMethod', p.Results.dfMethod);


%% ========================================================================
%  FINALIZE INFO & DISPLAY
%  ========================================================================

% Add Metadata (AIC, Formula, FitMethod, DFMethod)
lmeInfo.frml = frml;
lmeInfo.dfMethod = p.Results.dfMethod;
lmeInfo.fitMethod = fitMethod;
lmeInfo.aic = lmeMdl.ModelCriterion.AIC;

% Verbose Output
if verbose
    fprintf('\n_______________________________________________________\n');
    fprintf('  LME INFO\n');
    disp(lmeInfo);
    fprintf('\n\n');
    fprintf('  LME STATS\n');
    disp(removevars(lmeStats, {'Index', 'Type', 'pAdj', 'HVec'}));
    fprintf('\n_______________________________________________________\n');
end


end     % EOF




