function contrastTbl = mcu_lmeContrasts(lme, varargin)
% analyzes specific contrasts from a fitted lme model
% 
% INPUT
%   lme             fitted lme model
%   
% OPTIONS
%   'contrasts'     cell array of specific contrasts {all}
%   'correction'    type of multiple comparison correction {'bonf', 'fdr'}
%
% OUTPUT
%   contrastTbl     table with results
%
% CALLS
%   mafdr
%
% 07 jan 24 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% argument validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addRequired(p, 'lme');
addParameter(p, 'contrasts', [], @(x) isempty(x) || iscell(x));
addParameter(p, 'correction', 'bonf', @(x) ismember(x, {'bonf', 'fdr'}));
parse(p, lme, varargin{:});
args = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get model terms and prepare contrasts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get coefficient names and remove intercept
coefNames = lme.CoefficientNames;
coefNames = coefNames(~strcmp(coefNames, '(Intercept)'));

% separate main effects and interactions
isInteraction = contains(coefNames, ':');
mainEffects = coefNames(~isInteraction);
interactions = coefNames(isInteraction);

% get factor names from main effects (remove _X from names)
factorNames = cellfun(@(x) strtok(x, '_'), mainEffects, 'uni', false);
uniqueFactors = unique(factorNames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build contrast matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize
C = [];
contrastNames = {};

% add main effect contrasts
for ifact = 1 : length(uniqueFactors)
    % find all coefficients for this factor
    factCoefs = contains(coefNames, [uniqueFactors{ifact} '_']);
    
    if sum(factCoefs) > 0
        % make contrast vector
        cVec = zeros(1, length(coefNames) + 1);  % +1 for intercept
        cVec(find(factCoefs) + 1) = 1;  % +1 because first col is intercept
        
        % add to matrix
        C = [C; cVec];
        contrastNames = [contrastNames; uniqueFactors{ifact}];
    end
end

% add interaction contrasts
for iint = 1 : length(interactions)
    % get factors involved in interaction
    intFactors = strsplit(interactions{iint}, ':');
    
    % make contrast vector for interaction
    cVec = zeros(1, length(coefNames) + 1);
    
    % add main effects
    for ifact = 1:length(intFactors)
        factCoefs = contains(coefNames, [intFactors{ifact} '_']);
        cVec(find(factCoefs) + 1) = 1;
    end
    
    % add interaction term
    intCoef = strcmp(coefNames, interactions{iint});
    cVec(find(intCoef) + 1) = 1;
    
    % add to matrix
    C = [C; cVec];
    contrastNames = [contrastNames; strjoin(intFactors, ':')];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select contrasts if specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(args.contrasts)
    validContrasts = ismember(contrastNames, args.contrasts);
    C = C(validContrasts, :);
    contrastNames = contrastNames(validContrasts);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test each contrast separately
nContrasts = size(C, 1);
p = zeros(nContrasts, 1);
f = zeros(nContrasts, 1);
df = zeros(nContrasts, 1);

for ic = 1:nContrasts
    [p(ic), f(ic), df(ic)] = coefTest(lme, C(ic,:), 0);
end

% multiple comparison correction
switch args.correction
    case 'bonf'
        p_corr = min(p * nContrasts, 1);
    case 'fdr'
        p_corr = mafdr(p, 'BHFDR', true);
end

% calculate estimates
estimates = C * [lme.Coefficients.Estimate];  % add 0 for intercept

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contrastTbl = table(contrastNames, estimates, f, df, p, p_corr,...
    'VariableNames', {'Contrast', 'Estimate', 'FStat', 'DF', 'pValue', 'pCorrected'});

end

% EOF

% In Linear Mixed Effects models, contrast testing can be performed either
% individually or jointly through coefTest. When testing contrasts
% individually by passing one row of the contrast matrix at a time, each
% test evaluates a specific comparison (e.g., Group 1 vs Control for RS
% units) and returns its own p-value. However, when passing multiple rows
% simultaneously to coefTest, it performs a joint hypothesis test examining
% whether ALL specified contrasts are simultaneously equal to zero.