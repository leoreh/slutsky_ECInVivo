function [barIdx, barLbl] = lme_sigLines(lmeStats, lmeMdl, lmeData, varargin)
% GET_SIGLINES Generates bar indices and labels for plot_sigLines from LME results.
%
% This function automatically determines the correct bar pairings and
% significance labels for plotting simple and marginal effects from a linear
% mixed-effects model analysis. It parses the effects results, maps
% them to the bar plot structure, and generates the inputs required by
% plot_sigLines.
%
% INPUTS:
%   lmeStats    (Required) Table: The output from lme_analyse.
%   lmeMdl      (Required) LME Model: The fitted model object from fitlme/fitglme.
%   lmeData     (Required) Table: The data used for the LME model.
%
% VARARGIN (Name-Value Pairs):
%   'idxRow'    (Optional) Vector of integers specifying which rows of
%               lmeStats to generate lines for. If empty (default), it
%               processes all "Simple", "Marginal", and relevant "Coeff"
%               effects.
%   'pThresh'   (Optional) Vector: p-value thresholds for significance stars.
%               Default is [0.05, 0.01, 0.001, 0.0001].
%   'sigLbls'   (Optional) Cell Array: Labels for significance levels.
%               Default is {'*', '**', '***', '****'}.
%   'grpVar'    (Optional) Char: Name of the variable to use for grouping.
%               Must match the grpVar used in lme_plot.
%
% OUTPUTS:
%   barIdx      Cell array of 2-element (or 1.5-style) vectors specifying
%               bar pairs or midpoints to connect with significance lines.
%   barLbl      Cell array of strings ('*', 'NS', etc.) corresponding to
%               the significance of each comparison in barIdx.
%
% See also: plot_sigLines, lme_analyse, lme_plot

% 260524 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addParameter(p, 'idxRow', [], @isnumeric);
addParameter(p, 'pThresh', [0.05, 0.01, 0.001, 0.0001], @isnumeric);
addParameter(p, 'sigLbls', {'*', '**', '***', '****'}, @iscell);
addParameter(p, 'grpVar', 'Group', @ischar);
parse(p, varargin{:});
opts = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
barIdx = {};
barLbl = {};

% Filter for effects based on idxRow or Type
if ~isempty(opts.idxRow)
    relStats = lmeStats(opts.idxRow, :);
else
    isRelevant = ismember(lmeStats.Type, ["Simple", "Marginal", "Coeff"]);
    if ~any(isRelevant)
        warning('No "Simple", "Marginal", or "Coeff" effects found in lmeStats. No significance lines generated.');
        return;
    end
    relStats = lmeStats(isRelevant, :);
end

% Determine which p-value column to use (adjusted or unadjusted)
if ismember('pAdj', relStats.Properties.VariableNames)
    pVals = relStats.pAdj;
    % For coefficients, pAdj is NaN, so fall back to pVal
    pVals(isnan(pVals)) = relStats.pVal(isnan(pVals));
else
    pVals = relStats.pVal;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET PLOT STRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frml = lmeMdl.Formula.char;
[varsFxd, varRsp] = lme_frml2vars(frml);

if length(varsFxd) < 2
    return;
end

[~, xVals, grpLbls] = lme_tbl2cell(lmeData, varsFxd, varRsp, 'grpVar', opts.grpVar);

xLevels = xVals{1};
grpLevels = grpLbls;
nXVals = length(xLevels);
nGrps = length(grpLevels);

% Map from (X_Level, Group_Level) to a unique bar index (1 to N)
% Assumes bars are clustered by x-level, then ordered by group-level.
barMap = containers.Map;
for iX = 1:nXVals
    for iGrp = 1:nGrps
        key = sprintf('%s__%s', xLevels{iX}, grpLevels{iGrp});
        idx = (iX - 1) * nGrps + iGrp;
        barMap(key) = idx;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARSE EFFECTS AND GENERATE OUTPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:height(relStats)
    desc = relStats.Description{i};
    p = pVals(i);
    
    idxPair = [];
    
    % Try parsing as a simple effect: (A vs B) at C
    tokensSimple = regexp(desc, '\(([^ ]+) vs ([^)]+)\) at (.+)', 'tokens');
    % Try parsing as a marginal effect: (A vs B) over C
    tokensMarginal = regexp(desc, '\(([^ ]+) vs ([^)]+)\) over (.+)', 'tokens');

    if ~isempty(tokensSimple)
        refLvl = tokensSimple{1}{1}; cmpLvl = tokensSimple{1}{2}; atLvl = tokensSimple{1}{3};
        
        if ismember(refLvl, xLevels) && ismember(cmpLvl, xLevels) && ismember(atLvl, grpLevels)
            key1 = sprintf('%s__%s', refLvl, atLvl); key2 = sprintf('%s__%s', cmpLvl, atLvl);
            if isKey(barMap, key1) && isKey(barMap, key2); idxPair = [barMap(key1), barMap(key2)]; end
        elseif ismember(refLvl, grpLevels) && ismember(cmpLvl, grpLevels) && ismember(atLvl, xLevels)
            key1 = sprintf('%s__%s', atLvl, refLvl); key2 = sprintf('%s__%s', atLvl, cmpLvl);
            if isKey(barMap, key1) && isKey(barMap, key2); idxPair = [barMap(key1), barMap(key2)]; end
        end

    elseif ~isempty(tokensMarginal)
        refLvl = tokensMarginal{1}{1}; cmpLvl = tokensMarginal{1}{2}; overVar = tokensMarginal{1}{3};

        % Find all bar indices for the reference and comparison levels
        refIdx = []; cmpIdx = [];
        if ismember(refLvl, xLevels) && ismember(cmpLvl, xLevels)
            for iGrp = 1:nGrps
                key1 = sprintf('%s__%s', refLvl, grpLevels{iGrp}); refIdx(end+1) = barMap(key1); %#ok<AGROW>
                key2 = sprintf('%s__%s', cmpLvl, grpLevels{iGrp}); cmpIdx(end+1) = barMap(key2); %#ok<AGROW>
            end
        elseif ismember(refLvl, grpLevels) && ismember(cmpLvl, grpLevels)
             for iX = 1:nXVals
                key1 = sprintf('%s__%s', xLevels{iX}, refLvl); refIdx(end+1) = barMap(key1); %#ok<AGROW>
                key2 = sprintf('%s__%s', xLevels{iX}, cmpLvl); cmpIdx(end+1) = barMap(key2); %#ok<AGROW>
            end
        end
        
        if ~isempty(refIdx) && ~isempty(cmpIdx)
            idxPair = [min(refIdx) + (max(refIdx)-min(refIdx))/2, ...
                       min(cmpIdx) + (max(cmpIdx)-min(cmpIdx))/2];
        end
    end
    
    if isempty(idxPair)
        continue;
    end
    
    % Determine significance label
    currentLbl = 'NS';
    for iThresh = length(opts.pThresh):-1:1
        if p < opts.pThresh(iThresh)
            currentLbl = opts.sigLbls{iThresh};
            break;
        end
    end

    barIdx{end+1} = sort(idxPair); %#ok<AGROW>
    barLbl{end+1} = currentLbl; %#ok<AGROW>
end

end 