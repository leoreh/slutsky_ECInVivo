function lmeStats = lme_effects(mdl, varargin)
% LME_EFFECTS Comprehensive post-hoc analysis for LME/GLME models.
%
%   STATS = LME_EFFECTS(MDL, ...) analyzes a fitted Mixed-Effects Model
%   (LME or GLME). It generates an ANOVA table, extracts coefficients, and
%   computes Simple and Marginal effects for two-way interactions.
%   Multiple comparison correction is applied to derived effects.
%
%   INPUTS:
%       mdl         - (object) Fitted LinearMixedModel or GeneralizedLinearMixedModel.
%       varargin    - (param/value) Optional parameters:
%                     'contrasts' : (var) Filter for output rows.
%                                   'all' (default), numeric indices, or logical vector.
%                     'correction': (char) Multiple comparison correction method.
%                                   'holm' (default), 'bonferroni', 'fdr', 'none'.
%                     'dfMethod'  : (char) DF method for t-tests.
%                                   'Satterthwaite' (default for LME), 'Residual'.
%
%   OUTPUTS:
%       lmeStats    - (table) Results table with columns:
%                     .Type, .Description, .Estimate, .SE, .CI95, .Statistic,
%                     .DF, .pVal, .pAdj, .HVec.
%
%   See also: LME_FIT, FITLME, FITGLME

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'mdl', @(x) isa(x, 'LinearMixedModel') || isa(x, 'GeneralizedLinearMixedModel'));
addParameter(p, 'contrasts', 'all');
addParameter(p, 'correction', 'holm', @ischar);
addParameter(p, 'dfMethod', 'Satterthwaite', @ischar);

parse(p, mdl, varargin{:});
contrastReq = p.Results.contrasts;
correctionMethod = lower(p.Results.correction);
dfMethod = p.Results.dfMethod;

% GLME does not support Satterthwaite
if isa(mdl, 'GeneralizedLinearMixedModel') && strcmpi(dfMethod, 'Satterthwaite')
    dfMethod = 'Residual';
end

%% ========================================================================
%  EXTRACT MODEL DETAILS
%  ========================================================================

coefNames = mdl.CoefficientNames;
nCoefs = length(coefNames);
coefEst = mdl.Coefficients.Estimate;
coefMap = containers.Map(coefNames, 1:nCoefs);
coefCov = mdl.CoefficientCovariance;

% Recalculate fixed effects details (for DF estimation if applicable)
[~, ~, lmeCoefTbl] = fixedEffects(mdl, 'DFMethod', dfMethod);

% Reconstruct Factor Information from Model Variables
% This is crucial to know levels and reference levels used in Dummy Coding.
factorInfo = struct();
dataProps = mdl.Variables.Properties.VariableNames;
frmlFactors = {};
fixedTerms = mdl.Formula.FELinearFormula.TermNames;

for iPred = 1:length(mdl.PredictorNames)
    varName = mdl.PredictorNames{iPred};
    if ismember(varName, dataProps) && iscategorical(mdl.Variables.(varName))
        cats = categories(mdl.Variables.(varName));
        if ~isempty(cats)
            factorInfo.(varName).Levels = cats;
            factorInfo.(varName).RefLevel = cats{1}; % Assuming default dummy coding

            % Check if this factor is part of Fixed Effects
            if any(strcmp(varName, fixedTerms)) || ...
                    any(contains(fixedTerms, [varName ':'])) || ...
                    any(contains(fixedTerms, [':' varName]))
                frmlFactors{end+1} = varName; %#ok<AGROW>
            end
        end
    end
end
frmlFactors = unique(frmlFactors);


%% ========================================================================
%  DEFINE TESTS (ANOVA, COEFF, SIMPLE, MARGINAL)
%  ========================================================================

intDefList = struct(...
    'DefIdx', {}, 'Type', {}, 'Description', {}, ...
    'CoefCombo', {},'HVec', {}, 'OrigCoefIdx', {});
defIdxCounter = 0;

% --- 1. ANOVA (Interaction Terms) ---
% Find highest-order interaction terms using simple colon count
anovaTbl = anova(mdl, 'DFMethod', dfMethod);
nColons = count(string(anovaTbl.Term), ':');
if ~isempty(nColons)
    maxColons = max(nColons);
    idxInteractions = find(nColons == maxColons & nColons > 0);

    for iAnova = 1:length(idxInteractions)
        idxTerm = idxInteractions(iAnova);
        defIdxCounter = defIdxCounter + 1;
        termName = anovaTbl.Term{idxTerm};
        intDefList(defIdxCounter) = struct(...
            'DefIdx', defIdxCounter, ...
            'Type', "ANOVA", ...
            'Description', sprintf('ANOVA: %s', termName),...
            'CoefCombo', {{termName}}, ...
            'HVec', {NaN}, ...
            'OrigCoefIdx', idxTerm);
    end
end

% --- 2. Coefficients (Standard Output) ---
for iCoef = 1:nCoefs
    defIdxCounter = defIdxCounter + 1;
    coefName = coefNames{iCoef};

    if strcmp(coefName, '(Intercept)')
        desc = '(Intercept) Mean at Ref Levels';

    elseif ~contains(coefName, ':')
        % Main Effect
        parts = strsplit(coefName, '_');
        factor = parts{1};
        level = strjoin(parts(2:end),'_');

        if isfield(factorInfo, factor)
            refLvl = factorInfo.(factor).RefLevel;
            baseDesc = sprintf('(%s vs %s)', refLvl, level);

            % Add info about other factors being at Ref
            otherFactors = setdiff(frmlFactors, factor);
            atRefDesc = '';
            if ~isempty(otherFactors)
                refLvlParts = cellfun(@(f) factorInfo.(f).RefLevel, otherFactors, 'Uni', false);
                atRefDesc = [' at ' strjoin(refLvlParts, ', ')];
            end
            desc = [baseDesc, atRefDesc];
        else
            desc = coefName; % Fallback for continuous vars
        end

    else
        % Interaction Coefficient
        terms = strsplit(coefName, ':');
        termDescs = cell(size(terms));
        for iTerm = 1:length(terms)
            parts = strsplit(terms{iTerm}, '_');
            factor = parts{1};
            level = strjoin(parts(2:end),'_');
            if isfield(factorInfo, factor)
                termDescs{iTerm} = sprintf('(%s vs %s)', factorInfo.(factor).RefLevel, level);
            else
                termDescs{iTerm} = terms{iTerm};
            end
        end
        desc = strjoin(termDescs, ' * ');
    end

    intDefList(defIdxCounter) = struct(...
        'DefIdx', defIdxCounter, ...
        'Type', "Coeff", ...
        'Description', desc,...
        'CoefCombo', {{coefName}}, ...
        'HVec', {NaN}, ...
        'OrigCoefIdx', iCoef);
end

% --- 3. Simple & Marginal Effects (For 2-way Interactions) ---
interactionTerms = mdl.Formula.FELinearFormula.TermNames(...
    contains(mdl.Formula.FELinearFormula.TermNames, ':'));
processedInteractions = struct();

for iInt = 1:length(interactionTerms)
    factors = strsplit(interactionTerms{iInt}, ':');
    if length(factors) ~= 2, continue; end

    % Only process if both are categorical factors we tracked
    if ~ismember(factors{1}, frmlFactors) || ~ismember(factors{2}, frmlFactors)
        continue;
    end

    intKey = strjoin(sort(factors),'_x_');
    if isfield(processedInteractions, intKey), continue; end
    processedInteractions.(intKey) = true;

    factorA = factors{1};
    factorB = factors{2};
    lvlsA = factorInfo.(factorA).Levels; refA = factorInfo.(factorA).RefLevel; nA = length(lvlsA);
    lvlsB = factorInfo.(factorB).Levels; refB = factorInfo.(factorB).RefLevel; nB = length(lvlsB);

    % > Simple Effects
    % Effect of A at each level of B
    for iLvlB = 1:nB
        lvlB = lvlsB{iLvlB};
        for iLvlA = 1:nA
            lvlA = lvlsA{iLvlA};
            if strcmp(lvlA, refA), continue; end

            [hVec, coefNamesInv, validH] = contrast_simple(factorA, lvlA, refA, factorB, lvlB, refB, coefMap);

            if validH
                defIdxCounter = defIdxCounter + 1;
                desc = sprintf('(%s vs %s) at %s', refA, lvlA, lvlB);
                intDefList(defIdxCounter) = struct('DefIdx', defIdxCounter, ...
                    'Type', "Simple", 'Description', desc, ...
                    'CoefCombo', {unique(coefNamesInv)}, 'HVec', {hVec}, 'OrigCoefIdx', NaN);
            end
        end
    end

    % Effect of B at each level of A
    for iLvlA = 1:nA
        lvlA = lvlsA{iLvlA};
        for iLvlB = 1:nB
            lvlB = lvlsB{iLvlB};
            if strcmp(lvlB, refB), continue; end

            [hVec, coefNamesInv, validH] = contrast_simple(factorB, lvlB, refB, factorA, lvlA, refA, coefMap);

            if validH
                defIdxCounter = defIdxCounter + 1;
                desc = sprintf('(%s vs %s) at %s', refB, lvlB, lvlA);
                intDefList(defIdxCounter) = struct('DefIdx', defIdxCounter, ...
                    'Type', "Simple", 'Description', desc, ...
                    'CoefCombo', {unique(coefNamesInv)}, 'HVec', {hVec}, 'OrigCoefIdx', NaN);
            end
        end
    end

    % > Marginal Effects
    wB = 1/nB; wA = 1/nA;

    % Marginal Effect of A over B
    for iLvlA = 1:nA
        lvlA = lvlsA{iLvlA};
        if strcmp(lvlA, refA), continue; end

        [hVec, coefNamesInv, validH] = contrast_marginal(factorA, lvlA, factorB, lvlsB, refB, wB, coefMap);

        if validH
            defIdxCounter = defIdxCounter + 1;
            desc = sprintf('(%s vs %s) over %s', refA, lvlA, factorB);
            intDefList(defIdxCounter) = struct('DefIdx', defIdxCounter, ...
                'Type', "Marginal", 'Description', desc, ...
                'CoefCombo', {unique(coefNamesInv)}, 'HVec', {hVec}, 'OrigCoefIdx', NaN);
        end
    end

    % Marginal Effect of B over A
    for iLvlB = 1:nB
        lvlB = lvlsB{iLvlB};
        if strcmp(lvlB, refB), continue; end

        [hVec, coefNamesInv, validH] = contrast_marginal(factorB, lvlB, factorA, lvlsA, refA, wA, coefMap);

        if validH
            defIdxCounter = defIdxCounter + 1;
            desc = sprintf('(%s vs %s) over %s', refB, lvlB, factorA);
            intDefList(defIdxCounter) = struct('DefIdx', defIdxCounter, ...
                'Type', "Marginal", 'Description', desc, ...
                'CoefCombo', {unique(coefNamesInv)}, 'HVec', {hVec}, 'OrigCoefIdx', NaN);
        end
    end
end


%% ========================================================================
%  CALCULATE STATISTICS
%  ========================================================================

nDefs = length(intDefList);
resData = cell(nDefs, 8); % Type, Desc, Est, SE, CI95, Stat, DF, pVal
hVecsToStore = cell(nDefs, 1);

for iDef = 1:nDefs
    rowDef = intDefList(iDef);
    est=NaN; se=NaN; ci95={NaN}; stat=NaN; dfOut={NaN}; pVal=NaN; hVec=NaN;

    switch rowDef.Type
        case "ANOVA"
            anovaRowIdx = rowDef.OrigCoefIdx;
            stat    = anovaTbl.FStat(anovaRowIdx);
            df1     = round(anovaTbl.DF1(anovaRowIdx));
            df2     = round(anovaTbl.DF2(anovaRowIdx));
            dfOut   = {mat2str([df1, df2])};
            pVal    = anovaTbl.pValue(anovaRowIdx);

        case "Coeff"
            origCoefIdx = rowDef.OrigCoefIdx;
            est     = lmeCoefTbl.Estimate(origCoefIdx);
            se      = lmeCoefTbl.SE(origCoefIdx);
            ci95    = [lmeCoefTbl.Lower(origCoefIdx), lmeCoefTbl.Upper(origCoefIdx)];
            ci95    = {mat2str(round(ci95, 2))};
            stat    = lmeCoefTbl.tStat(origCoefIdx);
            dfOut   = {mat2str(round(lmeCoefTbl.DF(origCoefIdx)))};
            pVal    = lmeCoefTbl.pValue(origCoefIdx);

        case {"Simple", "Marginal"}
            hVec = rowDef.HVec;
            [pVal, FVal, ~, dfTest] = coefTest(mdl, hVec, 0, 'DFMethod', dfMethod);

            est = hVec * coefEst(:);
            varContrast = hVec * coefCov * hVec';
            if varContrast < 0 && abs(varContrast) < 1e-10
                varContrast = 0;
            end
            se = sqrt(varContrast);
            stat = sqrt(FVal) * sign(est); % t-stat approximation
            dfOut = {mat2str(round(dfTest))};
    end
    resData(iDef,:) = {rowDef.Type, rowDef.Description, est, se, ci95, stat, dfOut, pVal};
    hVecsToStore{iDef} = hVec;
end


%% ========================================================================
%  FORMAT & FILTER TABLE
%  ========================================================================

lmeStats = cell2table(resData, 'VariableNames', ...
    {'Type','Description','Estimate','SE','CI95','Statistic','DF','pVal'});

% Format H-Vectors
formattedHVecs = cell(nDefs, 1);
for iDef = 1:nDefs
    hCont = hVecsToStore{iDef};
    if isnumeric(hCont) && isvector(hCont)
        formattedHVecs{iDef} = mat2str(round(hCont,3));
    else
        formattedHVecs{iDef} = NaN;
    end
end
lmeStats.HVec = formattedHVecs;

% Add Index
lmeStats = addvars(lmeStats, (1:height(lmeStats))', 'Before', 1, 'NewVariableNames', 'Index');

% Apply Filtering
if isnumeric(contrastReq) || islogical(contrastReq)
    lmeStats = lmeStats(contrastReq, :);
end

% Rounding
vars2 = {'Estimate', 'SE', 'Statistic'};
for i = 1:length(vars2), lmeStats.(vars2{i}) = round(lmeStats.(vars2{i}), 2); end
lmeStats.pVal = round(lmeStats.pVal, 4);


%% ========================================================================
%  MULTIPLE COMPARISON CORRECTION
%  ========================================================================

idxCntrsts = find(ismember(lmeStats.Type, ["Simple", "Marginal"]));
pValCntrsts = lmeStats.pVal(idxCntrsts);
nCntrsts = length(pValCntrsts);

pAdjFull = nan(height(lmeStats), 1);

if nCntrsts > 0 && ~strcmpi(correctionMethod, 'none')
    [pVal_sorted, sortIdx] = sort(pValCntrsts);
    restoreIdx(sortIdx) = 1:nCntrsts;

    switch correctionMethod
        case 'bonferroni'
            pAdj_sorted = min(1, pVal_sorted * nCntrsts);
        case 'holm'
            adjFactor = (nCntrsts:-1:1)';
            pAdj_sorted = min(1, cummax(pVal_sorted .* adjFactor));
        case 'fdr'
            iRank = (1:nCntrsts)';
            pAdj_sorted = min(1, cummin(pVal_sorted .* nCntrsts ./ iRank, 'reverse'));
        otherwise
            pAdj_sorted = pVal_sorted;
    end
    pAdjFull(idxCntrsts) = round(pAdj_sorted(restoreIdx), 4);
end

lmeStats = addvars(lmeStats, pAdjFull, 'After', 'pVal', 'NewVariableNames', 'pAdj');

end % EOF


%% ========================================================================
%  HELPERS
%  ========================================================================

function [hVec, coefNamesInv, validH] = contrast_simple(fMain, lvlMain, ~, fCond, lvlCond, refCond, coefMap)
% Effect of Main at LevelMain, conditioned on Cond at LevelCond

hVec = zeros(1, coefMap.Count);
coefNamesInv = {};
validH = false;

coefMainStr = sprintf('%s_%s', fMain, lvlMain);
if isKey(coefMap, coefMainStr)
    hVec(coefMap(coefMainStr)) = 1;
    coefNamesInv{end+1} = coefMainStr;
    validH = true;

    if ~strcmp(lvlCond, refCond)
        % Add Interaction Term
        coefInt = sprintf('%s_%s:%s_%s', fMain, lvlMain, fCond, lvlCond);
        coefIntAlt = sprintf('%s_%s:%s_%s', fCond, lvlCond, fMain, lvlMain);

        if isKey(coefMap, coefInt)
            hVec(coefMap(coefInt)) = 1;
            coefNamesInv{end+1} = coefInt;
        elseif isKey(coefMap, coefIntAlt)
            hVec(coefMap(coefIntAlt)) = 1;
            coefNamesInv{end+1} = coefIntAlt;
        else
            validH = false; % Interaction required but missing
        end
    end
end
end

function [hVec, coefNamesInv, validH] = contrast_marginal(fMain, lvlMain, fAvg, lvlsAvg, refAvg, wAvg, coefMap)
% Marginal Effect of Main at LevelMain, averaged over fAvg

hVec = zeros(1, coefMap.Count);
coefNamesInv = {};
validH = false;

coefMainStr = sprintf('%s_%s', fMain, lvlMain);
if isKey(coefMap, coefMainStr)
    hVec(coefMap(coefMainStr)) = 1;
    coefNamesInv{end+1} = coefMainStr;
    validH = true;

    for i = 1:length(lvlsAvg)
        lAvg = lvlsAvg{i};
        if strcmp(lAvg, refAvg), continue; end

        coefInt = sprintf('%s_%s:%s_%s', fMain, lvlMain, fAvg, lAvg);
        coefIntAlt = sprintf('%s_%s:%s_%s', fAvg, lAvg, fMain, lvlMain);

        if isKey(coefMap, coefInt)
            hVec(coefMap(coefInt)) = hVec(coefMap(coefInt)) + wAvg;
            coefNamesInv{end+1} = coefInt;
        elseif isKey(coefMap, coefIntAlt)
            hVec(coefMap(coefIntAlt)) = hVec(coefMap(coefIntAlt)) + wAvg;
            coefNamesInv{end+1} = coefIntAlt;
        else
            validH = false;
            break;
        end
    end
end
end





