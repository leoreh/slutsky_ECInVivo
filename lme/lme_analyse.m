function [lmeStats, lmeMdl] = lme_analyse(lmeData, frml, lmeCfg, varargin)

% Fits a user-specified Linear Mixed-Effects (LME) or Generalized
% Linear Mixed-Effects (GLME) model using `fitlme` or `fitglme`.
% The chosen model is then analyzed: overall F-test for the highest-order interaction
% term(s) with anova, Model coefficients, Simple and marginal Effects (for
% 2-way interactions) are generated and tested while adjusting for multiple
% comparisons. See end of function for details.
%
% INPUTS:
%   lmeData         (Required) Table: Data for LME/GLME model.
%   frml            (Required) Char/string: Formula for `fitlme`/`fitglme`.
%   lmeCfg          (Optional) Struct: Configuration with fields (all optional):
%                     .contrasts - Controls which "Coeff", "Simple", "Marginal" effects are in `LME_RESULTS`.
%                     .correction - Multiple comparison correction for Simple/Marginal effects.
%                     .fitMethod - Method for `fitlme`/`fitglme`.
%                     .dfMethod - Method for degrees of freedom estimation.
%                     .distribution - Distribution of the response variable.
%                     .flgG - Force GLME (true/false). If present, overrides distribution logic.
%   ...             (Optional) Name-value pairs for any of the above fields (varargin takes precedence over lmeCfg).
%
% OUTPUTS:
%   LME_RESULTS     Table: Comprehensive analysis results:
%     .Index        - Row index.
%     .Type         - 'ANOVA', 'Coeff', 'Simple', or 'Marginal'.
%     .Description  - Description of the test/coefficient.
%     .Estimate     - Estimate (NaN for 'ANOVA').
%     .SE           - Standard Error (NaN for 'ANOVA').
%     .CI95         - [Lower Upper] 95% Confidence Interval (for 'Coeff' only,
%                     [NaN NaN] for others).
%     .Statistic    - F-statistic ('ANOVA') or t-statistic (others).
%     .DF           - Degrees of freedom: '[NumDF, DenDF]' string for F-tests;
%                     scalar Satterthwaite DF for t-tests (in cell).
%     .pVal         - Uncorrected p-value.
%     .pAdj         - Corrected p-value ('Simple'/'Marginal' only; NaN others).
%     .HVec         - String of H-vector ('Simple'/'Marginal' only; NaN others).
%
%   LME_MDL         Fitted LME/GLME model object.
%
% DEPENDENCIES:
%   MATLAB Statistics and Machine Learning Toolbox.
%
% 250513 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PARSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addRequired(p, 'lmeData', @istable);
addRequired(p, 'frml', @(x) ischar(x) || isstring(x));
addOptional(p, 'lmeCfg', struct(), @isstruct);
addParameter(p, 'contrasts', [], @(x) true);
addParameter(p, 'correction', [], @(x) true);
addParameter(p, 'fitMethod', [], @(x) true);
addParameter(p, 'dfMethod', [], @(x) true);
addParameter(p, 'distribution', [], @(x) true);
addParameter(p, 'flgG', [], @(x) true);
parse(p, lmeData, frml, lmeCfg, varargin{:});

% Simple helper function
function val = getVal(field, default)
    if ~isempty(p.Results.(field))
        val = p.Results.(field);
    elseif isfield(lmeCfg, field) && ~isempty(lmeCfg.(field))
        val = lmeCfg.(field);
    else
        val = default;
    end
end

% Extract values
frml = char(p.Results.frml); % Ensure frml is char
contrastReq = getVal('contrasts', 'all');
correctionMethod = getVal('correction', 'holm');
correctionMethod = lower(char(correctionMethod));  % Ensure lowercase for switch statement
fitMethod = getVal('fitMethod', 'REML');
dfMethod = getVal('dfMethod', 'Satterthwaite');
distribution = getVal('distribution', 'normal');
flgG = getVal('flgG', []);

% Determine if using GLME
if ~isempty(flgG)
    flgG = logical(flgG);
else
    flgG = ~isempty(distribution) && ~strcmpi(distribution, 'normal');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LME/GLME MODEL FITTING & SELECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check response variable distribution fit
resVar = strtrim(strsplit(frml, '~'));
resData = lmeData.(resVar{1});
resSkew = skewness(resData);
thrSkew = 1.5;
if resSkew > thrSkew && strcmp(distribution, 'normal')
    warning('Response variable "%s" does not fit %s distribution', ...
        resVar{1}, distribution);
end

% Fits a single LME/GLME model using the formula provided in frml.
if flgG
    fitMethod = 'REMPL';
    dfMethod = 'Residual';
    lmeMdl = fitglme(lmeData, frml, 'Distribution', distribution, 'FitMethod', fitMethod);
else
    lmeMdl = fitlme(lmeData, frml, 'FitMethod', fitMethod);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTRACT DETAILS FROM CHOSEN MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Retrieves coefficients, factor information, and other details from the mdl.
coefNames = lmeMdl.CoefficientNames;
nCoefs = length(coefNames);
coefEst = lmeMdl.Coefficients.Estimate;
coefMap = containers.Map(coefNames, 1:nCoefs);
coefCov = lmeMdl.CoefficientCovariance;

% recalculate coefficients to implement DF estimation
[~, ~, lmeCoefTbl] = fixedEffects(lmeMdl, 'DFMethod', dfMethod);

factorInfo = struct();
frmlFactors = {};
fixedTerms = lmeMdl.Formula.FELinearFormula.TermNames;
for iPred = 1:length(lmeMdl.PredictorNames)
    varName = lmeMdl.PredictorNames{iPred};
    if ismember(varName, lmeData.Properties.VariableNames) && iscategorical(lmeData.(varName))
        cats = categories(lmeData.(varName));
        if ~isempty(cats)
            factorInfo.(varName).Levels = cats;
            factorInfo.(varName).RefLevel = cats{1};
            if any(strcmp(varName, fixedTerms)) || ...
               any(contains(fixedTerms, [varName ':'])) || ...
               any(contains(fixedTerms, [':' varName]))
                frmlFactors{end+1} = varName;
            end
        end
    end
end
frmlFactors = unique(frmlFactors);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE ANOVA INTERACTION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If the model formula includes interactions, performs ANOVA and extracts
% F-test results for the highest-order interaction term(s).
intDefList = struct(...
    'DefIdx', {}, 'Type', {}, 'Description', {}, ...
    'CoefCombo', {},'HVec', {}, 'OrigCoefIdx', {});
defIdxCounter = 0;

if contains(frml, '*')
    
    % Find highest-order interaction terms (simplistic: those with most colons)
    anovaTbl = anova(lmeMdl, 'DFMethod', dfMethod);
    nColons = count(string(anovaTbl.Term), ':');
    maxColons = max(nColons);
    idxInteractions = find(nColons == maxColons);
    
    for iAnova = 1:length(idxInteractions)
        idxTerm = idxInteractions(iAnova);
        defIdxCounter = defIdxCounter + 1;
        termName = anovaTbl.Term{idxTerm};
        coefDscrpt = sprintf('ANOVA: %s', termName);
        intDefList(defIdxCounter) = struct(...
            'DefIdx', defIdxCounter, 'Type', "ANOVA", 'Description', coefDscrpt,...
            'CoefCombo', {{termName}}, 'HVec', {NaN}, 'OrigCoefIdx', idxTerm);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE COEFFICIENT, SIMPLE & MARGINAL EFFECT DEFINITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Appends definitions for "Coeff" effects
for iCoef = 1:nCoefs
    defIdxCounter = defIdxCounter + 1;
    coefName = coefNames{iCoef};
    if strcmp(coefName, '(Intercept)')
        coefDscrpt = '(Intercept) Mean at Ref Levels';
        
        % fixed effect
    elseif ~contains(coefName, ':')
        parts = strsplit(coefName, '_'); factor = parts{1}; level = strjoin(parts(2:end),'_');
        
        if isfield(factorInfo, factor)
            refLvl = factorInfo.(factor).RefLevel;
            baseDesc = sprintf('(%s vs %s)', refLvl, level);
            otherFactors = setdiff(frmlFactors, factor); atRefDesc = '';
            if ~isempty(otherFactors)
                refLvlParts = cellfun(@(f) factorInfo.(f).RefLevel, otherFactors, 'Uni', false);
                atRefDesc = [' at ' strjoin(refLvlParts, ', ')];
            end
            coefDscrpt = [baseDesc, atRefDesc];
        end
    
        % interaction
    else        
        terms = strsplit(coefName, ':'); termDescs = cell(size(terms));
        for iTerm = 1:length(terms)
            parts = strsplit(terms{iTerm}, '_'); factor = parts{1}; level = strjoin(parts(2:end),'_');
            termDescs{iTerm} = sprintf('(%s vs %s)', factorInfo.(factor).RefLevel, level);
        end
        coefDscrpt = strjoin(termDescs, ' * ');
    end
    
    intDefList(defIdxCounter) = struct(...
        'DefIdx', defIdxCounter, 'Type', "Coeff", 'Description', coefDscrpt,...
        'CoefCombo', {{coefName}}, 'HVec', {NaN}, 'OrigCoefIdx', iCoef);
end

% Appends definitions for "Simple" and "Marginal" effects.
interactionTerms = lmeMdl.Formula.FELinearFormula.TermNames(...
    contains(lmeMdl.Formula.FELinearFormula.TermNames, ':'));
processedInteractions = struct();
for iInt = 1:length(interactionTerms)
    factors = strsplit(interactionTerms{iInt}, ':');
    if length(factors) ~= 2
        continue
    end

    intKey = strjoin(sort(factors),'_x_');
    if isfield(processedInteractions, intKey) || ...
       ~ismember(factors{1}, frmlFactors) || ~ismember(factors{2}, frmlFactors)
        continue;
    end
    processedInteractions.(intKey) = true;

    factorA = factors{1};
    factorB = factors{2};
    lvlsA = factorInfo.(factorA).Levels;
    refA = factorInfo.(factorA).RefLevel;
    nA = length(lvlsA);
    lvlsB = factorInfo.(factorB).Levels;
    refB = factorInfo.(factorB).RefLevel;
    nB = length(lvlsB);

    % Simple Effects
    for iLvlB = 1:nB % Effect of A at each level of B
        lvlB = lvlsB{iLvlB};
        for iLvlA = 1:nA
            lvlA = lvlsA{iLvlA};
            if strcmp(lvlA, refA)
                continue
            end
            hVec = zeros(1, nCoefs);
            coefNamesInv = {};
            validH = false;
            coefMainA = sprintf('%s_%s', factorA, lvlA);

            if isKey(coefMap, coefMainA)
                hVec(coefMap(coefMainA)) = 1;
                coefNamesInv{end+1} = coefMainA;
                validH = true;
                coefDscrpt = sprintf('(%s vs %s) at %s', refA, lvlA, lvlB);
                if ~strcmp(lvlB, refB)
                    coefInt = sprintf('%s_%s:%s_%s', factorA, lvlA, factorB, lvlB);
                    coefIntAlt = sprintf('%s_%s:%s_%s', factorB, lvlB, factorA, lvlA);
                    if isKey(coefMap, coefInt)
                        hVec(coefMap(coefInt)) = 1;
                        coefNamesInv{end+1} = coefInt;
                    elseif isKey(coefMap, coefIntAlt)
                        hVec(coefMap(coefIntAlt)) = 1;
                        coefNamesInv{end+1} = coefIntAlt;
                    else
                        validH = false;
                    end
                end
            end
            if validH && ~(sum(hVec==1)==1 && sum(hVec==0)==length(hVec)-1) % Not duplicate of a Coeff
                defIdxCounter = defIdxCounter + 1;
                intDefList(defIdxCounter) = struct('DefIdx', defIdxCounter,...
                    'Type', "Simple", 'Description', coefDscrpt,...
                    'CoefCombo', {unique(coefNamesInv)}, 'HVec', {hVec},...
                    'OrigCoefIdx', NaN);
            end
        end
    end

    for iLvlA = 1:nA % Effect of B at each level of A
        lvlA = lvlsA{iLvlA};
        for iLvlB = 1:nB
            lvlB = lvlsB{iLvlB}; if strcmp(lvlB, refB); continue; end
            hVec = zeros(1, nCoefs); coefNamesInv = {}; validH = false;
            coefMainB = sprintf('%s_%s', factorB, lvlB);
            if isKey(coefMap, coefMainB)
                hVec(coefMap(coefMainB)) = 1; coefNamesInv{end+1} = coefMainB; validH = true;
                coefDscrpt = sprintf('(%s vs %s) at %s', refB, lvlB, lvlA);
                if ~strcmp(lvlA, refA)
                    coefInt = sprintf('%s_%s:%s_%s', factorA, lvlA, factorB, lvlB);
                    coefIntAlt = sprintf('%s_%s:%s_%s', factorB, lvlB, factorA, lvlA);
                    if isKey(coefMap, coefInt); hVec(coefMap(coefInt)) = 1; coefNamesInv{end+1} = coefInt;
                    elseif isKey(coefMap, coefIntAlt)
                        hVec(coefMap(coefIntAlt)) = 1;
                        coefNamesInv{end+1} = coefIntAlt;
                    else
                        validH = false;
                    end
                end
            end
            if validH && ~(sum(hVec==1)==1 && sum(hVec==0)==length(hVec)-1)
                defIdxCounter = defIdxCounter+1;
                intDefList(defIdxCounter) = struct('DefIdx', defIdxCounter,...
                    'Type', "Simple", 'Description', coefDscrpt,...
                    'CoefCombo', {unique(coefNamesInv)},...
                    'HVec', {hVec}, 'OrigCoefIdx', NaN);
            end
        end
    end
    
    % Marginal Effects
    wB = 1/nB; wA = 1/nA;
    for iLvlA = 1:nA % Marginal Effect of A over B
        lvlA = lvlsA{iLvlA};
        if strcmp(lvlA, refA)
            continue
        end
        hVec = zeros(1, nCoefs);
        coefNamesInv = {};
        validH = false;
        coefMainA = sprintf('%s_%s', factorA, lvlA);
        if isKey(coefMap, coefMainA)
            hVec(coefMap(coefMainA)) = 1;
            coefNamesInv{end+1} = coefMainA;
            validH = true;

            for iLvlB = 1:nB
                lvlB = lvlsB{iLvlB};
                if strcmp(lvlB, refB)
                    continue
                end
                coefInt = sprintf('%s_%s:%s_%s', factorA, lvlA, factorB, lvlB);
                coefIntAlt = sprintf('%s_%s:%s_%s', factorB, lvlB, factorA, lvlA);
                if isKey(coefMap, coefInt)
                    hVec(coefMap(coefInt)) = hVec(coefMap(coefInt)) + wB;
                    coefNamesInv{end+1} = coefInt;
                elseif isKey(coefMap, coefIntAlt)
                    hVec(coefMap(coefIntAlt)) = hVec(coefMap(coefIntAlt)) + wB;
                    coefNamesInv{end+1} = coefIntAlt;
                else
                    validH = false;
                    break
                end
            end
        end
        if validH
            defIdxCounter = defIdxCounter + 1;
            coefDscrpt = sprintf('(%s vs %s) over %s', refA, lvlA, factorB);
            intDefList(defIdxCounter) = struct('DefIdx', defIdxCounter,...
                'Type', "Marginal", 'Description', coefDscrpt,...
                'CoefCombo', {unique(coefNamesInv)}, 'HVec', {hVec},...
                'OrigCoefIdx', NaN);
        end
    end

    for iLvlB = 1:nB % Marginal Effect of B over A
        lvlB = lvlsB{iLvlB}; 
        if strcmp(lvlB, refB)
            continue
        end
        hVec = zeros(1, nCoefs);
        coefNamesInv = {};
        validH = false;
        coefMainB = sprintf('%s_%s', factorB, lvlB);
        if isKey(coefMap, coefMainB)
            hVec(coefMap(coefMainB)) = 1;
            coefNamesInv{end+1} = coefMainB;
            validH = true;
            for iLvlA = 1:nA
                lvlA = lvlsA{iLvlA};
                if strcmp(lvlA, refA)
                    continue
                end
                coefInt = sprintf('%s_%s:%s_%s', factorA, lvlA, factorB, lvlB);
                coefIntAlt = sprintf('%s_%s:%s_%s', factorB, lvlB, factorA, lvlA);
                if isKey(coefMap, coefInt)
                    hVec(coefMap(coefInt)) = hVec(coefMap(coefInt)) + wA;
                    coefNamesInv{end+1} = coefInt;
                elseif isKey(coefMap, coefIntAlt)
                    hVec(coefMap(coefIntAlt)) = hVec(coefMap(coefIntAlt)) + wA;
                    coefNamesInv{end+1} = coefIntAlt;
                else
                    validH = false;
                    break
                end
            end
        end
        if validH
            defIdxCounter = defIdxCounter + 1;
            coefDscrpt = sprintf('(%s vs %s) over %s', refB, lvlB, factorA);
            intDefList(defIdxCounter) = struct('DefIdx', defIdxCounter,...
                'Type', "Marginal", 'Description', coefDscrpt,...
                'CoefCombo', {unique(coefNamesInv)}, 'HVec', {hVec},...
                'OrigCoefIdx', NaN);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE STATISTICS FOR ALL DEFINITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iterates through all generated definitions, calculates relevant statistics,
% and stores them. CI95 is added for coefficients.
nDefs = length(intDefList);
resData = cell(nDefs, 8); % Type,Desc,Est,SE,CI95,Stat,DF,pVal
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
            [pVal, FVal, ~, dfTest] = coefTest(lmeMdl, hVec,...
                0, 'DFMethod', dfMethod);
            est = hVec * coefEst(:);
            varContrast = hVec * coefCov * hVec';
            if varContrast < 0 && abs(varContrast) < 1e-10
                varContrast = 0;
            end
            se = sqrt(varContrast);
            stat = sqrt(FVal) * sign(est);       % same as est / se
            dfOut = {mat2str(round(dfTest))};
    end
    resData(iDef,:) = {rowDef.Type,rowDef.Description,est,se,ci95,stat,dfOut,pVal};
    hVecsToStore{iDef} = hVec;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSEMBLE AND FILTER FULL RESULTS TABLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converts data to table, formats, reorders, and applies user's contrast
% selection.

%**************************************************************************
% Create Table
lmeStats = cell2table(resData, 'VariableNames', ...
    {'Type','Description','Estimate','SE','CI95','Statistic','DF','pVal'});

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

%**************************************************************************
% Round numeric columns
clmnRnd2 = {'Estimate', 'SE', 'Statistic'};
clmnRnd4 = {'pVal'};
for iCol = 1:length(clmnRnd2)
    col = clmnRnd2{iCol};
    lmeStats.(col) = round(lmeStats.(col), 2);
end
for iCol = 1:length(clmnRnd4)
    col = clmnRnd4{iCol};
    lmeStats.(col) = round(lmeStats.(col), 4);
end
if isnumeric(lmeStats.CI95)
    lmeStats.CI95 = round(lmeStats.CI95, 2);
end

% Add index for comparisons before applying use selection
lmeStats = addvars(lmeStats, (1:height(lmeStats))', 'Before', 1, 'NewVariableNames', 'Index');

% Apply user's contrast selection ('all', numeric indices, or logical
% vector). Refers to indices in the full generated list.
if isnumeric(contrastReq) || islogical(contrastReq)
    lmeStats = lmeStats(contrastReq, :);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% APPLY MULTIPLE COMPARISON CORRECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collects p-values for "Simple" and "Marginal" effects from the filtered
% lmeStats table and applies the specified correction method.

% Identify 'Simple' or 'Marginal' contrasts
idxCntrsts = find(ismember(lmeStats.Type, ["Simple", "Marginal"]));
pValCntrsts = lmeStats.pVal(idxCntrsts);
nCntrsts = length(pValCntrsts);

if nCntrsts > 0 && ~strcmpi(correctionMethod, 'none')
        [pVal_sorted, sortIdx] = sort(pValCntrsts);
        restoreIdx(sortIdx) = 1:nCntrsts;
        
        % calculate corrected pVal
        switch correctionMethod
            case 'bonferroni'
                pAdj_sorted = min(1, pVal_sorted * nCntrsts);
            case 'holm'
                adjFactor = (nCntrsts:-1:1)';
                adjPSortedRaw = pVal_sorted(:) .* adjFactor;
                pAdj_sorted = min(1,cummax(adjPSortedRaw));
            case 'fdr'
                iRank = (1:nCntrsts)';
                adjPRaw = pVal_sorted(:).*nCntrsts./iRank;
                pAdj_sorted = min(1,cummin(adjPRaw,'reverse'));
            otherwise
                pAdj_sorted = pVal_sorted;
        end
        
        % oragnize output
        pAdj = nan(height(lmeStats), 1);
        pAdj(idxCntrsts) = round(pAdj_sorted(restoreIdx), 4);
        lmeStats = addvars(lmeStats, pAdj, 'After', 'pVal', 'NewVariableNames', 'pAdj');
end
end % EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THEORETICAL CONSIDERATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear Mixed-Effects (LME) models are used for hierarchically structured
% data (e.g., repeated measures within subjects) to account for
% non-independence of observations. They partition variance into fixed and
% random effects. This function fits models using Restricted Maximum
% Likelihood (REML) by default, as opposed to Maximum Likelihood (ML),
% which provides unbiased estimates of variance components.
%
% Fixed Effects', denoted 'beta', represent the average effects of
% predictor variables across the entire population being studied. These are
% usually the primary focus of the analysis. 
%
% 'Random Effects', denoted 'b', capture the variability *between* the
% groups or subjects. They model how individual subjects systematically
% deviate from the average population. For instance, a 'random intercept'
% allows each subject to have their own unique baseline response level
% (e.g., specified as `(1 | Subject)` in the formula). A 'random slope'
% allows the effect of a predictor, like time, to differ across subjects
% (e.g., `(Time | Subject)`). By incorporating random effects, the model
% properly accounts for the data's dependency structure when evaluating the
% fixed effects.
%
% 'Residual Error', denoted as 'epsilon', accounts for the
% remaining, unexplained variability or random noise in the data that isn't
% captured by either the fixed or the random effects structure.

%**************************************************************************
% G/LME
%**************************************************************************
% The choice between using LME and GLME is determined exclusively by
% the distribution of the response variable. LME is used when the
% response variable is continuous and follows a normal distribution, or can
% be transformed to approximate one. GLME (Generalized Linear
% Mixed-Effects Model) is required when the response variable does not
% follow a normal distribution. 

% The decision to transform predictor variables is a separate consideration
% from the choice of model and is applicable to both LME and GLME.
% Predictors should be log-transformed when their distribution is highly
% skewed. The purpose of this transformation is twofold: it makes the
% predictor's distribution more symmetric, and more importantly, it can
% linearize a non-linear relationship between the predictor and the
% response. By converting an exponential or power-law relationship into a
% linear one, the log-transformed predictor better satisfies the
% fundamental assumption of linear models. This practice reduces the
% disproportionate influence of extreme data points and often leads to a
% more accurate and stable model fit.
% 
% Z-scoring predictors is performed not to address distribution shape but
% to standardize variables onto a common scale. This is particularly useful
% when predictors are measured in different units (e.g., firing rate in Hz,
% time in seconds). After z-scoring, the predictors model coefficients are
% directly comparable, as a larger coefficient now unambiguously indicates
% a stronger effect per standard deviation, regardless of the original
% units.

% TLDR: Z-score and log-transform contineous variables only when they are
% predictors. As response variables, adjust the model to their
% distribution.

%**************************************************************************
% Dummy Coding and Coefficient Types
%**************************************************************************
% MATLAB's `fitlme` and `fitglme` functions use 'Dummy Coding' by default for categorical
% predictors. (also known as 'treatment contrasts'). This means one level
% of each categorical factor is chosen as the 'reference level' (often the
% first alphabetically or numerically). All other levels of that factor are
% then compared directly to this reference level.

% The '(Intercept)' coefficient represents the estimated average value of
% the response variable under a specific baseline condition: when all
% categorical predictors in the model are at their reference levels, and
% simultaneously, all continuous predictors included in the model formula
% are equal to zero. A statistically significant intercept test suggests
% this baseline mean is different from zero.
%
% A 'Main Effect' coefficient, such as one named 'FactorA_LevelX',
% estimates the *difference* in the average response between 'LevelX' of
% 'FactorA' and the reference level of 'FactorA'. However, a critical point
% arises when interactions are present in the model: this comparison is
% valid *only* when all other factors that interact with 'FactorA' are held
% at their respective reference levels. Therefore, in the presence of
% interactions, this coefficient does not represent the overall main effect
% averaged across the levels of other factors.
%
% An 'Interaction' Coefficient, for example
% 'FactorA_LevelX:FactorB_LevelY', assesses how the effect of one factor
% changes depending on the level of another. It quantifies the *additional
% difference* in the response associated with FactorA being at LevelX (vs.
% its reference) when FactorB is specifically at LevelY (compared to when
% FactorB is at its reference level). It essentially captures a "difference
% of differences". A significant interaction term indicates that the effect
% of FactorA is not constant but depends on the specific level of FactorB
% (and vice-versa).

%**************************************************************************
% Contrast Types 
%**************************************************************************
% This function facilitates hypothesis testing about specific comparisons
% or combinations of the fixed-effect coefficients (the betas).
% Mathematically, each test evaluates a hypothesis like H0: H * beta = 0,
% where 'H' is a row vector that defines the linear combination of
% coefficients forming the contrast of interest. The results table
% (`lme_results`) categorizes these tests into distinct types.
%
% The "Coefficient" type refers to the standard tests for individual model
% coefficients equalling zero. The results for these are taken directly
% from the standard output of the LME/GLME model. No separate contrast
% vector 'H' is constructed or tested using `coefTest` here and these
% standard coefficient tests are *not* included in the multiple comparison
% correction procedures; the `pAdj` column will always be NaN for
% rows of Type "Coefficient".
%
% The "Simple Effect" type is generated specifically when two factors
% interact. It tests the effect of one factor at a single, fixed level of
% the other factor. For instance, it might test the difference between
% LevelX and the reference level of FactorA, specifically *only* when
% FactorB is at LevelY. These tests require constructing a specific
% contrast vector 'H' that typically combines the main effect coefficient
% for FactorA_LevelX with the relevant interaction coefficient (e.g.,
% FactorA_LevelX:FactorB_LevelY). The function runs `coefTest` using this H
% vector. The Estimate (calculated as H * beta_hat) and its Standard Error
% (SE, calculated as sqrt(H * Cov(beta_hat) * H')) are computed. These
% simple effect tests undergo multiple comparison correction if they are
% selected for inclusion in the final results table (via the 'contrasts'
% option). The correction is applied collectively across all selected
% simple and marginal effects.
%
% The "Marginal Effect" type is also generated for two-way interactions. It
% tests the effect of one factor averaged across all the levels of the
% other factor it interacts with. This often aligns more closely with the
% intuitive notion of a "main effect" in the presence of an interaction.
% For example, it might test the difference between LevelX and the
% reference level of FactorA, averaged across all levels of FactorB.
% Constructing the H vector for this involves combining the main effect
% coefficient with a weighted sum of the relevant interaction coefficients,
% where the weights ensure averaging (typically 1 divided by the number of
% levels of the factor being averaged over). Again, `coefTest` is used with
% this H vector, the Estimate and SE are computed, and the H vector is
% stored. Like simple effects, these marginal effect tests undergo multiple
% comparison correction if selected for the output.
%
% All results produced by this function – whether standard coefficients,
% simple effects, or marginal effects – pertain strictly to the *fixed
% effects* part of the LME model. They allow inferences about
% *population averages*, representing overall trends or differences in the
% population from which the sample was drawn. These results do not describe
% effects specific to individual subjects or clusters (often called
% conditional effects). Analyzing subject-specific deviations would require
% examining the estimated random effects ('b') themselves, e.g. using
% `randomEffects(lme)`.

%**************************************************************************
% ANOVA 
%**************************************************************************
% The `anova(lme)` method provides F-tests for each term. For a model
% fitted with dummy coding:
%
% Interaction F-test (e.g., for Group:Day): This is a robust test for the
% overall significance of the interaction. It should be examined first.
%
% Main Effect F-tests (e.g., for Group, for Day): These test hypotheses
% about the simple effects at reference levels (due to dummy coding). If the
% interaction is significant, these main effect F-tests are less directly
% interpretable as "averaged effects" and the focus should shift to
% simple/marginal effects that dissect the interaction. If the interaction
% is *not* significant, these F-tests become more interpretable as overall
% main effects.
%
% Note: To obtain Type III F-tests where main effects are tested more akin to
% being "averaged over" other factors in the presence of interactions, `fitlme`/`fitglme`
% would need to be called with `'DummyVarCoding', 'effects'`. This function
% currently uses default dummy coding to simplify H-vector construction for
% simple/marginal effects.

%**************************************************************************
% Estimate +/ SE
%**************************************************************************
% The `Estimate` column provides the model's calculated magnitude for each
% specific term or contrast. For the Intercept, the `Estimate` represents
% the predicted mean value of the response variable when all categorical
% fixed factors are at their designated reference levels and any continuous
% fixed predictors are held at zero. For other coefficients (main effects
% or interactions under dummy coding), the `Estimate` represents the
% calculated difference compared to the relevant reference level(s). For
% instance, a main effect coefficient's Estimate is the difference between
% that factor level and its reference level, specifically evaluated when
% other interacting factors are at their reference levels. An interaction
% coefficient's Estimate quantifies how the effect of one factor changes
% across levels of another. For derived contrasts (like simple or marginal
% effects), the Estimate is the calculated value of the specific linear
% combination of coefficients being tested (e.g., the estimated mean
% difference between two groups at a specific condition). The `SE`
% (Standard Error) accompanies each Estimate and quantifies the precision
% of that estimate; it reflects the expected standard deviation of the
% estimate if the study were repeated many times. The `SE` is crucial for
% assessing statistical significance, as it forms the denominator in the
% t-statistic calculation (t = Estimate / SE) and is used to construct
% confidence intervals around the `Estimate`.
 
%**************************************************************************
% Statistic (t or F)
%**************************************************************************
% The `Statistic` column contains the test statistic associated with the
% hypothesis test for that row. For rows of type 'Coeff', 'Simple', or
% 'Marginal', this is the t-statistic (Estimate / SE), testing whether the
% specific coefficient or contrast is significantly different from zero.
% For rows of type 'ANOVA', this is the F-statistic obtained from the
% `anova(lme)` function, testing the overall significance of the
% interaction term in the model.

%**************************************************************************
% Degrees of Freedom (DF)
%**************************************************************************
% The `DF` column reports the degrees of freedom associated with the test
% statistic. For t-tests ('Coeff', 'Simple', 'Marginal'), this is a single
% value representing the denominator degrees of freedom. For F-tests
% ('ANOVA'), this column displays the numerator and denominator degrees of
% freedom as a string '[DF1, DF2]'. DF1 relates to the number of parameters
% associated with the effect being tested, and DF2 relates to the residual
% or error degrees of freedom. The default in this function is to estimate
% DFs using Satterthwaite approximations. This is important for
% hierarchical data, because observations within the same group (e.g.,
% subject) are not independent, reducing the effective degrees of freedom.
% The default method in `fitlme` often overestimates DFs by treating random
% effects as if they were fixed effects, and calculating DFs as the number
% of observations minus the number of parameters, leading to overly
% optimistic p-values. Satterthwaite's method provides a more accurate
% approximation by accounting for the hierarchical structure of the data,
% the uncertainty in variance component estimates, and the actual
% dependencies between observations. This is particularly important in
% unbalanced designs or when the number of higher-level units (e.g.,
% subjects) is small relative to the total number of observations. While
% Satterthwaite's method may yield smaller DFs than the default method, it
% provides more reliable statistical inference by better reflecting the
% true uncertainty in the parameter estimates.

%**************************************************************************
% Model Fit Statistics 
%**************************************************************************
% LogLikelihood (LL): This reflects the probability of observing the actual
% data given the estimated model parameters. A higher LL indicates a better
% fit, meaning the model makes the observed data more likely. However, LL
% naturally increases as more parameters are added to a model, even if
% those parameters don't meaningfully improve the fit, making it less
% suitable for comparing models with different numbers of parameters on its
% own. 
% 
% Deviance: Calculated as -2 times the LogLikelihood (-2\*LL). Lower
% deviance values correspond to higher log-likelihoods and thus indicate a
% better model fit. Deviance is primarily useful for comparing nested
% models (where one model is a simplified version of the other) using a
% Likelihood Ratio Test (LRT), where the difference in deviance follows a
% chi-squared distribution.
% 
% AIC (Akaike Information Criterion): AIC provides a way to select models
% by estimating the prediction error and thereby the relative quality of
% statistical models for a given set of data. It balances goodness-of-fit
% (related to LL) with model parsimony by adding a penalty for the number
% of estimated parameters (k): AIC = -2*LL + 2*k. Lower AIC values suggest
% a preferred model, indicating a better trade-off between fit and
% complexity, often favoring models with better predictive performance.
% 
% BIC (Bayesian Information Criterion): Similar to AIC, BIC also balances
% fit and complexity but applies a larger penalty for the number of
% parameters, particularly in larger datasets (n): BIC = -2*LL + k*ln(n).
% Lower BIC values are preferred. BIC tends to favor simpler, more
% parsimonious models compared to AIC and is sometimes considered better
% for selecting the "true" model if one is assumed to exist within the set
% of candidates.

%**************************************************************************
% Multiple Comparisons
%**************************************************************************
% Applied only to the family of selected "Simple Effect" and "Marginal
% Effect". Not applied to ANOVA F-tests or individual model coefficients.
%
% To counteract this inflation of the error rate across the set of tests,
% correction methods are applied to the calculated p-values. The
% 'correction' parameter allows choosing between different strategies:
%
% Bonferroni ('bonferroni'): This is the simplest method, controlling the
% Family-Wise Error Rate (FWER) – the probability of making one or more
% Type I errors across all tests performed. It achieves this by multiplying
% each individual p-value by the total number of tests conducted (m) or,
% equivalently, by dividing the target significance level (e.g., 0.05) by
% m. Bonferroni is often overly conservative, significantly reducing
% statistical power and increasing the risk of failing to detect true
% effects (Type II errors).
%
% Holm-Bonferroni ('holm', default): This method also controls the FWER but
% is uniformly more powerful (less conservative) than the standard
% Bonferroni procedure. It works sequentially: p-values are ordered from
% smallest to largest, and the significance threshold is adjusted
% step-by-step. The smallest p-value is compared against alpha/m, the next
% smallest against alpha/(m-1), and so on, stopping when a p-value fails to
% meet its adjusted threshold. This provides the same strong FWER control
% as Bonferroni but with a better chance of detecting true effects.
%
% False Discovery Rate ('fdr', specifically Benjamini-Hochberg): Instead of
% controlling the probability of making any Type I error (FWER), the FDR
% approach controls the expected proportion of rejected null hypotheses
% that are actually false positives. For example, setting FDR control at 5%
% aims to ensure that, among all the effects declared significant, no more
% than 5% are expected to be false discoveries. This method is less
% stringent than FWER control, particularly when a large number of tests
% are performed. Consequently, it offers considerably more power to detect
% true effects, making it suitable for more exploratory analyses where
% controlling the proportion of false findings is deemed acceptable, rather
% than strictly preventing any single false positive.