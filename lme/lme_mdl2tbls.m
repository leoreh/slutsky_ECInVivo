function lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)
% LME_MDL2TBLS Extracts key information from an LME model into a structured 
% array of formatted tables suitable for Excel export.
%
% INPUT:
%   lmeMdl   - A LinearMixedModel object (from fitlme).
%   lmeStats - Statistics table from lme_postHoc.
%   lmeInfo  - Information structure from lme_analyse.
%
% OUTPUT:
%   lmeTbls  - A structured array with fields Title and Table.

lmeTbls = struct('Title', {}, 'Table', {});
tblIdx = 1;

% Helper function to add to structural array
    function add_tbl(titleStr, tbl, tblType)
        if nargin < 3, tblType = 'Other'; end
        if isempty(tbl) || height(tbl) == 0, return; end
        
        tbl = format_table_for_pub(tbl, tblType);
        % Convert arrays to string for export
        tbl = stringify_arrays(tbl);

        lmeTbls(tblIdx).Title = titleStr;
        lmeTbls(tblIdx).Table = tbl;
        tblIdx = tblIdx + 1;
    end

% --- 1. Model Information ---
infoNames = {};
infoVals = {};
if ~isempty(lmeInfo)
    if isfield(lmeInfo, 'frml'), infoNames{end+1}='Formula'; infoVals{end+1}=lmeInfo.frml; end
    if isfield(lmeInfo, 'distFinal'), infoNames{end+1}='Distribution'; infoVals{end+1}=lmeInfo.distFinal; end
    if isfield(lmeInfo, 'fitMethod'), infoNames{end+1}='Fit Method'; infoVals{end+1}=lmeInfo.fitMethod; end
    if isfield(lmeInfo, 'dfMethod'), infoNames{end+1}='DF Method'; infoVals{end+1}=lmeInfo.dfMethod; end
    if isfield(lmeInfo, 'aic'), infoNames{end+1}='AIC'; infoVals{end+1}=mat2str(round(lmeInfo.aic, 2)); end
end

if ~isempty(lmeMdl) && isprop(lmeMdl, 'NumObservations')
    infoNames{end+1}='Observations'; infoVals{end+1}=mat2str(lmeMdl.NumObservations);
end
if ~isempty(infoNames)
    Property = infoNames';
    Value = infoVals';
    add_tbl('MODEL INFORMATION', table(Property, Value), 'Other');
end

% --- 2. Continuous Predictors Table ---
if ~isempty(lmeInfo) && isfield(lmeInfo, 'transParams') && isfield(lmeInfo.transParams, 'varsTrans')
    predNames = fieldnames(lmeInfo.transParams.varsTrans);
    numPredictors = length(predNames);
    
    if numPredictors > 0
        isStnd = false(numPredictors, 1);
        transMethod = cell(numPredictors, 1);
        
        for i = 1:numPredictors
            vName = predNames{i};
            pVar = lmeInfo.transParams.varsTrans.(vName);
            
            % Check standardization status
            if isfield(pVar, 'flgZ')
                isStnd(i) = pVar.flgZ;
            elseif isfield(lmeInfo, 'standardized')
                isStnd(i) = lmeInfo.standardized; % Fallback
            end
            
            % Check Transformation status
            if isfield(pVar, 'logBase')
                lb = pVar.logBase;
                if ischar(lb) && strcmpi(lb, 'logit')
                    tMethod = 'Logit';
                elseif ischar(lb) && strcmpi(lb, 'e')
                    tMethod = 'Natural Log';
                elseif isnumeric(lb) && ~isempty(lb)
                    tMethod = sprintf('Log%g', lb);
                else
                    tMethod = 'None';
                end
            else
                tMethod = 'None';
            end
            
            % Check Offset
            if isfield(pVar, 'offset') && pVar.offset > 0
                if strcmp(tMethod, 'None')
                    tMethod = sprintf('Offset: +%g', round(pVar.offset, 4));
                else
                    tMethod = sprintf('%s (Offset: +%g)', tMethod, round(pVar.offset, 4));
                end
            end
            
            transMethod{i} = tMethod;
        end
        
        keepPred = true(numPredictors, 1);
        if ~isempty(lmeMdl) && isprop(lmeMdl, 'ResponseName')
            keepPred = ~strcmp(predNames, lmeMdl.ResponseName);
        end
        
        predNames = predNames(keepPred);
        isStnd = isStnd(keepPred);
        transMethod = transMethod(keepPred);
        
        if ~isempty(predNames)
            predTbl = table(predNames, isStnd, transMethod, ...
                'VariableNames', {'Predictor', 'Standardized', 'Transformation'});
            add_tbl('CONTINUOUS PREDICTORS', predTbl, 'Other');
        end
    end
end

% --- 3. Fixed Effects / Coefficients (Type == "Coeff") ---
if ~isempty(lmeStats) && ismember('Type', lmeStats.Properties.VariableNames)
    idxCoef = lmeStats.Type == "Coeff";
    if any(idxCoef)
        coefTbl = lmeStats(idxCoef, :);
        coefTbl = removevars_safely(coefTbl, {'Index', 'Type', 'HVec'});
        add_tbl('FIXED EFFECTS (COEFFICIENTS)', coefTbl, 'Coeff');
    end
end

% --- 4. Post-Hoc (Type == "Simple" or "Marginal") ---
if ~isempty(lmeStats) && ismember('Type', lmeStats.Properties.VariableNames)
    idxPost = ismember(lmeStats.Type, ["Simple", "Marginal"]);
    if any(idxPost)
        postTbl = lmeStats(idxPost, :);
        postTbl = removevars_safely(postTbl, {'Index', 'Type', 'HVec'});
        add_tbl('POST-HOC EFFECTS', postTbl, 'PostHoc');
    end
end

% --- 5. Conditional ANOVA ---
if ~isempty(lmeStats) && ismember('Type', lmeStats.Properties.VariableNames)
    idxAnova = find(lmeStats.Type == "ANOVA");
    if ~isempty(idxAnova) && ~isempty(lmeMdl)
        includeAnova = false;
        dataProps = lmeMdl.Variables.Properties.VariableNames;
        anovaTerms = lmeStats.Description(idxAnova);
        
        for iT = 1:length(anovaTerms)
            termStr = char(anovaTerms{iT});
            termStr = strrep(termStr, 'ANOVA: ', '');
            factors = strtrim(strsplit(termStr, ':'));
            
            for iF = 1:length(factors)
                varName = factors{iF};
                if ismember(varName, dataProps)
                    varData = lmeMdl.Variables.(varName);
                    if iscategorical(varData)
                        numLvls = numel(categories(varData));
                        if numLvls > 2, includeAnova = true; break; end
                    end
                end
            end
            if includeAnova, break; end 
        end
        
        if includeAnova
            anovaTbl = lmeStats(idxAnova, :);
            anovaTbl = removevars_safely(anovaTbl, {'Index', 'Type', 'Estimate', 'SE', 'CI95', 'HVec'});
            add_tbl('ANOVA (OMNIBUS INTERACTION)', anovaTbl, 'ANOVA');
        end
        
        % Other Unclassified Effects (not ANOVA, Coef, Simple, Marginal)
        idxOther = ~(lmeStats.Type == "ANOVA" | lmeStats.Type == "Coeff" | ismember(lmeStats.Type, ["Simple", "Marginal"]));
        if any(idxOther)
            otherTbl = lmeStats(idxOther, :);
            otherTbl = removevars_safely(otherTbl, {'Index', 'HVec'});
            add_tbl('OTHER LME EFFECTS', otherTbl, 'Other');
        end
    elseif ~ismember('Type', lmeStats.Properties.VariableNames)
        add_tbl('ALL LME STATISTICS', lmeStats, 'Other');
    end
end

% --- 6. Covariance & Error Parameters ---
if ~isempty(lmeMdl)
    try
        [~, ~, resStats] = covarianceParameters(lmeMdl);
        add_tbl('COVARIANCE PARAMETERS', dataset2table(resStats{1}), 'Other');
        add_tbl('ERROR PARAMETERS', dataset2table(resStats{2}), 'Other');
    catch
    end
end

end % EOF

%% ========================================================================
%  HELPERS
%  ========================================================================

function tbl = stringify_arrays(tbl)
% STRINGIFY_ARRAYS converts cell arrays of numeric vectors to string
tblVars = tbl.Properties.VariableNames;
for v = 1:length(tblVars)
    colName = tblVars{v};
    colData = tbl.(colName);
    if iscell(colData)
        for i = 1:length(colData)
            val = colData{i};
            if isnumeric(val) && isvector(val) && length(val) > 1
                colData{i} = mat2str(round(val, 2));
            elseif isnumeric(val) && length(val) == 1
                colData{i} = mat2str(round(val, 2)); % Or keep as numeric depending on needs
            end
        end
        % Re-assign the column if all elements are strings
        allStrings = all(cellfun(@ischar, colData) | cellfun(@isstring, colData));
        if allStrings
            tbl.(colName) = string(colData);
        else
            tbl.(colName) = colData;
        end
    end
end
end

function tbl = format_table_for_pub(tbl, tableType)
% FORMAT_TABLE_FOR_PUB Rounds continuous columns appropriately and 
% renames them to present a neat, publication-ready table depending on Type.

clmnRnd2 = {'Estimate', 'SE', 'SEPred', 'Statistic', 'tStat', 'FStat', 'Lower', 'Upper'};
clmnRnd4 = {'pValue', 'pVal', 'pAdj'};

for iCol = 1:length(clmnRnd2)
    col = clmnRnd2{iCol};
    if ismember(col, tbl.Properties.VariableNames)
        tbl.(col) = round(tbl.(col), 2);
    end
end

for iCol = 1:length(clmnRnd4)
    col = clmnRnd4{iCol};
    if ismember(col, tbl.Properties.VariableNames)
        tbl.(col) = round(tbl.(col), 4);
    end
end

% Global rename
renameMap = {
    'pVal',      'P-value';
    'pValue',    'P-value';
    'pAdj',      'Adj P-value';
};

for iR = 1:size(renameMap, 1)
    oldName = renameMap{iR, 1};
    newName = renameMap{iR, 2};
    if ismember(oldName, tbl.Properties.VariableNames)
        tbl.Properties.VariableNames{strcmp(tbl.Properties.VariableNames, oldName)} = newName;
    end
end

% Specific modifications based on Table Type requested by user
switch tableType
    case 'ANOVA'
        if ismember('Adj P-value', tbl.Properties.VariableNames)
            tbl = removevars_safely(tbl, {'Adj P-value'});
        end
        if ismember('Statistic', tbl.Properties.VariableNames)
            tbl.Properties.VariableNames{strcmp(tbl.Properties.VariableNames, 'Statistic')} = 'F-statistic';
        end
    case 'Coeff'
        if ismember('Adj P-value', tbl.Properties.VariableNames)
            tbl = removevars_safely(tbl, {'Adj P-value'});
        end
        if ismember('Statistic', tbl.Properties.VariableNames)
            tbl.Properties.VariableNames{strcmp(tbl.Properties.VariableNames, 'Statistic')} = 't-statistic';
        end
    case 'PostHoc'
        if ismember('P-value', tbl.Properties.VariableNames)
            tbl = removevars_safely(tbl, {'P-value'});
        end
        if ismember('Adj P-value', tbl.Properties.VariableNames)
            tbl.Properties.VariableNames{strcmp(tbl.Properties.VariableNames, 'Adj P-value')} = 'P-value (Adj.)';
        end
        if ismember('Statistic', tbl.Properties.VariableNames)
            tbl.Properties.VariableNames{strcmp(tbl.Properties.VariableNames, 'Statistic')} = 't-statistic';
        end
    case 'Other'
        if ismember('Statistic', tbl.Properties.VariableNames)
            tbl.Properties.VariableNames{strcmp(tbl.Properties.VariableNames, 'Statistic')} = 't-statistic';
        end
end
end

function tblOut = removevars_safely(tblIn, varsToRemove)
% REMOVEVARS_SAFELY Removes variables from a table, ignoring those that don't exist.
varsToRemove = varsToRemove(ismember(varsToRemove, tblIn.Properties.VariableNames));
if ~isempty(varsToRemove)
    tblOut = removevars(tblIn, varsToRemove);
else
    tblOut = tblIn;
end
end