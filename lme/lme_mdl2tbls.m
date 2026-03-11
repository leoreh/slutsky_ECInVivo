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
        % Reorder columns globally before stringifying
        tbl = reorder_columns(tbl);
        
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
    if isfield(lmeInfo, 'distSelected') && ~isempty(lmeInfo.distSelected)
        distStr = regexprep(char(lmeInfo.distSelected), '(\<[a-z])', '${upper($1)}');
        infoNames{end+1}='Distribution'; infoVals{end+1}=distStr; 
    elseif isfield(lmeInfo, 'distFinal')
        infoNames{end+1}='Distribution'; infoVals{end+1}=lmeInfo.distFinal; 
    end
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

% --- 3. Conditional ANOVA ---
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
            
            % Inject empty columns for alignment
            anovaTbl.Estimate = nan(height(anovaTbl), 1);
            anovaTbl.CI95 = repmat({''}, height(anovaTbl), 1);
            anovaTbl.SE = nan(height(anovaTbl), 1);
            
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

% --- 4. Fixed Effects / Coefficients (Type == "Coeff") ---
if ~isempty(lmeStats) && ismember('Type', lmeStats.Properties.VariableNames)
    idxCoef = lmeStats.Type == "Coeff";
    if any(idxCoef)
        coefTbl = lmeStats(idxCoef, :);
        coefTbl = removevars_safely(coefTbl, {'Index', 'Type', 'HVec'});
        add_tbl('FIXED EFFECTS (COEFFICIENTS)', coefTbl, 'Coeff');
    end
end

% --- 5. Post-Hoc (Type == "Simple" or "Marginal") ---
if ~isempty(lmeStats) && ismember('Type', lmeStats.Properties.VariableNames)
    idxPost = ismember(lmeStats.Type, ["Simple", "Marginal"]);
    if any(idxPost)
        postTbl = lmeStats(idxPost, :);
        postTbl = removevars_safely(postTbl, {'Index', 'Type', 'HVec'});
        add_tbl('POST-HOC EFFECTS', postTbl, 'PostHoc');
    end
end

% --- 6. Covariance & Error Parameters ---
if ~isempty(lmeMdl)
    try
        [~, ~, resStats] = covarianceParameters(lmeMdl);
        
        % Format Covariance
        covTbl = dataset2table(resStats{1});
        if ~isempty(covTbl)
            descStrs = repmat({''}, height(covTbl), 1);
            for r = 1:height(covTbl)
                grp = ''; if ismember('Group', covTbl.Properties.VariableNames), grp = char(covTbl.Group(r)); end
                n1 = ''; if ismember('Name1', covTbl.Properties.VariableNames), n1 = char(covTbl.Name1(r)); end
                n2 = ''; if ismember('Name2', covTbl.Properties.VariableNames), n2 = char(covTbl.Name2(r)); end
                
                if strcmp(grp, 'DummyVar') || isempty(grp), grp = 'Fixed'; end
                descStr = sprintf('Cov: %s', grp);
                if ~isempty(n1), descStr = [descStr, ' (', n1]; end
                if ~isempty(n2) && ~strcmp(n1, n2), descStr = [descStr, ' x ', n2]; end
                if ~isempty(n1), descStr = [descStr, ')']; end
                
                descStrs{r} = descStr;
            end
            covTbl = addvars(covTbl, descStrs, 'Before', 1, 'NewVariableNames', 'Description');
            covTbl = removevars_safely(covTbl, {'Group', 'Name1', 'Name2', 'Type'});
            if ismember('StandardError', covTbl.Properties.VariableNames)
                covTbl.Properties.VariableNames{strcmp(covTbl.Properties.VariableNames, 'StandardError')} = 'SE';
            end
            
            if ismember('Lower', covTbl.Properties.VariableNames) && ismember('Upper', covTbl.Properties.VariableNames)
                ci95List = cell(height(covTbl), 1);
                for r = 1:height(covTbl)
                    ci95List{r} = [covTbl.Lower(r), covTbl.Upper(r)];
                end
                covTbl = addvars(covTbl, ci95List, 'NewVariableNames', 'CI95');
                covTbl = removevars_safely(covTbl, {'Lower', 'Upper'});
            end
            
            add_tbl('COVARIANCE PARAMETERS', covTbl, 'Other');
        end
        
        % Format Error parameters
        errTbl = dataset2table(resStats{2});
        if ~isempty(errTbl)
            descStrs = repmat({''}, height(errTbl), 1);
            for r = 1:height(errTbl)
                n = ''; if ismember('Name', errTbl.Properties.VariableNames), n = char(errTbl.Name(r)); end
                if isempty(n), n = 'Residual'; end
                descStrs{r} = sprintf('Err: %s', n);
            end
            errTbl = addvars(errTbl, descStrs, 'Before', 1, 'NewVariableNames', 'Description');
            errTbl = removevars_safely(errTbl, {'Name', 'Group', 'Type'});
            if ismember('StandardError', errTbl.Properties.VariableNames)
                errTbl.Properties.VariableNames{strcmp(errTbl.Properties.VariableNames, 'StandardError')} = 'SE';
            end
            
            if ismember('Lower', errTbl.Properties.VariableNames) && ismember('Upper', errTbl.Properties.VariableNames)
                ci95List = cell(height(errTbl), 1);
                for r = 1:height(errTbl)
                    ci95List{r} = [errTbl.Lower(r), errTbl.Upper(r)];
                end
                errTbl = addvars(errTbl, ci95List, 'NewVariableNames', 'CI95');
                errTbl = removevars_safely(errTbl, {'Lower', 'Upper'});
            end
            
            add_tbl('ERROR PARAMETERS', errTbl, 'Other');
        end
        
    catch
    end
end

end % EOF

%% ========================================================================
%  HELPERS
%  ========================================================================

function tblOut = reorder_columns(tblIn)
% REORDER_COLUMNS Enforces a strict global column ordering.
globalOrder = {'Description', 'Estimate', 'CI95', 'SE', 'Statistic', 'DF', 'pVal'};

currentVars = tblIn.Properties.VariableNames;
orderedVars = {};

% 1. Extract columns present in the global priority list
for i = 1:length(globalOrder)
    col = globalOrder{i};
    if ismember(col, currentVars)
        orderedVars{end+1} = col;
    end
end

% 2. Extract remaining columns not in the priority list
remainingVars = setdiff(currentVars, orderedVars, 'stable');
finalVars = [orderedVars, remainingVars];

% 3. Reorder
tblOut = tblIn(:, finalVars);
end

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