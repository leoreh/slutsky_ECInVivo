function lme_save(sheetName, lmeMdl, lmeStats, lmeInfo, varargin)
% LME_SAVE Exports LME analysis results to an organized Excel spreadsheet.
%
%   LME_SAVE(SHEETNAME, LMEMDL, LMESTATS, LMEINFO) 
%   saves the results of the LME analysis pipeline into a single, 
%   well-organized worksheet named SHEETNAME.
%
%   Optional Parameters:
%       pathName - Directory path (default: 'D:\OneDrive - Tel-Aviv University\PhD\Slutsky\Manuscripts\MCU\Results')
%       xlsName  - Excel file name (default: 'mcu_suppTbl.xlsx')
%       verbose  - Print progress (default: true)

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'sheetName', @(x) ischar(x) || isstring(x));
addRequired(p, 'lmeMdl', @(x) isempty(x) || ...
    isa(x, 'classreg.regr.LinearLikeMixedModel') || isobject(x));
addRequired(p, 'lmeStats', @(x) isempty(x) || istable(x));
addRequired(p, 'lmeInfo', @(x) isempty(x) || isstruct(x));

% Default path and name as requested
defaultPath = 'D:\OneDrive - Tel-Aviv University\PhD\Slutsky\Manuscripts\MCU\Results';
defaultXls = 'mcu_suppTbl.xlsx';

addParameter(p, 'pathName', defaultPath, @(x) ischar(x) || isstring(x));
addParameter(p, 'xlsName', defaultXls, @(x) ischar(x) || isstring(x));
addParameter(p, 'verbose', true, @islogical);

parse(p, sheetName, lmeMdl, lmeStats, lmeInfo, varargin{:});

sheetName = char(p.Results.sheetName);
pathName  = char(p.Results.pathName);
xlsName   = char(p.Results.xlsName);
verbose   = p.Results.verbose;

if isempty(pathName), pathName = pwd; end
if ~endsWith(xlsName, '.xlsx', 'IgnoreCase', true)
    xlsName = [char(xlsName), '.xlsx'];
end

% Construct full file path
fullXlsPath = fullfile(pathName, xlsName);

if verbose
    fprintf('[LME_SAVE] Exporting results to %s (Sheet: %s)...\n', fullXlsPath, sheetName);
end

% Ensure directory exists
if ~isfolder(pathName)
    mkdir(pathName);
end

%% ========================================================================
%  PREPARE SPREADSHEET CONTENT
%  ========================================================================

exportData = {};

% --- Helper Function to Append Formatted Tables ---
    function appendTable(tbl, titleStr, tableType)
        if nargin < 3, tableType = 'Other'; end
        if isempty(tbl) || height(tbl) == 0
            return;
        end
        
        % Add Title
        exportData{end+1} = {titleStr};
        
        % Format columns specifically based on table type
        tbl = format_table_for_pub(tbl, tableType);
        
        % Add Headers
        exportData{end+1} = tbl.Properties.VariableNames;
        
        % Add Data
        tblCell = table2cell(tbl);
        
        % Convert datasets, arrays, or strings for cleaner writing
        for r = 1:size(tblCell, 1)
            for c = 1:size(tblCell, 2)
                val = tblCell{r,c};
                if isnumeric(val) && (length(val) > 1)
                    tblCell{r,c} = char(join(string(val), ', ')); % cleaner array representation
                elseif isnumeric(val) && isnan(val)
                    tblCell{r,c} = '';
                elseif islogical(val)
                    if val, tblCell{r,c} = 'true'; else, tblCell{r,c} = 'false'; end
                elseif isstring(val) && ismissing(val)
                    tblCell{r,c} = '';
                elseif iscategorical(val)
                    tblCell{r,c} = char(val);
                elseif iscell(val) && numel(val) == 1
                    tblCell{r,c} = val{1};
                end
            end
        end
        
        % Append into exportData
        for r = 1:size(tblCell, 1)
            exportData{end+1} = tblCell(r, :);
        end
        exportData{end+1} = {''}; % Empty row separator
    end

% --- 1. Model Information ---
exportData{end+1} = {'MODEL INFORMATION'};
if ~isempty(lmeInfo)
    if isfield(lmeInfo, 'frml'), exportData{end+1} = {'Formula', lmeInfo.frml}; end
    if isfield(lmeInfo, 'distFinal'), exportData{end+1} = {'Distribution', lmeInfo.distFinal}; end
    if isfield(lmeInfo, 'fitMethod'), exportData{end+1} = {'Fit Method', lmeInfo.fitMethod}; end
    if isfield(lmeInfo, 'dfMethod'), exportData{end+1} = {'DF Method', lmeInfo.dfMethod}; end
    if isfield(lmeInfo, 'aic'), exportData{end+1} = {'AIC', mat2str(round(lmeInfo.aic, 2))}; end
end

if ~isempty(lmeMdl) && isprop(lmeMdl, 'NumObservations')
    exportData{end+1} = {'Observations', mat2str(lmeMdl.NumObservations)};
end
exportData{end+1} = {''}; % Empty row separator

% --- 2. Continuous Predictors Table ---
if ~isempty(lmeInfo) && isfield(lmeInfo, 'transParams') && isfield(lmeInfo.transParams, 'varsInc') && ~isempty(lmeInfo.transParams.varsInc)
    varsInc = lmeInfo.transParams.varsInc(:);
    
    numPredictors = length(varsInc);
    
    if numPredictors > 0
        predNames = cell(numPredictors, 1);
        isStnd = false(numPredictors, 1);
        transMethod = cell(numPredictors, 1);
        
        for i = 1:numPredictors
            vName = varsInc{i};
            predNames{i} = vName;
            
            % Check standardization status
            if isfield(lmeInfo, 'standardized')
                isStnd(i) = lmeInfo.standardized;
            else
                isStnd(i) = false; % fallback
            end
            
            % Check Transformation status inside varsTrans struct
            if isfield(lmeInfo.transParams, 'varsTrans') && isfield(lmeInfo.transParams.varsTrans, vName)
                vInfo = lmeInfo.transParams.varsTrans.(vName);
                if isfield(vInfo, 'funType')
                    t = lower(vInfo.funType);
                    if strcmp(t, 'log') && isfield(vInfo, 'logBase')
                        if strcmpi(vInfo.logBase, 'e')
                            transMethod{i} = 'Natural Log';
                        elseif strcmpi(vInfo.logBase, 'logit')
                            transMethod{i} = 'Logit';
                        else
                            transMethod{i} = sprintf('Log%s', mat2str(vInfo.logBase));
                        end
                    elseif ~isempty(t)
                        transMethod{i} = t; % Capitalize appropriately if needed
                    else
                        transMethod{i} = 'None';
                    end
                else
                    transMethod{i} = 'None';
                end
            else
                transMethod{i} = 'None';
            end
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
            appendTable(predTbl, 'CONTINUOUS PREDICTORS', 'Other');
        end
    end
end


% --- 3. Fixed Effects / Coefficients (Type == "Coeff") ---
if ~isempty(lmeStats) && ismember('Type', lmeStats.Properties.VariableNames)
    idxCoef = lmeStats.Type == "Coeff";
    if any(idxCoef)
        coefTbl = lmeStats(idxCoef, :);
        coefTbl = removevars_safely(coefTbl, {'Index', 'Type', 'HVec'});
        appendTable(coefTbl, 'FIXED EFFECTS (COEFFICIENTS)', 'Coeff');
    end
end


% --- 4. Post-Hoc (Type == "Simple" or "Marginal") ---
if ~isempty(lmeStats) && ismember('Type', lmeStats.Properties.VariableNames)
    idxPost = ismember(lmeStats.Type, ["Simple", "Marginal"]);
    if any(idxPost)
        postTbl = lmeStats(idxPost, :);
        % Explicitly remove Type here so it aligns with fixed effects columns
        postTbl = removevars_safely(postTbl, {'Index', 'Type', 'HVec'});
        appendTable(postTbl, 'POST-HOC EFFECTS', 'PostHoc');
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
                        if numLvls > 2
                            includeAnova = true;
                            break;
                        end
                    end
                end
            end
            if includeAnova, break; end 
        end
        
        if includeAnova
            anovaTbl = lmeStats(idxAnova, :);
            anovaTbl = removevars_safely(anovaTbl, {'Index', 'Type', 'Estimate', 'SE', 'CI95', 'HVec'});
            appendTable(anovaTbl, 'ANOVA (OMNIBUS INTERACTION)', 'ANOVA');
        end
        
        % Other Unclassified Effects (not ANOVA, Coef, Simple, Marginal)
        idxOther = ~(lmeStats.Type == "ANOVA" | lmeStats.Type == "Coeff" | ismember(lmeStats.Type, ["Simple", "Marginal"]));
        if any(idxOther)
            otherTbl = lmeStats(idxOther, :);
            otherTbl = removevars_safely(otherTbl, {'Index', 'HVec'});
            appendTable(otherTbl, 'OTHER LME EFFECTS', 'Other');
        end
        
    elseif ~ismember('Type', lmeStats.Properties.VariableNames)
        % Fallback if 'Type' isn't properly formatted internally for some reason
        appendTable(lmeStats, 'ALL LME STATISTICS', 'Other');
    end
end


% --- 6. Covariance & Error Parameters ---
if ~isempty(lmeMdl)
    try
        [~, ~, resStats] = covarianceParameters(lmeMdl);
        appendTable(dataset2table(resStats{1}), 'COVARIANCE PARAMETERS', 'Other');
        appendTable(dataset2table(resStats{2}), 'ERROR PARAMETERS', 'Other');
    catch
        % Silent catch if model type doesn't support covarianceParameters
    end
end

% --- 7. Append Logic for Existing Sheets ---
existingData = {};
if isfile(fullXlsPath)
    try
        % Safely get sheet names 
        sheets = sheetnames(fullXlsPath);
        if ismember(string(sheetName), sheets)
            % Read existing sheet
            existingData = readcell(fullXlsPath, 'Sheet', sheetName, 'DateType', 'text');
            
            % Sanitize existingData to replace 'missing' with empty strings
            % (prevents writecell from writing #N/A errors)
            if ~isempty(existingData)
                for r = 1:size(existingData, 1)
                    for c = 1:size(existingData, 2)
                        if ismissing(existingData{r,c})
                            existingData{r,c} = '';
                        end
                    end
                end
            end
        end
    catch
        % File might be locked or unreadable
    end
end

%% ========================================================================
%  WRITE TO EXCEL
%  ========================================================================

% Create a uniform rectangular cell array for the newly generated data
maxCols = max(cellfun(@length, exportData));
numRows = length(exportData);

newDataBlock = cell(numRows, maxCols);
newDataBlock(:) = {''};

for iRow = 1:numRows
    rowLen = length(exportData{iRow});
    if rowLen > 0
        newDataBlock(iRow, 1:rowLen) = exportData{iRow};
    end
end

% Combine existing and new data if extending
if ~isempty(existingData)
    % Add a 2-row gap
    gapCols = max(size(existingData, 2), size(newDataBlock, 2));
    gapBlock = cell(2, gapCols);
    gapBlock(:) = {''};
    
    % Align columns of existing data
    if size(existingData, 2) < gapCols
        existingData{1, gapCols} = ''; 
    end
    
    % Align columns of new data
    if size(newDataBlock, 2) < gapCols
        newDataBlock{1, gapCols} = ''; 
    end
    
    finalBuffer = [existingData; gapBlock; newDataBlock];
else
    finalBuffer = newDataBlock;
end

% Write via writecell
% By providing the fully concatenated array, we explicitly overwrite the sheet 
% with our properly preserved and appended data, entirely bypassing writecell bugs.
writecell(finalBuffer, fullXlsPath, 'Sheet', sheetName);

if verbose
    fprintf('[LME_SAVE] Successfully saved Excel file to Sheet: %s\n', sheetName);
end

end % EOF

%% ========================================================================
%  HELPERS
%  ========================================================================

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
        % Remove Adj P-value
        if ismember('Adj P-value', tbl.Properties.VariableNames)
            tbl = removevars_safely(tbl, {'Adj P-value'});
        end
        % Statistic is an F-statistic for ANOVA
        if ismember('Statistic', tbl.Properties.VariableNames)
            tbl.Properties.VariableNames{strcmp(tbl.Properties.VariableNames, 'Statistic')} = 'F-statistic';
        end
    case 'Coeff'
        % Remove Adj P-value
        if ismember('Adj P-value', tbl.Properties.VariableNames)
            tbl = removevars_safely(tbl, {'Adj P-value'});
        end
        % Statistic is a t-statistic
        if ismember('Statistic', tbl.Properties.VariableNames)
            tbl.Properties.VariableNames{strcmp(tbl.Properties.VariableNames, 'Statistic')} = 't-statistic';
        end
    case 'PostHoc'
        % Remove original P-value, keep only Adj P-value but named P-value (Adj.)
        if ismember('P-value', tbl.Properties.VariableNames)
            tbl = removevars_safely(tbl, {'P-value'});
        end
        if ismember('Adj P-value', tbl.Properties.VariableNames)
            tbl.Properties.VariableNames{strcmp(tbl.Properties.VariableNames, 'Adj P-value')} = 'P-value (Adj.)';
        end
        % Statistic is a t-statistic
        if ismember('Statistic', tbl.Properties.VariableNames)
            tbl.Properties.VariableNames{strcmp(tbl.Properties.VariableNames, 'Statistic')} = 't-statistic';
        end
    case 'Other'
        % Default handling for t-statistic assumption unless explicitly F
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