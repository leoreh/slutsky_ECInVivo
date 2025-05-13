function lmeTbl = lme_mdl2tbl(lmeMdl)
% ORGANIZE_LME_FOR_EXCEL Extracts key information from an LME model
% into a single wide table suitable for Excel export.
%
% INPUT:
%   lmeMdl - A LinearMixedModel object (from fitlme).
%
% OUTPUT:
%   outputTable - A MATLAB table with combined LME model information.

% Extract Model-Level Information 
formulaStr = lmeMdl.Formula.char;
fitMet = lmeMdl.FitMethod;
numObs = lmeMdl.NumObservations;

% Extract model fit params
fitCrit = lmeMdl.ModelCriterion;
fitStats = dataset2table(fitCrit);

% Extract Random Effects and Covariance Parameters
[~, ~, randStats] = randomEffects(lmeMdl, 'DFMethod', 'satterthwaite');
randStats = dataset2table(randStats);
fixedStats = dataset2table(lmeMdl.Coefficients);

[~, ~, resStats] = covarianceParameters(lmeMdl);
covStats = dataset2table(resStats{1});
errStats = dataset2table(resStats{2});

% Define all unique column names including a source indicator column
allVarNames = {'SourceTable', 'Group', 'Level', 'Name', 'Estimate', 'SE', 'SEPred', 'tStat', 'DF', 'pValue', 'Lower', 'Upper', 'Name1', 'Name2', 'Type'};
% Define which of these columns are expected to be numeric for default filling
numericCols = {'Estimate', 'SE', 'SEPred', 'tStat', 'DF', 'pValue', 'Lower', 'Upper'};

processedTables = {}; % Cell array to hold standardized tables before final concatenation

% List of the original statistics tables and their corresponding source names
tablesToProcess = {randStats, fixedStats, covStats, errStats};
sourceNames = {'RandomEffects', 'FixedEffects', 'CovarianceParameters', 'ErrorStats'};

for k_table = 1:length(tablesToProcess)
    currentStatsTable = tablesToProcess{k_table};
    sourceName = sourceNames{k_table};
    
    % Skip if the current table is empty or has no rows
    if isempty(currentStatsTable) || height(currentStatsTable) == 0
        continue;
    end
    numRows = height(currentStatsTable);

    % Create a temporary cell array to build the data for the current table block
    % This allows for mixed data types before converting to a table
    tempTableDataAsCell = cell(numRows, length(allVarNames));
    
    % Find the index of 'SourceTable' column for efficient assignment
    sourceTableColIdx = find(strcmp(allVarNames, 'SourceTable'), 1);

    for r_idx = 1:numRows % Iterate over each row of the current input stats table
        % Assign the source name to the 'SourceTable' column
        tempTableDataAsCell{r_idx, sourceTableColIdx} = sourceName;

        % Populate other columns based on allVarNames
        for c_idx = 1:length(allVarNames)
            colName = allVarNames{c_idx};
            
            if strcmp(colName, 'SourceTable')
                continue; % This column was already handled
            end

            if ismember(colName, currentStatsTable.Properties.VariableNames)
                % Column exists in the source table, copy its data for the current row
                val = currentStatsTable.(colName)(r_idx);
                if iscell(val) % If the data in the source table is a cell (e.g., a cell containing a string)
                    tempTableDataAsCell{r_idx, c_idx} = val{1}; % Extract the content of the cell
                else
                    tempTableDataAsCell{r_idx, c_idx} = val; % Assign the value directly (e.g. numeric)
                end
            else
                % Column does not exist in the source table, fill with a default value
                if any(strcmp(colName, numericCols))
                    tempTableDataAsCell{r_idx, c_idx} = NaN; % Default for numeric columns
                else 
                    tempTableDataAsCell{r_idx, c_idx} = {''}; % Default for text-like columns (cell containing an empty char)
                end
            end
        end
    end
    % Convert the populated cell array for this block of rows into a table
    processedTables{end+1} = cell2table(tempTableDataAsCell, 'VariableNames', allVarNames);
end

% Combine all processed tables.
% vertcat can robustly combine tables that have the same column names and types.
if ~isempty(processedTables)
    lmeTbl = vertcat(processedTables{:});
else
    % If all input tables were empty, create an empty lmeTbl with the defined structure
    varTypes = cell(1, length(allVarNames)); % Determine variable types for the empty table
    for v_idx = 1:length(allVarNames)
        colName = allVarNames{v_idx};
        if any(strcmp(colName, numericCols))
            varTypes{v_idx} = 'double'; % Standard type for numeric data
        else 
            % Covers 'SourceTable' and other string-like columns (e.g., 'Group', 'Name')
            varTypes{v_idx} = 'cell'; % Cell type for strings or mixed content
        end
    end
    lmeTbl = table('Size', [0, length(allVarNames)], 'VariableTypes', varTypes, 'VariableNames', allVarNames);
end
% The variable 'lmeTbl' is now populated (or correctly empty) and will be returned by the function.

end