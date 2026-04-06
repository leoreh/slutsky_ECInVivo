function mcu_xlsFormat(xlsName, varargin)
% MCU_XLSFORMAT Formats the supplementary Excel file via ActiveX server.
%
%   MCU_XLSFORMAT(XLSNAME) opens the specified Excel file and applies
%   the following formatting across all sheets: left alignment, asserting
%   hyperlink formulas, and highlighting formula rows (containing ' ~ ')
%   with bold font, gray background, and merged centered columns A-H.
%
%   Optional Parameters:
%       pathName - Directory path (default: current directory)
%       verbose  - Print progress (default: true)

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'xlsName', @(x) ischar(x) || isstring(x));

defaultPath = pwd;

addParameter(p, 'pathName', defaultPath, @(x) ischar(x) || isstring(x));
addParameter(p, 'verbose', true, @islogical);

parse(p, xlsName, varargin{:});

fullXlsPath = fullfile(p.Results.pathName, p.Results.xlsName);

%% ========================================================================
%  FORMAT VIA COM SERVER
%  ========================================================================

try
    excelApp = actxserver('Excel.Application');
    excelApp.Visible = false;
    excelApp.DisplayAlerts = false;
    
    workbook = excelApp.Workbooks.Open(fullXlsPath);
    numSheets = workbook.Sheets.Count;
    
    for iSheet = 1:numSheets
        sheet = get(workbook.Sheets, 'Item', iSheet);
        if p.Results.verbose
            fprintf('[MCU_XLSFORMAT] Formatting sheet: %s\n', sheet.Name);
        end
        
        usedRange = sheet.UsedRange;
        
        % Align everything to the left
        xlHAlignLeft = -4131;
        usedRange.HorizontalAlignment = xlHAlignLeft;
        
        numRows = usedRange.Rows.Count;
        numCols = usedRange.Columns.Count;
        
        xlHAlignCenter = -4108;
        colorIndexLightGray = 15;

        for iRow = 1:numRows

            % Check first cell for formula rows (contain ' ~ ')
            cellA = get(usedRange, 'Item', iRow, 1);
            cellAVal = cellA.Value;
            if ischar(cellAVal) && contains(cellAVal, ' ~ ')
                absRow = cellA.Row;
                mergeRange = sheet.Range(sprintf('A%d:H%d', absRow, absRow));
                mergeRange.Merge;
                mergeRange.HorizontalAlignment = xlHAlignCenter;
                mergeRange.Font.Bold = true;
                mergeRange.Interior.ColorIndex = colorIndexLightGray;
                continue;
            end

            % Cell-by-cell processing for non-formula rows
            for iCol = 1:numCols
                cellObj = get(usedRange, 'Item', iRow, iCol);

                % Assert hyperlinks by forcing Excel to re-evaluate the local formula
                cellFormula = cellObj.Formula;
                if ischar(cellFormula) && startsWith(cellFormula, '=HYPERLINK')
                    cellObj.FormulaLocal = cellFormula;
                end
            end
        end
    end
    
    workbook.Save();
    workbook.Close();
    excelApp.Quit();
    delete(excelApp);
    
    if p.Results.verbose
        fprintf('[MCU_XLSFORMAT] Formatting completed successfully.\n');
    end
    
catch ME
    if exist('workbook', 'var')
        workbook.Close(false);
    end
    if exist('excelApp', 'var')
        excelApp.Quit();
        delete(excelApp);
    end
    rethrow(ME);
end
end