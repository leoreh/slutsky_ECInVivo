function mcu_xlsFormat(xlsName, varargin)
% MCU_XLSFORMAT Formats the supplementary Excel file via ActiveX server.
%
%   MCU_XLSFORMAT(XLSNAME) opens the specified Excel file and applies
%   the following formatting across all sheets: left alignment, asserting
%   hyperlink formulas, highlighting rows with 'MODEL INFORMATION', and
%   applying bold font to cells containing text in all capital letters.
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
        
        for iRow = 1:numRows
            for iCol = 1:numCols
                cellObj = get(usedRange, 'Item', iRow, iCol);
                cellVal = cellObj.Value;
                
                % Assert hyperlinks by forcing Excel to re-evaluate the local formula
                cellFormula = cellObj.Formula;
                if ischar(cellFormula) && startsWith(cellFormula, '=HYPERLINK')
                    cellObj.FormulaLocal = cellFormula;
                end
                
                if ischar(cellVal)
                    % Highlight entire row in light gray
                    if contains(cellVal, 'MODEL INFORMATION')
                        rowObj = cellObj.EntireRow;
                        colorIndexLightGray = 15;
                        rowObj.Interior.ColorIndex = colorIndexLightGray;
                    end
                    
                    % Apply bold to cells that are entirely uppercase
                    % if isequal(upper(cellVal), cellVal) && all(isletter(cellVal))
                    %     rowObj = cellObj.EntireRow;
                    %     rowObj.Font.Bold = true;
                    % end
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