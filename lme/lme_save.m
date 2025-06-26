function lme_save(varargin)

% Saves the figure handle in specified graphic formats and also saves LME
% analysis data (lmeData, lmeStats, lmeMdl) to XLSX and MAT files.
%
% Assumes all inputs are correct and valid when used in the body of the function.
% lmeStats is assumed to be a table.
%
% Graphic formats supported: jpg, fig, ai (Adobe Illustrator), svg.
% Data formats supported: xlsx, mat.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments using inputParser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'hFig', [], @(x) isa(x, 'matlab.ui.Figure') || isempty(x));
addOptional(p, 'savepath', 'D:\OneDrive - Tel-Aviv University\PhD\Slutsky\Manuscripts\Heim2025\Graphs', @ischar);
addOptional(p, 'fpath', [], @(x) ischar(x) || isempty(x));
addOptional(p, 'fname', [], @(x) ischar(x) || isempty(x));
addOptional(p, 'frmt', {'ai', 'jpg', 'svg'}, @(x) iscell(x) || ischar(x)); 
addOptional(p, 'lmeData', [], @(x) istable(x) || isempty(x));
addOptional(p, 'lmeStats', [], @(x) istable(x) || isempty(x));
addOptional(p, 'lmeMdl', [], @(x) isobject(x) || isempty(x));
addOptional(p, 'sheetNames', [], @(x) iscell(x) || isempty(x));

parse(p, varargin{:});
hFig = p.Results.hFig;
savepath = p.Results.savepath;
fpath = p.Results.fpath;
fname = p.Results.fname;
frmt = p.Results.frmt;
lmeData = p.Results.lmeData;
lmeStats = p.Results.lmeStats; 
lmeMdl = p.Results.lmeMdl;
sheetNames = p.Results.sheetNames;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set default sheet names if not provided
if isempty(sheetNames)
    sheetNames = {'LME_Data', 'LME_Model', 'LME_Stats'};
end

% simplify fname using frml2char 
if isempty(fname)
    fname = lme_frml2char(fname);
end

% create new folder based on response variable if fpath is empty
if isempty(fpath)
    resName = regexp(fname, '^[^~]+', 'match', 'once');
    resName = strtrim(resName);
    fpath = fullfile(savepath, resName);
end
mkdir(fpath);

% construct full file path base (without extension)
saveName = fullfile(fpath, fname);

% assert figure format is cell
if ~iscell(frmt)
    frmt = {char(frmt)}; % Assumes frmt is char or string if not cell
end
nfrmts = length(frmt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save based on format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ifrmt = 1 : nfrmts
    currentFormat = lower(frmt{ifrmt});
    switch currentFormat
        case 'jpg'
            print(hFig, saveName, '-djpeg', '-r300');
        case 'fig'
            savefig(hFig, [saveName '.fig']);
        case 'ai'
            print(hFig, saveName, '-depsc', '-vector'); % Using .eps for AI compatibility
        case 'svg'
            print(hFig, saveName, '-dsvg', '-r300'); 
        case 'xlsx'
            xlsName = [saveName '.xlsx'];
            
            % Save lmeData to Sheet 1
            if ~isempty(lmeData)
                writetable(lmeData, xlsName, 'Sheet', sheetNames{1});
            end
                       
            % Save lmeMdl to Sheet 2
            if ~isempty(lmeMdl) 
                % Get LME tables from the model
                lmeTbls = lme_mdl2tbls(lmeMdl);
                
                % Write formula in first cell of Sheet 2
                writematrix(lmeMdl.Formula.char, xlsName, 'Sheet', sheetNames{2}, 'Range', 'A1');
                
                % Get field names of lmeTbl
                tblFlds = fieldnames(lmeTbls);
                
                % Start after formula (add 2 rows for gap)
                currentRow = 3;
                
                % Loop through each table and write to Excel
                for iFld = 1:length(tblFlds)
                    % Get current table
                    currTbl = lmeTbls.(tblFlds{iFld});
                    
                    % Write table name as header
                    writematrix(tblFlds{iFld}, xlsName, 'Sheet', sheetNames{2}, 'Range', ['A' num2str(currentRow)]);
                    currentRow = currentRow + 1;
                    
                    % Write table data
                    writetable(currTbl, xlsName, 'Sheet', sheetNames{2}, 'Range', ['A' num2str(currentRow)]);
                    
                    % Update row counter (add table height + 2 for gap)
                    currentRow = currentRow + height(currTbl) + 2;
                end
            end
            
            % Save lmeStats to Sheet 3
            if ~isempty(lmeStats)
                writetable(lmeStats, xlsName, 'Sheet', sheetNames{3});
            end

        case 'mat'
            matName = [saveName '.mat'];
            varsToSave = struct();
            varsToSave.hFig = hFig;
            varsToSave.lmeData = lmeData;
            varsToSave.lmeStats = lmeStats;
            varsToSave.lmeMdl = lmeMdl;
            
            save(matName, '-struct', 'varsToSave');
    end
end

end

% EOF 