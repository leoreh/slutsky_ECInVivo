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
addOptional(p, 'fh', gcf, @(x) isa(x, 'matlab.ui.Figure') || isempty(x));
addOptional(p, 'basepath', 'D:\OneDrive - Tel-Aviv University\PhD\Slutsky\Manuscripts\Heim2025\Graphs', @ischar);
addOptional(p, 'fname', 'figure', @ischar);
addOptional(p, 'frmt', {'ai', 'jpg', 'svg'}, @(x) iscell(x) || ischar(x)); % Default formats, inputParser handles cell/char validation implicitly to some extent
addOptional(p, 'lmeData', [], @(x) istable(x) || isempty(x));
addOptional(p, 'lmeStats', [], @(x) istable(x) || isempty(x)); % Simplified: assume table or empty
addOptional(p, 'lmeMdl', [], @(x) isobject(x) || isempty(x));

parse(p, varargin{:});
fh = p.Results.fh;
basepath = p.Results.basepath;
fname = p.Results.fname;
frmt = p.Results.frmt;
lmeData = p.Results.lmeData;
lmeStats = p.Results.lmeStats; 
lmeMdl = p.Results.lmeMdl;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% simplify fname using frml2char (assumed to exist and work)
fname = lme_frml2char(fname);

% create new folder based on response variable
resName = regexp(fname, '^[^~]+', 'match', 'once');
resName = strtrim(resName);

fpath = fullfile(basepath, resName);
mkdir(fpath); % Assumes basepath is valid and mkdir will succeed or path already exists

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
            print(fh, saveName, '-djpeg', '-r300');
        case 'fig'
            savefig(fh, [saveName '.fig']);
        case 'ai'
            print(fh, saveName, '-depsc', '-vector'); % Using .eps for AI compatibility
        case 'svg'
            print(fh, saveName, '-dsvg', '-r300'); 
        case 'xlsx'
            xlsName = [saveName '.xlsx'];
            
            % Save lmeData to Sheet 1
            writetable(lmeData, xlsName, 'Sheet', 'LME_Data');
            
            % Save lmeStats to Sheet 2
            writetable(lmeStats, xlsName, 'Sheet', 'LME_Stats');
            
            % Get additional tables from the model if lmeMdl is provided
            if ~isempty(lmeMdl) 
                % Get LME tables from the model
                lmeTbls = lme_mdl2tbls(lmeMdl);
                
                % Get field names of lmeTbl
                tblFlds = fieldnames(lmeTbls);
                
                % Start after lmeStats (add 2 rows for gap)
                currentRow = height(lmeStats) + 4;
                
                % Loop through each table and write to Excel
                for iFld = 1:length(tblFlds)
                    % Get current table
                    currTbl = lmeTbls.(tblFlds{iFld});
                    
                    % Write table name as header
                    writematrix(tblFlds{iFld}, xlsName, 'Sheet', 'LME_Stats', 'Range', ['A' num2str(currentRow)]);
                    currentRow = currentRow + 1;
                    
                    % Write table data
                    writetable(currTbl, xlsName, 'Sheet', 'LME_Stats', 'Range', ['A' num2str(currentRow)]);
                    
                    % Update row counter (add table height + 2 for gap)
                    currentRow = currentRow + height(currTbl) + 3;
                end
            end

        case 'mat'
            matName = [saveName '.mat'];
            varsToSave = struct();
            varsToSave.fh = fh;
            varsToSave.lmeData = lmeData;
            varsToSave.lmeStats = lmeStats;
            varsToSave.lmeMdl = lmeMdl;
            
            save(matName, '-struct', 'varsToSave');
    end
end

end

% EOF 