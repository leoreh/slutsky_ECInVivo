function lme_save(varargin)

% Saves the figure handle in specified graphic formats and also saves LME
% analysis data (lme_tbl, lme_results, lme_cfg) to XLSX and MAT files.
%
% Assumes all inputs are correct and valid when used in the body of the function.
% lme_results is assumed to be a table.
%
% Graphic formats supported: jpg, fig, ai (Adobe Illustrator), svg.
% Data formats supported: xlsx, mat.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments using inputParser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'fh', gcf, @(x) isa(x, 'matlab.ui.Figure') || isempty(x));
addOptional(p, 'fpath', 'D:\\OneDrive - Tel-Aviv University\\PhD\\Slutsky\\Manuscripts\\Heim2025\\Graphs', @ischar);
addOptional(p, 'fname', 'figure', @ischar);
addOptional(p, 'frmt', {'ai', 'jpg', 'svg'}, @(x) iscell(x) || ischar(x)); % Default formats, inputParser handles cell/char validation implicitly to some extent
addOptional(p, 'lme_tbl', [], @(x) istable(x) || isempty(x));
addOptional(p, 'lme_results', [], @(x) istable(x) || isempty(x)); % Simplified: assume table or empty
addOptional(p, 'lme_cfg', [], @(x) isstruct(x) || isempty(x));

parse(p, varargin{:});
fh = p.Results.fh;
fpath = p.Results.fpath;
fname = p.Results.fname;
frmt = p.Results.frmt;
lme_tbl = p.Results.lme_tbl;
lme_results = p.Results.lme_results; 
lme_cfg = p.Results.lme_cfg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ensure path exists
mkdir(fpath); % Assumes fpath is valid and mkdir will succeed or path already exists

% simplify fname using frml2char (assumed to exist and work)
fname = lme_frml2char(fname);

% construct full file path base (without extension)
fullPathBase = fullfile(fpath, fname);

% assert figure format is cell
if ~iscell(frmt)
    frmt = {char(frmt)}; % Assumes frmt is char or string if not cell
end
nfrmts = length(frmt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save based on format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ifrmt = 1 : nfrmts
    current_format = lower(frmt{ifrmt});
    switch current_format
        case 'jpg'
            print(fh, fullPathBase, '-djpeg', '-r300');
        case 'fig'
            savefig(fh, [fullPathBase '.fig']);
        case 'ai'
            print(fh, fullPathBase, '-depsc', '-vector'); % Using .eps for AI compatibility
        case 'svg'
            print(fh, fullPathBase, '-dsvg', '-r300'); 
        case 'xlsx'
            xlsx_fullPath = [fullPathBase '.xlsx'];

            % Save lme_tbl to Sheet 1
            writetable(lme_tbl, xlsx_fullPath, 'Sheet', 'LME_Data');
            
            % Save lme_results to Sheet 2
            writetable(lme_results, xlsx_fullPath, 'Sheet', 'LME_Analysis'); 
     
        case 'mat'
            mat_fullPath = [fullPathBase '.mat'];
            varsToSave = struct();
            varsToSave.fh = fh;
            varsToSave.lme_tbl = lme_tbl;
            varsToSave.lme_results = lme_results;
            varsToSave.lme_cfg = lme_cfg;
            
            save(mat_fullPath, '-struct', 'varsToSave');
    end
end

end

% EOF 