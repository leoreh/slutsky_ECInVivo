function [outMat, rowLabels, colLabels] = tbl2prism(tbl, varargin)

% TBL2PRISM converts a MATLAB table into a matrix format suitable for 
% GraphPad Prism. It handles two modes:
% 1. One-Way (Column Table): Data arranged vertically in columns.
% 2. Two-Way (Grouped Table): Data arranged horizontally in subcolumns,
%    padded to the global maximum replicate count.
%
% INPUT (Required):
%   tbl         - Input table (tidy format).
%   yVar        - String name of the data variable (measure).
%
% INPUT (Name-Value):
%   grpVar      - String name of the column grouping variable (e.g., 'Genotype').
%                 In 2-Way mode, these form the main column blocks.
%   rowVar      - String name of the row variable (e.g., 'UnitType').
%                 If provided, triggers 2-Way (Grouped) mode.
%   flgClip     - Boolean. Copy to clipboard {true}
%
% OUTPUT:
%   outMat      - The formatted numeric matrix (NaN padded).
%
% EXAMPLE:
%   % 1-Way (Column Table):
%   tbl2prism(myTbl, 'yVar', 'fr', 'grpVar', 'Group');
%   
%   % 2-Way (Grouped Table - Horizontal with subcolumns):
%   tbl2prism(myTbl, 'yVar', 'fr', 'grpVar', 'Group', 'rowVar', 'UnitType');

% 10 Dec 2025

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addParameter(p, 'yVar', '', @ischar);
addParameter(p, 'grpVar', 'Group', @ischar);
addParameter(p, 'rowVar', '', @ischar); % Defines rows in 2-way
addParameter(p, 'flgClip', true, @islogical);
addParameter(p, 'flgSort', true, @islogical);

parse(p, tbl, varargin{:});
tbl     = p.Results.tbl;
yVar    = p.Results.yVar;
grpVar  = p.Results.grpVar;
rowVar  = p.Results.rowVar;
flgClip = p.Results.flgClip;
flgSort = p.Results.flgSort;

%% ========================================================================
%  PREPARE LABELS
%  ========================================================================

% Get Group Columns (e.g., WT, KO)
uGrps = unique(tbl.(grpVar));
if flgSort && exist('natsort', 'file')
    uGrps = string(natsort(cellstr(string(uGrps))));
else
    uGrps = string(uGrps);
end
nGrps = length(uGrps);

% Get Row Labels (e.g., RS, FS) if applicable
if ~isempty(rowVar)
    uRows = unique(tbl.(rowVar));
    if flgSort && exist('natsort', 'file')
        uRows = string(natsort(cellstr(string(uRows))));
    else
        uRows = string(uRows);
    end
    nRows = length(uRows);
    isTwoWay = true;
else
    uRows = "Data"; % Dummy label for 1-way
    nRows = 1;
    isTwoWay = false;
end

%% ========================================================================
%  1-WAY MODE (COLUMN TABLE - VERTICAL)
%  ========================================================================
if ~isTwoWay
    fprintf('Mode: 1-Way (Column Table). Vertical columns.\n');
    
    % Collect data
    dataCell = cell(1, nGrps);
    for i = 1:nGrps
        idx = tbl.(grpVar) == uGrps(i);
        dataCell{i} = tbl.(yVar)(idx);
    end
    
    % Pad to max length
    nMax = max(cellfun(@length, dataCell));
    outMat = nan(nMax, nGrps);
    for i = 1:nGrps
        outMat(1:length(dataCell{i}), i) = dataCell{i};
    end
    
    rowLabels = [];
    colLabels = uGrps;
    
    % Clipboard String
    strOut = strjoin(colLabels, '\t') + newline;
    for r = 1:nMax
        rowStr = "";
        for c = 1:nGrps
            val = outMat(r,c);
            if ~isnan(val), rowStr = rowStr + sprintf('%.6g', val); end
            if c < nGrps, rowStr = rowStr + sprintf('\t'); end
        end
        strOut = strOut + rowStr + newline;
    end

%% ========================================================================
%  2-WAY MODE (GROUPED TABLE - HORIZONTAL BLOCKS)
%  ========================================================================
else
    fprintf('Mode: 2-Way (Grouped Table). Horizontal subcolumns.\n');
    
    % 1. Determine Global Max Replicates (Subcolumns)
    % We need the maximum N across ALL conditions to align blocks
    maxRep = 0;
    for iR = 1:nRows
        for iG = 1:nGrps
            idx = tbl.(rowVar) == uRows(iR) & tbl.(grpVar) == uGrps(iG);
            maxRep = max(maxRep, sum(idx));
        end
    end
    fprintf('Global Max Replicates: %d (Padding all groups to this width)\n', maxRep);

    % 2. Initialize Matrix
    % Rows = uRows, Cols = uGrps * maxRep
    outMat = nan(nRows, nGrps * maxRep);
    
    % 3. Fill Matrix
    for iG = 1:nGrps
        % Calculate column offset for this group block
        colStart = (iG-1) * maxRep + 1;
        
        for iR = 1:nRows
            % Extract Data
            idx = tbl.(rowVar) == uRows(iR) & tbl.(grpVar) == uGrps(iG);
            dat = tbl.(yVar)(idx);
            
            % Fill horizontally in the specific block
            if ~isempty(dat)
                % Ensure vector is row
                dat = reshape(dat, 1, []);
                % Place in matrix
                outMat(iR, colStart : colStart + length(dat) - 1) = dat;
            end
        end
    end
    
    rowLabels = uRows;
    colLabels = uGrps;
    
    % 4. Clipboard String (With Headers)
    % Header: Empty | Group1 | (tabs...) | Group2 | ...
    strOut = "" + sprintf('\t'); % Skip row label col
    for iG = 1:nGrps
        strOut = strOut + uGrps(iG);
        % Add tabs for the rest of the subcolumns in this group
        if iG < nGrps
            strOut = strOut + repmat(sprintf('\t'), 1, maxRep);
        end
    end
    strOut = strOut + newline;
    
    % Data Rows
    for r = 1:nRows
        strOut = strOut + uRows(r) + sprintf('\t'); % Row Label
        for c = 1:(nGrps * maxRep)
            val = outMat(r,c);
            if ~isnan(val)
                strOut = strOut + sprintf('%.6g', val);
            end
            % Add tab unless it's the very last element
            if c < (nGrps * maxRep)
                strOut = strOut + sprintf('\t');
            end
        end
        strOut = strOut + newline;
    end
end

%% ========================================================================
%  COPY TO CLIPBOARD
%  ========================================================================

if flgClip
    clipboard('copy', char(strOut));
    fprintf('Copied to clipboard.\n');
end

end

% EOF