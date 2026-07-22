function [outMat, hdr] = tbl2prismSum(tbl, varargin)

% TBL2PRISMSUM Summarises a table variable as per-group Mean / SD / N for
% GraphPad Prism and copies the block to the clipboard. Companion to
% tbl2prism, which instead lays out the raw replicates.
%
% Two shapes are handled automatically from the class of yVar:
%   Scalar  (one value per row)          -> COLUMN block [3 x nGrp]:
%           rows {Mean; SD; N}, one column per group. Paste into a Prism
%           Column table set to 'Enter Mean, SD, N'.
%   Vector  (matrix, or cell of vectors) -> GROUPED block [nX x (1+3*nGrp)]:
%           first column is X, then a {Mean, SD, N} sub-block per group.
%           Paste into a Prism Grouped table set to 'Enter Mean, SD, N'.
%
% INPUT (Required):
%   tbl         (table) Tidy data table.
%   yVar        (char)  Data variable to summarise.
%
% INPUT (Name-Value):
%   grpVar      (char)  Grouping variable. {'genotype'}
%   xVec        (vec)   X axis for the vector shape. Defaults to 1:nX.
%   flgClip     (log)   Copy the block to the clipboard. {true}
%
% OUTPUT:
%   outMat      (numeric) The Mean/SD/N block (see shapes above).
%   hdr         (struct)  .grps and .stats label the block.
%
% EXAMPLE:
%   tbl2prismSum(tblRipp, 'yVar', 'freq');                  % bars
%   tbl2prismSum(tblMaps, 'yVar', 't_lfp', 'xVec', xVec);   % waveform
%
% See also: TBL2PRISM

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addParameter(p, 'yVar', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'grpVar', 'genotype', @(x) ischar(x) || isstring(x));
addParameter(p, 'xVec', [], @isnumeric);
addParameter(p, 'flgClip', true, @islogical);

parse(p, tbl, varargin{:});
yVar    = char(p.Results.yVar);
grpVar  = char(p.Results.grpVar);
xVec    = p.Results.xVec;
flgClip = p.Results.flgClip;

% Group list. A categorical keeps its declared level order (Control, MCU-KO,
% CAG-MCU-KO), so drop only the absent ones; anything else is natsorted.
grp = tbl.(grpVar);
if ~iscategorical(grp), grp = categorical(grp); end
grp = removecats(grp);
uGrps = categories(grp);
if ~iscategorical(tbl.(grpVar)) && exist('natsort', 'file')
    uGrps = natsort(uGrps);
end
nGrp = numel(uGrps);

% Scalar (one value per row) vs vector (matrix / cell of vectors) shape.
raw = tbl.(yVar);
isVec = iscell(raw) || size(raw, 2) > 1;

%% ========================================================================
%  SCALAR - COLUMN BLOCK  [Mean; SD; N] x nGrp
%  ========================================================================

if ~isVec
    outMat = nan(3, nGrp);
    for iGrp = 1:nGrp
        vals = raw(grp == uGrps{iGrp});
        vals = vals(~isnan(vals));
        outMat(:, iGrp) = [mean(vals); std(vals); numel(vals)];
    end
    hdr.stats = {'Mean'; 'SD'; 'N'};
    hdr.grps  = uGrps;

    % Clipboard: header row of group names, then one row per statistic.
    lines = {strjoin([{''}, uGrps(:)'], char(9))};
    for iStat = 1:3
        vals = arrayfun(@(c) num2str(outMat(iStat, c), '%.6g'), ...
            1:nGrp, 'uni', false);
        lines{end + 1} = strjoin([hdr.stats(iStat), vals], char(9)); %#ok<AGROW>
    end

%% ========================================================================
%  VECTOR - GROUPED BLOCK  X | {Mean SD N} per group
%  ========================================================================

else
    if iscell(raw)
        nX = numel(raw{find(~cellfun(@isempty, raw), 1)});
    else
        nX = size(raw, 2);
    end
    if isempty(xVec), xVec = (1:nX)'; end

    blocks = nan(nX, 3 * nGrp);
    for iGrp = 1:nGrp
        idx = grp == uGrps{iGrp};
        if iscell(raw)
            mat = cell2mat(cellfun(@(x) x(:)', raw(idx), 'uni', false));
        else
            mat = raw(idx, :);
        end
        col = (iGrp - 1) * 3 + (1:3);
        blocks(:, col) = [mean(mat, 1, 'omitnan')', ...
            std(mat, [], 1, 'omitnan')', sum(~isnan(mat), 1)'];
    end
    outMat = [xVec(:), blocks];
    hdr.stats = {'Mean', 'SD', 'N'};
    hdr.grps  = uGrps;

    % Clipboard: two header rows (group name over its block, then Mean/SD/N).
    h1 = {'X'};
    h2 = {''};
    for iGrp = 1:nGrp
        h1 = [h1, uGrps{iGrp}, '', '']; %#ok<AGROW>
        h2 = [h2, hdr.stats{:}];        %#ok<AGROW>
    end
    lines = {strjoin(h1, char(9)), strjoin(h2, char(9))};
    for iX = 1:nX
        vals = arrayfun(@(c) num2str(outMat(iX, c), '%.6g'), ...
            1:size(outMat, 2), 'uni', false);
        lines{end + 1} = strjoin(vals, char(9)); %#ok<AGROW>
    end
end

%% ========================================================================
%  CLIPBOARD
%  ========================================================================

if flgClip
    clipboard('copy', strjoin(lines, newline));
    fprintf('tbl2prismSum: %s by %s copied (%d groups).\n', yVar, grpVar, nGrp);
end

end     % EOF
