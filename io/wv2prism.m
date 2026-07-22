function blocks = wv2prism(tbl, xVec, varargin)
% WV2PRISM Waveform table -> GraphPad Prism XY blocks (Mean, SD, N per column).
%
%   blocks = WV2PRISM(tbl, xVec, varargin)
%
%   SUMMARY:
%       Turns a one-row-per-event waveform table (ed_wvTbl) into the layout a
%       Prism XY "Mean, SD, N" table wants: the X column is time, each COLUMN
%       is a mouse, and every mouse carries three subcolumns - Mean, SD, N -
%       computed across that mouse's events at each time point. Paste it into
%       an XY table formatted "Enter and plot error... Mean, SD, N" and Prism
%       draws the mean waveform with its spread.
%
%       ONE BLOCK PER SPLIT LEVEL. With 'splitVar' (genotype), each level is a
%       separate block with its own mice as columns - a Prism table each, which
%       is how the figure is built, one graph per genotype. The block named by
%       'copy' goes to the clipboard; the rest are returned to paste in turn.
%
%       It only summarises. Detrending and alignment are ed_wvTbl's job and are
%       already in the waveforms it hands over, so a mouse's mean here is the
%       trace its guiTbl_xy tile shows.
%
%       N is counted PER TIME POINT (finite samples), so an event too near a
%       recording edge to reach a given time simply is not counted there rather
%       than dragging the mean. A mouse with one event has an undefined SD; that
%       cell is left blank, and Prism still plots the mean.
%
%   INPUTS:
%       tbl      - (Table) one row per event; needs a matrix waveform column
%                          and the grouping column(s) below.
%       xVec     - (Vec)   [1 x nSamp] x axis shared by all rows, e.g. tstamps
%                          in ms. Length must match the waveform width.
%       varargin - Parameter/Value:
%           'yVar'     - (Char) waveform column. {'lfp'}
%           'scale'    - (Num)  multiply the waveforms by this before
%                               summarising. The default 1/1000 converts the LFP
%                               from uV to mV (wv2prism's origin); pass 1 for
%                               data that is already unitless (a normalised PETH).
%           'grpVar'   - (Char) column -> data columns (one mouse each).
%                               {'sbjID'}
%           'splitVar' - (Char) column -> one block each (one per genotype).
%                               '' = a single block over the whole table. {''}
%           'xLbl'     - (Char) header for the X column. {'x'}
%           'xLim'     - (Vec)  [lo hi] to crop the x axis before export, in
%                               xVec units. [] = keep all. {[]}
%           'copy'     - (Char) split level to put on the clipboard now.
%                               '' = the first. {''}
%           'flgSort'  - (Log)  natsort the columns and the blocks. {true}
%
%   OUTPUT:
%       blocks - (Struct array) one per split level:
%           .name - (Char)   the level ('all' when no splitVar).
%           .str  - (Char)   the tab-separated block, ready to paste.
%           .cols - (Cellstr) mouse ids, in column order.
%
%   DEPENDENCIES:
%       natsort (optional; only for ordering).
%
%   HISTORY:
%       260722 created for the per-genotype waveform figure in mcu_ed.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'tbl', @istable);
addRequired(p, 'xVec', @isnumeric);
addParameter(p, 'yVar', 'lfp', @ischar);
addParameter(p, 'scale', 1/1000, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'grpVar', 'sbjID', @ischar);
addParameter(p, 'splitVar', '', @ischar);
addParameter(p, 'xLbl', 'x', @ischar);
addParameter(p, 'xLim', [], @(x) isempty(x) || numel(x) == 2);
addParameter(p, 'copy', '', @ischar);
addParameter(p, 'flgSort', true, @islogical);
addParameter(p, 'flgClip', true, @islogical);
parse(p, tbl, xVec, varargin{:});
yVar     = p.Results.yVar;
scale    = p.Results.scale;
grpVar   = p.Results.grpVar;
splitVar = p.Results.splitVar;
xLbl     = p.Results.xLbl;
xLim     = p.Results.xLim;
whichCp  = p.Results.copy;
flgSort  = p.Results.flgSort;
flgClip  = p.Results.flgClip;

x = xVec(:)';
if size(tbl.(yVar), 2) ~= numel(x)
    error('wv2prism:xLen', 'xVec (%d) does not match %s width (%d)', ...
        numel(x), yVar, size(tbl.(yVar), 2));
end

% crop the x axis once, up front, so every block shares it
iKeep = true(1, numel(x));
if ~isempty(xLim)
    iKeep = x >= xLim(1) & x <= xLim(2);
end
x = x(iKeep);

%% ========================================================================
%  ONE BLOCK PER SPLIT LEVEL
%  ========================================================================
if isempty(splitVar)
    levels = {'all'};
    splitOf = ones(height(tbl), 1);
else
    [levels, ~, splitOf] = uniqueSorted(tbl.(splitVar), flgSort);
end

blocks = struct('name', {}, 'str', {}, 'cols', {});
for iL = 1 : numel(levels)
    sub = tbl(splitOf == iL, :);
    [mice, ~, grpOf] = uniqueSorted(sub.(grpVar), flgSort);

    % [nSamp x nMouse] mean / sd / n, each mouse reduced over its own events
    nS = numel(x);
    mu = nan(nS, numel(mice));
    sd = nan(nS, numel(mice));
    nn = zeros(nS, numel(mice));
    for iM = 1 : numel(mice)
        W = double(sub.(yVar)(grpOf == iM, iKeep)) * scale; % uV->mV by default
        mu(:, iM) = mean(W, 1, 'omitnan');
        sd(:, iM) = std(W, 0, 1, 'omitnan');
        nn(:, iM) = sum(isfinite(W), 1);
    end
    % std of a single sample is 0, not NaN, which Prism would draw as a
    % zero-length error bar - a fake certainty. Blank it where n < 2.
    sd(nn < 2) = NaN;

    blocks(iL).name = char(levels(iL));
    blocks(iL).cols = cellstr(string(mice(:)))';
    blocks(iL).str  = buildStr(x, mu, sd, nn, blocks(iL).cols, xLbl);
end

%% ========================================================================
%  CLIPBOARD
%  ========================================================================
if flgClip && ~isempty(blocks)
    iCp = 1;
    if ~isempty(whichCp)
        hit = find(strcmp({blocks.name}, whichCp), 1);
        if ~isempty(hit), iCp = hit; end
    end
    clipboard('copy', blocks(iCp).str);
    fprintf('wv2prism: copied "%s" (%d mice). Others: %s\n', ...
        blocks(iCp).name, numel(blocks(iCp).cols), ...
        strjoin(setdiff({blocks.name}, blocks(iCp).name), ', '));
end

end     % EOF


% =========================================================================
%  LOCAL
% =========================================================================
function [u, ia, ic] = uniqueSorted(col, flgSort)
% Unique levels of a column as a string array, natsorted when asked.
[u, ia, ic] = unique(col, 'stable');
u = string(u);
if flgSort && exist('natsort', 'file')
    [~, ord] = ismember(natsort(cellstr(u)), cellstr(u));
    u = u(ord);
    remap(ord) = 1 : numel(ord);
    ic = remap(ic)';
    ia = ia(ord);
end

end     % uniqueSorted


function s = buildStr(x, mu, sd, nn, cols, xLbl)
% Tab-separated for a Prism XY "Mean, SD, N" paste: ONE title row - the X
% label, then each mouse name over the first of its three subcolumns - and one
% row per time point (x, then mean/sd/n per mouse). Only the names are a
% header; Prism labels the Mean/SD/N subcolumns itself once the table is
% formatted that way, so a second header row here would misalign the paste.
nl = newline;
tab = sprintf('\t');

h1 = xLbl;
for iM = 1 : numel(cols)
    h1 = [h1, tab, cols{iM}, tab, tab]; %#ok<AGROW>
end

lines = cell(numel(x), 1);
for iX = 1 : numel(x)
    row = num2str(x(iX), '%.6g');
    for iM = 1 : size(mu, 2)
        row = [row, tab, fmt(mu(iX, iM)), tab, fmt(sd(iX, iM)), ...
            tab, fmt(nn(iX, iM))]; %#ok<AGROW>
    end
    lines{iX} = row;
end

s = [h1, nl, strjoin(lines, nl), nl];

end     % buildStr


function t = fmt(v)
% A blank for a non-finite cell (e.g. SD of a single event), else %.6g.
if ~isfinite(v)
    t = '';
else
    t = num2str(v, '%.6g');
end

end     % fmt
