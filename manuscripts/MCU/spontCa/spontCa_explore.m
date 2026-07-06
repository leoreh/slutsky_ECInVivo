function hFig = spontCa_explore(tblEvent, tblCell, fs, varargin)
% SPONTCA_EXPLORE  Two-panel viewer for spontCa pipeline tables.
%
%   hFig = spontCa_explore(tblEvent, tblCell, fs, ...) opens a figure
%   with a thin top control bar driving two side-by-side
%   guiTbl_scatHist widgets.
%
%       Panel 1 (left)  : within-compartment scatter+hist. X / Y / group
%                         pickers live inside the inner widget.
%       Panel 2 (right) : cross-compartment scatter+hist with fixed
%                         compartment axes (X = Cyto.<P2-X-Metric>,
%                         Y = Mito.<P2-Y-Metric>). At event level one
%                         row per paired cyto event (partner via
%                         pairIdx). At cell level one row per sbjID
%                         (Cyto and Mito rows of tblCell joined on
%                         sbjID).
%
%   Top-bar controls:
%       Level        : Event | Cell
%       Compartment  : Cyto  | Mito                (Panel 1 only)
%       Pair filter  : All | Paired | Unpaired     (global state)
%       P2 X-Metric  : Cyto-axis metric on Panel 2
%       P2 Y-Metric  : Mito-axis metric on Panel 2
%
%   Built on the shared graphics/gui layer (uifigure); embeds
%   guiTbl_scatHist into the two panels.
%
%   INPUTS
%       tblEvent - (table) event-level table from spontCa2_metrics. Must
%                  carry `paired` and `pairIdx`.
%       tblCell  - (table) cell-level aggregate from spontCa2_metrics.
%       fs       - (scalar Hz) sampling rate.
%
%   OPTIONAL KEY-VALUE PAIRS
%       'level'/'compartment'/'pairFilter'/'xMetric'/'yMetric'/'clr'/'figPos'.
%
%   See also GUITBL_SCATHIST, MCU_SPONTCA, SPONTCA2_METRICS.


%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tblEvent', @istable);
addRequired(p, 'tblCell',  @istable);
addRequired(p, 'fs',       @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'level',       'event', @(x) any(strcmpi(x, {'event', 'cell'})));
addParameter(p, 'compartment', 'cyto',  @(x) any(strcmpi(x, {'cyto', 'mito'})));
addParameter(p, 'pairFilter',  'all',   @(x) any(strcmpi(x, {'all', 'paired', 'unpaired'})));
addParameter(p, 'xMetric',     'amp',   @(x) ischar(x) || isstring(x));
addParameter(p, 'yMetric',     'amp',   @(x) ischar(x) || isstring(x));
addParameter(p, 'clr',         [],      @(x) isempty(x) || (isnumeric(x) && size(x,2)==3));
addParameter(p, 'figPos',      [50, 50, 1800, 950], @(x) isnumeric(x) && numel(x)==4);
parse(p, tblEvent, tblCell, fs, varargin{:});

clr = p.Results.clr;
if isempty(clr)
    try
        cfg = mcu_cfg();
        clr = cfg.clr.grp;
    catch
        clr = [0.2 0.2 0.2; 0.75 0.55 0.35];
    end
end


%% ========================================================================
%  FIGURE
%  ========================================================================

hFig = uifigure('Name', 'spontCa explore', 'Position', p.Results.figPos);

gMain = uigridlayout(hFig, [2, 1], 'RowHeight', {38, '1x'}, ...
    'Padding', [4 4 4 4], 'RowSpacing', 4);

gCtrl = uigridlayout(gMain, [1, 11], ...
    'ColumnWidth', {40, 75, 95, 70, 35, 85, 45, 95, 45, 95, '1x'}, ...
    'Padding', [2 2 2 2], 'ColumnSpacing', 4);
gCtrl.Layout.Row = 1;

gPanels = uigridlayout(gMain, [1, 2], 'ColumnWidth', {'1x', '1x'}, ...
    'Padding', [0 0 0 0], 'ColumnSpacing', 6);
gPanels.Layout.Row = 2;
hPnlLeft  = uipanel(gPanels, 'BorderType', 'none');
hPnlLeft.Layout.Column = 1;
hPnlRight = uipanel(gPanels, 'BorderType', 'none');
hPnlRight.Layout.Column = 2;


%% ========================================================================
%  STATE + CONTROLS
%  ========================================================================

state              = struct();
state.tblEvent     = tblEvent;
state.tblCellAll   = tblCell;
state.fs           = fs;
state.clr          = clr;
state.hPnlLeft     = hPnlLeft;
state.hPnlRight    = hPnlRight;
state.cellTblCache = struct('all', tblCell, 'paired', [], 'unpaired', []);

% ItemsData carry the lowercase option keys, so dropdown .Value reads return
% them directly (no index bookkeeping).
uilabel(gCtrl, 'Text', 'Level:', 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
state.ddLevel = uidropdown(gCtrl, 'Items', {'Event', 'Cell'}, ...
    'ItemsData', {'event', 'cell'}, 'ValueChangedFcn', @(~, ~) onLevelChange(hFig));

uilabel(gCtrl, 'Text', 'Compartment:', 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
state.ddComp = uidropdown(gCtrl, 'Items', {'Cyto', 'Mito'}, ...
    'ItemsData', {'cyto', 'mito'}, 'ValueChangedFcn', @(~, ~) onCompChange(hFig));

uilabel(gCtrl, 'Text', 'Pair:', 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
state.ddPair = uidropdown(gCtrl, 'Items', {'All', 'Paired', 'Unpaired'}, ...
    'ItemsData', {'all', 'paired', 'unpaired'}, 'ValueChangedFcn', @(~, ~) onPairChange(hFig));

uilabel(gCtrl, 'Text', 'P2 X:', 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
state.ddXMetric = uidropdown(gCtrl, 'Items', {' '}, 'ValueChangedFcn', @(~, ~) onMetricChange(hFig));

uilabel(gCtrl, 'Text', 'P2 Y:', 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
state.ddYMetric = uidropdown(gCtrl, 'Items', {' '}, 'ValueChangedFcn', @(~, ~) onMetricChange(hFig));

state.hStatus = uilabel(gCtrl, 'Text', '', 'FontWeight', 'bold');

state.ddLevel.Value = lower(p.Results.level);
state.ddComp.Value  = lower(p.Results.compartment);
state.ddPair.Value  = lower(p.Results.pairFilter);

state.pendingX = char(p.Results.xMetric);
state.pendingY = char(p.Results.yMetric);

hFig.UserData = state;
refreshPanels(hFig, 'all');

end % EOF spontCa_explore


%% ========================================================================
%  CALLBACKS
%  ========================================================================

function onLevelChange(hFig)
    refreshPanels(hFig, 'all');
end

function onCompChange(hFig)
    refreshPanels(hFig, 'panel1');
end

function onPairChange(hFig)
    refreshPanels(hFig, 'all');
end

function onMetricChange(hFig)
    refreshPanels(hFig, 'panel2');
end


%% ========================================================================
%  REFRESH
%  ========================================================================

function refreshPanels(hFig, which)
% which: 'all' | 'panel1' | 'panel2'

    state    = hFig.UserData;
    level    = state.ddLevel.Value;
    comp     = state.ddComp.Value;
    pairMode = state.ddPair.Value;

    % Resolve current source tables. Panel 1 uses the filter-aware source.
    % Panel 2 at event level always uses the full tblEvent so pairIdx values
    % (absolute row indices) resolve correctly; at cell level it uses the
    % filter-aware cell table.
    if strcmpi(level, 'event')
        srcEvent = applyPairFilter(state.tblEvent, pairMode);
        srcCell  = [];
    else
        srcEvent = [];
        srcCell  = getCellTbl(hFig, pairMode);
        state    = hFig.UserData;   % refresh after possible cache write
    end

    if strcmpi(which, 'all')
        if strcmpi(level, 'event')
            populateMetricLists(hFig, state.tblEvent);
        else
            populateMetricLists(hFig, srcCell);
        end
        state = hFig.UserData;
    end

    xMetric = pickMetric(state.ddXMetric);
    yMetric = pickMetric(state.ddYMetric);

    % --- Panel 1 ---
    n1 = NaN;
    if any(strcmpi(which, {'all', 'panel1'}))
        snap   = readPanelState(state.hPnlLeft);
        tbl1   = prepPanel1(srcEvent, srcCell, level, comp);
        defXY1 = pickDefaultsPanel1(level);
        excl1  = idColsToHide();
        [x1, y1, g1] = pickXYG(tbl1, snap.x, snap.y, snap.g, defXY1, excl1);

        delete(allchild(state.hPnlLeft));
        guiTbl_scatHist(tbl1, ...
            'Parent', state.hPnlLeft, ...
            'xVar', x1, 'yVar', y1, 'grpVar', g1, ...
            'xScale', snap.xs, 'yScale', snap.ys, 'fitType', snap.ft, ...
            'clr', state.clr, ...
            'varsExclude', excl1);
        n1 = height(tbl1);
    elseif isstruct(state.hPnlLeft.UserData) && isfield(state.hPnlLeft.UserData, 'tbl')
        n1 = height(state.hPnlLeft.UserData.tbl);
    end

    % --- Panel 2 ---
    n2 = NaN;
    if any(strcmpi(which, {'all', 'panel2'}))
        snap   = readPanelState(state.hPnlRight);
        if strcmpi(level, 'event')
            tbl2 = prepPanel2(state.tblEvent, [], level, xMetric, yMetric);
        else
            tbl2 = prepPanel2([], srcCell, level, xMetric, yMetric);
        end
        defXY2 = struct('x', 'cytoX', 'y', 'mitoY');
        excl2  = idColsToHide();
        [x2, y2, g2] = pickXYG(tbl2, snap.x, snap.y, snap.g, defXY2, excl2);

        delete(allchild(state.hPnlRight));
        guiTbl_scatHist(tbl2, ...
            'Parent', state.hPnlRight, ...
            'xVar', x2, 'yVar', y2, 'grpVar', g2, ...
            'xScale', snap.xs, 'yScale', snap.ys, 'fitType', snap.ft, ...
            'clr', state.clr, ...
            'varsExclude', excl2);
        n2 = height(tbl2);
    elseif isstruct(state.hPnlRight.UserData) && isfield(state.hPnlRight.UserData, 'tbl')
        n2 = height(state.hPnlRight.UserData.tbl);
    end

    state.hStatus.Text = sprintf( ...
        'Level: %s | Comp: %s | Pair: %s | P2: cyto.%s vs mito.%s | n1=%d  n2=%d', ...
        level, comp, pairMode, xMetric, yMetric, n1, n2);

    hFig.UserData = state;
end


function populateMetricLists(hFig, src)
% Both X-metric and Y-metric dropdowns share the same numeric-metric list.
% Preserves the current selection if still available, else pendingX/pendingY
% (first render), 'amp', or first.
    state   = hFig.UserData;
    metrics = listMetrics(src);

    for ddField = {'ddXMetric', 'ddYMetric'}
        dd = state.(ddField{1});
        current = '';
        if iscell(dd.Items) && ~isempty(dd.Items)
            current = dd.Value;
        end
        if strcmpi(ddField{1}, 'ddXMetric'), pending = state.pendingX;
        else,                                pending = state.pendingY;
        end

        if ismember(current, metrics)
            pick = current;
        elseif ~isempty(pending) && ismember(pending, metrics)
            pick = pending;
        elseif ismember('amp', metrics)
            pick = 'amp';
        else
            pick = metrics{1};
        end
        dd.Items = metrics;
        dd.Value = pick;
    end

    state.pendingX = '';
    state.pendingY = '';
    hFig.UserData = state;
end


%% ========================================================================
%  HELPERS: STATE SNAPSHOT, DEFAULTS, EXCLUSIONS
%  ========================================================================

function snap = readPanelState(hPanel)
% Snap the embedded widget's current selections (X/Y/group/scales/fit).
% xs/ys default to 'log' before the first render so initial display uses log
% scales; once the widget exists, the user's choice is carried.
    snap = struct('x', '', 'y', '', 'g', '', 'xs', 'log', 'ys', 'log', 'ft', '');
    if ~isgraphics(hPanel), return; end
    d = hPanel.UserData;
    if ~isstruct(d) || ~isfield(d, 'ddX') || ~isgraphics(d.ddX)
        return;
    end
    try
        snap.x  = d.ddX.Value;
        snap.y  = d.ddY.Value;
        snap.g  = d.ddGrp.Value;
        if strcmp(snap.g, 'None'), snap.g = ''; end
        snap.xs = d.ddXScale.Value;
        snap.ys = d.ddYScale.Value;
        snap.ft = d.ddFit.Value;
    catch
        snap = struct('x', '', 'y', '', 'g', '', 'xs', 'log', 'ys', 'log', 'ft', '');
    end
end

function [x, y, g] = pickXYG(tbl, oldX, oldY, oldG, defXY, excl)
    numericMask = varfun(@isnumeric, tbl, 'OutputFormat', 'uniform');
    cols        = tbl.Properties.VariableNames;
    numCols     = setdiff(cols(numericMask), excl, 'stable');

    catMask = varfun(@(v) iscategorical(v) || isstring(v) || islogical(v), ...
        tbl, 'OutputFormat', 'uniform');
    catCols = setdiff(cols(catMask), excl, 'stable');

    x = pickOne(oldX, defXY.x, numCols);
    y = pickOne(oldY, defXY.y, numCols);
    g = pickOne(oldG, 'genotype', catCols);
end

function v = pickOne(carried, fallback, allowed)
    if ~isempty(carried) && ismember(carried, allowed)
        v = carried;
    elseif ismember(fallback, allowed)
        v = fallback;
    else
        v = '';
    end
end

function defXY = pickDefaultsPanel1(~)
% Non-pair defaults so unpaired events (NaN pairAmp / pairFlux under mutual
% NN) still appear in the scatter.
    defXY = struct('x', 'amp', 'y', 'flux');
end

function excl = idColsToHide()
    excl = {'sbjID', 'unitID', 'excluded', 'trace', ...
            'start', 'stop', 'pairIdx'};
end

function m = pickMetric(dd)
    if isempty(dd.Items) || (numel(dd.Items) == 1 && strcmp(dd.Items{1}, ' '))
        m = '';
    else
        m = dd.Value;
    end
end


%% ========================================================================
%  HELPERS: PAIR FILTER + CELL CACHE
%  ========================================================================

function evt = applyPairFilter(tblEvent, mode)
    switch lower(mode)
        case 'paired'
            evt = tblEvent(tblEvent.paired, :);
        case 'unpaired'
            evt = tblEvent(~tblEvent.paired, :);
        otherwise
            evt = tblEvent;
    end
end

function tbl = getCellTbl(hFig, mode)
% Lazy cache of cell-level aggregates per pair-filter mode.
    state = hFig.UserData;
    if isfield(state.cellTblCache, mode) && ~isempty(state.cellTblCache.(mode))
        tbl = state.cellTblCache.(mode);
        return;
    end
    evt = applyPairFilter(state.tblEvent, mode);
    [tbl, ~] = spontCa2_metrics(state.tblCellAll, evt, state.fs, ...
        'mode', 'cellOnly', 'aggFcn', 'mean');
    state.cellTblCache.(mode) = tbl;
    hFig.UserData = state;
end


%% ========================================================================
%  HELPERS: METRICS + TABLE PREPARATION
%  ========================================================================

function metrics = listMetrics(src)
    cols  = src.Properties.VariableNames;
    isNum = varfun(@isnumeric, src, 'OutputFormat', 'uniform');
    skip  = [idColsToHide(), {'compartment'}];
    metrics = setdiff(cols(isNum), skip, 'stable');
    if isempty(metrics), metrics = {' '}; end
end

function tbl1 = prepPanel1(srcEvent, srcCell, level, comp)
    if strcmpi(level, 'event')
        src = srcEvent;
    else
        src = srcCell;
    end
    if strcmpi(comp, 'cyto')
        keep = src.compartment == 'Cyto';
    else
        keep = src.compartment == 'Mito';
    end
    tbl1 = src(keep, :);
end

function tbl2 = prepPanel2(fullEvent, srcCell, level, xMetric, yMetric)
% Wide layout for cross-compartment scatter. X = Cyto.<xMetric>,
% Y = Mito.<yMetric>.
    if strcmpi(level, 'event')
        focal      = find(fullEvent.compartment == 'Cyto' & fullEvent.paired);
        partnerIdx = fullEvent.pairIdx(focal);
        tbl2 = table();
        tbl2.sbjID    = fullEvent.sbjID(focal);
        tbl2.genotype = fullEvent.genotype(focal);
        tbl2.cytoX    = fullEvent.(xMetric)(focal);
        tbl2.mitoY    = fullEvent.(yMetric)(partnerIdx);
    else
        isC = srcCell.compartment == 'Cyto';
        isM = srcCell.compartment == 'Mito';
        cytoTbl = srcCell(isC, :);
        mitoTbl = srcCell(isM, :);
        [~, idxM] = ismember(cytoTbl.sbjID, mitoTbl.sbjID);
        keep      = idxM > 0;
        cytoTbl   = cytoTbl(keep, :);
        mitoTbl   = mitoTbl(idxM(keep), :);
        tbl2 = table();
        tbl2.sbjID    = cytoTbl.sbjID;
        tbl2.genotype = cytoTbl.genotype;
        tbl2.cytoX    = cytoTbl.(xMetric);
        tbl2.mitoY    = mitoTbl.(yMetric);
    end
end
