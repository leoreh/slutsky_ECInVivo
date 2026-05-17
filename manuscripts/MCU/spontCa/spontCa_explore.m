function hFig = spontCa_explore(tblEvent, tblCell, fs, varargin)
% SPONTCA_EXPLORE  Two-panel viewer for spontCa pipeline tables.
%
%   hFig = spontCa_explore(tblEvent, tblCell, fs, ...) opens a figure
%   with a thin top control bar driving two side-by-side
%   tblGUI_scatHist widgets.
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
%   Pair filter:
%       Event level. Filters tblEvent by the upstream `paired` boolean
%       set by spontCa2_metrics. 'paired' keeps paired==true,
%       'unpaired' keeps paired==false, 'all' keeps everything.
%       Cell level. Re-aggregates per-cell metrics from the filtered
%       event set by calling spontCa2_metrics(mode='cellOnly'). Results
%       cached per mode for the lifetime of the figure. Cell pairs
%       themselves are not filtered; only the per-cell aggregates shift.
%
%   INPUTS
%       tblEvent - (table) event-level table from spontCa2_metrics. Must
%                  carry `paired` and `pairIdx`.
%       tblCell  - (table) cell-level aggregate from spontCa2_metrics.
%       fs       - (scalar Hz) sampling rate. Used for cell-level
%                  re-aggregation on pair-filter change.
%
%   OPTIONAL KEY-VALUE PAIRS
%       'level'       - 'event' (default) | 'cell'
%       'compartment' - 'cyto'  (default) | 'mito'
%       'pairFilter'  - 'all'   (default) | 'paired' | 'unpaired'
%       'xMetric'     - initial Panel-2 X metric. Default 'amp'.
%       'yMetric'     - initial Panel-2 Y metric. Default 'amp'.
%       'clr'         - 2x3 RGB for genotype groups.
%       'figPos'      - [x y w h] figure position.
%
%   See also TBLGUI_SCATHIST, MCU_SPONTCA, SPONTCA2_METRICS.


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

hFig = figure('Name', 'spontCa explore', 'NumberTitle', 'off', ...
    'Color', 'w', 'Position', p.Results.figPos, 'MenuBar', 'none', ...
    'ToolBar', 'figure');

ctrlH = 0.06;
hPnlCtrl = uipanel('Parent', hFig, 'BorderType', 'none', ...
    'Units', 'normalized', 'Position', [0, 1 - ctrlH, 1, ctrlH], ...
    'BackgroundColor', get(hFig, 'Color'));

gap = 0.005;
hPnlLeft  = uipanel('Parent', hFig, 'BorderType', 'none', ...
    'Units', 'normalized', 'Position', [0, 0, 0.5 - gap, 1 - ctrlH], ...
    'BackgroundColor', get(hFig, 'Color'));
hPnlRight = uipanel('Parent', hFig, 'BorderType', 'none', ...
    'Units', 'normalized', 'Position', [0.5 + gap, 0, 0.5 - gap, 1 - ctrlH], ...
    'BackgroundColor', get(hFig, 'Color'));


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
state.levelOpts    = {'event', 'cell'};
state.compOpts     = {'cyto',  'mito'};
state.pairOpts     = {'all',   'paired', 'unpaired'};
state.cellTblCache = struct('all', tblCell, 'paired', [], 'unpaired', []);

yRow = 0.25; rowH = 0.5;

% Level
uicontrol('Parent', hPnlCtrl, 'Style', 'text', 'String', 'Level:', ...
    'Units', 'normalized', 'Position', [0.01, yRow, 0.04, rowH], ...
    'HorizontalAlignment', 'right', 'FontWeight', 'bold', ...
    'BackgroundColor', get(hFig, 'Color'));
state.ddLevel = uicontrol('Parent', hPnlCtrl, 'Style', 'popupmenu', ...
    'String', {'Event', 'Cell'}, 'Units', 'normalized', ...
    'Position', [0.06, yRow, 0.06, rowH], ...
    'Callback', @(s, e) onLevelChange(hFig));

% Compartment
uicontrol('Parent', hPnlCtrl, 'Style', 'text', 'String', 'Compartment:', ...
    'Units', 'normalized', 'Position', [0.13, yRow, 0.07, rowH], ...
    'HorizontalAlignment', 'right', 'FontWeight', 'bold', ...
    'BackgroundColor', get(hFig, 'Color'));
state.ddComp = uicontrol('Parent', hPnlCtrl, 'Style', 'popupmenu', ...
    'String', {'Cyto', 'Mito'}, 'Units', 'normalized', ...
    'Position', [0.21, yRow, 0.06, rowH], ...
    'Callback', @(s, e) onCompChange(hFig));

% Pair filter
uicontrol('Parent', hPnlCtrl, 'Style', 'text', 'String', 'Pair:', ...
    'Units', 'normalized', 'Position', [0.28, yRow, 0.04, rowH], ...
    'HorizontalAlignment', 'right', 'FontWeight', 'bold', ...
    'BackgroundColor', get(hFig, 'Color'));
state.ddPair = uicontrol('Parent', hPnlCtrl, 'Style', 'popupmenu', ...
    'String', {'All', 'Paired', 'Unpaired'}, 'Units', 'normalized', ...
    'Position', [0.33, yRow, 0.07, rowH], ...
    'Callback', @(s, e) onPairChange(hFig));

% Panel 2 X-metric
uicontrol('Parent', hPnlCtrl, 'Style', 'text', 'String', 'P2 X:', ...
    'Units', 'normalized', 'Position', [0.41, yRow, 0.04, rowH], ...
    'HorizontalAlignment', 'right', 'FontWeight', 'bold', ...
    'BackgroundColor', get(hFig, 'Color'), ...
    'TooltipString', 'Cyto-axis metric on Panel 2');
state.ddXMetric = uicontrol('Parent', hPnlCtrl, 'Style', 'popupmenu', ...
    'String', {' '}, 'Units', 'normalized', ...
    'Position', [0.46, yRow, 0.09, rowH], ...
    'Callback', @(s, e) onMetricChange(hFig));

% Panel 2 Y-metric
uicontrol('Parent', hPnlCtrl, 'Style', 'text', 'String', 'P2 Y:', ...
    'Units', 'normalized', 'Position', [0.56, yRow, 0.04, rowH], ...
    'HorizontalAlignment', 'right', 'FontWeight', 'bold', ...
    'BackgroundColor', get(hFig, 'Color'), ...
    'TooltipString', 'Mito-axis metric on Panel 2');
state.ddYMetric = uicontrol('Parent', hPnlCtrl, 'Style', 'popupmenu', ...
    'String', {' '}, 'Units', 'normalized', ...
    'Position', [0.61, yRow, 0.09, rowH], ...
    'Callback', @(s, e) onMetricChange(hFig));

% Status
state.hStatus = uicontrol('Parent', hPnlCtrl, 'Style', 'text', 'String', '', ...
    'Units', 'normalized', 'Position', [0.71, yRow, 0.28, rowH], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold', ...
    'BackgroundColor', get(hFig, 'Color'));

set(state.ddLevel, 'Value', find(strcmpi(p.Results.level,       state.levelOpts)));
set(state.ddComp,  'Value', find(strcmpi(p.Results.compartment, state.compOpts)));
set(state.ddPair,  'Value', find(strcmpi(p.Results.pairFilter,  state.pairOpts)));

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
    level    = state.levelOpts{get(state.ddLevel, 'Value')};
    comp     = state.compOpts {get(state.ddComp,  'Value')};
    pairMode = state.pairOpts {get(state.ddPair,  'Value')};

    % Resolve current source tables. Panel 1 uses the filter-aware
    % source. Panel 2 at event level always uses the full tblEvent so
    % pairIdx values (absolute row indices) resolve correctly; at cell
    % level it uses the filter-aware cell table.
    if strcmpi(level, 'event')
        srcEvent = applyPairFilter(state.tblEvent, pairMode);
        srcCell  = [];
    else
        srcEvent = [];
        srcCell  = getCellTbl(hFig, pairMode);
        state    = hFig.UserData;   % refresh after possible cache write
    end

    % Refresh metric option lists on full rebuilds.
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
        tblGUI_scatHist(tbl1, ...
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
        tblGUI_scatHist(tbl2, ...
            'Parent', state.hPnlRight, ...
            'xVar', x2, 'yVar', y2, 'grpVar', g2, ...
            'xScale', snap.xs, 'yScale', snap.ys, 'fitType', snap.ft, ...
            'clr', state.clr, ...
            'varsExclude', excl2);
        n2 = height(tbl2);
    elseif isstruct(state.hPnlRight.UserData) && isfield(state.hPnlRight.UserData, 'tbl')
        n2 = height(state.hPnlRight.UserData.tbl);
    end

    set(state.hStatus, 'String', sprintf( ...
        'Level: %s | Comp: %s | Pair: %s | P2: cyto.%s vs mito.%s | n1=%d  n2=%d', ...
        level, comp, pairMode, xMetric, yMetric, n1, n2));

    hFig.UserData = state;
end


function populateMetricLists(hFig, src)
% Both X-metric and Y-metric dropdowns share the same numeric-metric
% list. Preserves the current selection if still available, else falls
% back to pendingX / pendingY (used on first render), 'amp', or first.
    state   = hFig.UserData;
    metrics = listMetrics(src);

    for ddField = {'ddXMetric', 'ddYMetric'}
        dd = state.(ddField{1});
        items = get(dd, 'String');
        current = '';
        if iscell(items) && ~isempty(items) && get(dd, 'Value') <= numel(items)
            current = items{get(dd, 'Value')};
        end
        if strcmpi(ddField{1}, 'ddXMetric'),  pending = state.pendingX;
        else,                                   pending = state.pendingY;
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
        set(dd, 'String', metrics, 'Value', find(strcmp(metrics, pick), 1));
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
% xs/ys default to 'log' before the first render so initial display
% uses log scales; once the widget exists, the user's choice is carried.
    snap = struct('x', '', 'y', '', 'g', '', 'xs', 'log', 'ys', 'log', 'ft', '');
    if ~isgraphics(hPanel), return; end
    d = hPanel.UserData;
    if ~isstruct(d) || ~isfield(d, 'ddX') || ~isgraphics(d.ddX)
        return;
    end
    try
        snap.x = d.numericVars{get(d.ddX, 'Value')};
        snap.y = d.numericVars{get(d.ddY, 'Value')};
        items  = get(d.ddGrp, 'String');
        snap.g = items{get(d.ddGrp, 'Value')};
        if strcmp(snap.g, 'None'), snap.g = ''; end

        xs = get(d.ddXScale, 'String'); snap.xs = xs{get(d.ddXScale, 'Value')};
        ys = get(d.ddYScale, 'String'); snap.ys = ys{get(d.ddYScale, 'Value')};
        ft = get(d.ddFit,    'String'); snap.ft = ft{get(d.ddFit,    'Value')};
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
% Non-pair defaults so unpaired events (NaN pairAmp / pairFlux under
% mutual NN) still appear in the scatter.
    defXY = struct('x', 'amp', 'y', 'flux');
end

function excl = idColsToHide()
    excl = {'sbjID', 'unitID', 'excluded', 'trace', ...
            'start', 'stop', 'pairIdx'};
end

function m = pickMetric(dd)
    items = get(dd, 'String');
    val   = get(dd, 'Value');
    if iscell(items) && val >= 1 && val <= numel(items)
        m = items{val};
    else
        m = '';
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
% Lazy cache of cell-level aggregates per pair-filter mode. 'all' is
% seeded at init from the input tblCell. 'paired' and 'unpaired' are
% computed on first access via spontCa2_metrics(mode='cellOnly') on the
% filtered event set.
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
% Numeric columns on src minus identifiers and the trace blob. Both
% Panel-2 X and Y dropdowns populate from this list.
    cols  = src.Properties.VariableNames;
    isNum = varfun(@isnumeric, src, 'OutputFormat', 'uniform');
    skip  = [idColsToHide(), {'compartment'}];
    metrics = setdiff(cols(isNum), skip, 'stable');
    if isempty(metrics), metrics = {' '}; end
end

function tbl1 = prepPanel1(srcEvent, srcCell, level, comp)
% Filter source table to the requested compartment. At event level the
% source is already pair-filtered upstream in refreshPanels.
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
%
% Event level: pass the FULL event table (state.tblEvent), not the
% pair-filtered subset. pairIdx values are absolute row indices into
% the original tblEvent and would not resolve against a subset. Panel
% 2 at event level is the paired-events view by construction.
%
% Cell level: srcCell is the (filter-aware) cell table; cyto and mito
% rows are joined by sbjID.

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
