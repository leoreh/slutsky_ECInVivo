function hFig = spontCa_explore(tblEvent, tblCell, varargin)
% SPONTCA_EXPLORE  Two-panel viewer for spontCa pipeline tables.
%
%   hFig = spontCa_explore(tblEvent, tblCell, ...) opens a figure with a
%   thin top control bar driving two side-by-side tblGUI_scatHist
%   widgets:
%
%       Panel 1 (left)  : within-compartment scatter+hist on the source
%                         table filtered to the selected compartment.
%                         All numeric metrics available, including pair-
%                         relationship metrics (pairFlux, pairAmp, tf,
%                         pairLag). Canonical view: amp vs pairAmp.
%
%       Panel 2 (right) : cross-compartment scatter+hist on a wide
%                         reshape. Each row is one paired (cyto, mito)
%                         observation; the table has exactly two metric
%                         columns - 'cyto' and 'mito' - both populated
%                         from the metric chosen in the top-bar Metric
%                         dropdown. Default X=cyto, Y=mito; changing
%                         Metric updates both axes together.
%
%   Top-bar controls:
%       Level       : Event | Cell                (table source)
%       Compartment : Cyto  | Mito                (Panel 1 + Panel 2)
%       Metric      : amp | dur | flux | ...      (Panel 2 only;
%                                                  symmetric metrics in
%                                                  the source table)
%
%   The Compartment selector drives both panels:
%       Level = Event, Comp = Cyto -> Panel 1 = Cyto events; Panel 2 =
%           cyto->mito pairs (focal = Cyto with valid pairIdx).
%       Level = Event, Comp = Mito -> Panel 1 = Mito events; Panel 2 =
%           mito->cyto pairs (focal = Mito with valid pairIdx).
%       Level = Cell , Comp = Cyto -> Panel 1 = Cyto rows of tblCell.
%       Level = Cell , Comp = Mito -> Panel 1 = Mito rows of tblCell.
%       At cell level, Panel 2 always pairs by sbjID (direction free) so
%       Compartment changes refresh Panel 1 only.
%
%   Symmetric-metric detection (Panel 2 + Metric dropdown):
%       A column is treated as asymmetric when its name starts with
%       'asymPrefix' (default 'pair') or appears in 'asymExtra' (default
%       {'tf', 'triggered'}). Adding a future pair_<x> metric auto-
%       flags; no per-name editing here.
%
%   State carried across control changes:
%       - X / Y / group selections inside each panel.
%       - Scale (linear/log) on each axis and fit type (None/Linear/
%         Ortho), independently per panel.
%       - Metric (across Compartment and Level changes, when still
%         available in the new source table).
%
%   INPUTS
%       tblEvent - (table) event-level table from PAIR & CROSS-FLUX in
%                  mcu_spontCa.m.
%       tblCell  - (table) cell-level aggregate from PER-CELL SUMMARY.
%
%   OPTIONAL KEY-VALUE PAIRS
%       'level'       - 'event' (default) | 'cell'
%       'compartment' - 'cyto'  (default) | 'mito'
%       'metric'      - initial Panel-2 metric (default: 'amp' if
%                       available, else the first symmetric metric).
%       'asymPrefix'  - char prefix flagging asymmetric columns. Default
%                       'pair'. Pass '' to disable prefix-based flagging.
%       'asymExtra'   - cellstr of extra asymmetric names. Default
%                       {'tf', 'triggered'}.
%       'clr'         - 2x3 RGB for genotype groups. Default pulled from
%                       mcu_cfg().clr.grp.
%       'figPos'      - [x y w h] figure position. Default
%                       [50 50 1800 950].
%
%   See also TBLGUI_SCATHIST, MCU_SPONTCA, MCU_CFG.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tblEvent', @istable);
addRequired(p, 'tblCell',  @istable);
addParameter(p, 'level',       'event', @(x) any(strcmpi(x, {'event', 'cell'})));
addParameter(p, 'compartment', 'cyto',  @(x) any(strcmpi(x, {'cyto', 'mito'})));
addParameter(p, 'metric',      '',      @(x) ischar(x) || isstring(x));
addParameter(p, 'asymPrefix',  'pair',  @(x) ischar(x) || isstring(x));
addParameter(p, 'asymExtra',   {'tf', 'triggered'}, @iscell);
addParameter(p, 'clr',         [],      @(x) isempty(x) || (isnumeric(x) && size(x,2)==3));
addParameter(p, 'figPos',      [50, 50, 1800, 950], @(x) isnumeric(x) && numel(x)==4);
parse(p, tblEvent, tblCell, varargin{:});

clr = p.Results.clr;
if isempty(clr)
    try
        cfg = mcu_cfg();
        clr = cfg.clr.grp;
    catch
        clr = [0.2 0.2 0.2; 0.75 0.55 0.35];
    end
end

% Defensive: attach 'triggered' to tblEvent if the caller skipped that step.
if ~ismember('triggered', tblEvent.Properties.VariableNames)
    isCyto = tblEvent.compartment == 'Cyto';
    trig = NaN(height(tblEvent), 1);
    trig(isCyto) = ~isnan(tblEvent.pairIdx(isCyto));
    tblEvent.triggered = categorical(trig, [0 1], {'noPair', 'triggered'});
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

state            = struct();
state.tblEvent   = tblEvent;
state.tblCell    = tblCell;
state.clr        = clr;
state.asymPrefix = char(p.Results.asymPrefix);
state.asymExtra  = p.Results.asymExtra;
state.hPnlLeft   = hPnlLeft;
state.hPnlRight  = hPnlRight;
state.levelOpts  = {'event', 'cell'};
state.compOpts   = {'cyto', 'mito'};

yRow = 0.25; rowH = 0.5;

% Level
uicontrol('Parent', hPnlCtrl, 'Style', 'text', 'String', 'Level:', ...
    'Units', 'normalized', 'Position', [0.01, yRow, 0.04, rowH], ...
    'HorizontalAlignment', 'right', 'FontWeight', 'bold', ...
    'BackgroundColor', get(hFig, 'Color'));
state.ddLevel = uicontrol('Parent', hPnlCtrl, 'Style', 'popupmenu', ...
    'String', {'Event', 'Cell'}, 'Units', 'normalized', ...
    'Position', [0.06, yRow, 0.07, rowH], ...
    'Callback', @(s, e) onLevelChange(hFig));

% Compartment
uicontrol('Parent', hPnlCtrl, 'Style', 'text', 'String', 'Compartment:', ...
    'Units', 'normalized', 'Position', [0.15, yRow, 0.08, rowH], ...
    'HorizontalAlignment', 'right', 'FontWeight', 'bold', ...
    'BackgroundColor', get(hFig, 'Color'));
state.ddComp = uicontrol('Parent', hPnlCtrl, 'Style', 'popupmenu', ...
    'String', {'Cyto', 'Mito'}, 'Units', 'normalized', ...
    'Position', [0.24, yRow, 0.07, rowH], ...
    'Callback', @(s, e) onCompChange(hFig));

% Panel-2 Metric. Options are computed lazily inside refreshPanels each
% time the source table changes (Level switch), so adding a new
% symmetric metric to tblEvent/tblCell surfaces here automatically.
uicontrol('Parent', hPnlCtrl, 'Style', 'text', 'String', 'P2 Metric:', ...
    'Units', 'normalized', 'Position', [0.33, yRow, 0.07, rowH], ...
    'HorizontalAlignment', 'right', 'FontWeight', 'bold', ...
    'BackgroundColor', get(hFig, 'Color'), ...
    'TooltipString', 'Metric loaded into Panel 2''s cyto and mito columns');
state.ddMetric = uicontrol('Parent', hPnlCtrl, 'Style', 'popupmenu', ...
    'String', {' '}, 'Units', 'normalized', ...
    'Position', [0.41, yRow, 0.10, rowH], ...
    'Callback', @(s, e) onMetricChange(hFig));

% Status
state.hStatus = uicontrol('Parent', hPnlCtrl, 'Style', 'text', 'String', '', ...
    'Units', 'normalized', 'Position', [0.53, yRow, 0.46, rowH], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold', ...
    'BackgroundColor', get(hFig, 'Color'));

set(state.ddLevel, 'Value', find(strcmpi(p.Results.level,       state.levelOpts)));
set(state.ddComp,  'Value', find(strcmpi(p.Results.compartment, state.compOpts)));

% Seed the metric list and selection. populateMetrics uses the level we
% just picked. The initial 'metric' input is honored if it's in the list.
state.pendingMetric = char(p.Results.metric);

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
    % At cell level, Panel 2 pairs by sbjID and doesn't depend on
    % compartment, so only Panel 1 needs a refresh.
    state = hFig.UserData;
    level = state.levelOpts{get(state.ddLevel, 'Value')};
    if strcmpi(level, 'cell')
        refreshPanels(hFig, 'panel1');
    else
        refreshPanels(hFig, 'all');
    end
end

function onMetricChange(hFig)
    % Metric drives Panel 2 only.
    refreshPanels(hFig, 'panel2');
end

%% ========================================================================
%  REFRESH
%  ========================================================================

function refreshPanels(hFig, which)
% which: 'all' rebuilds both panels;
%        'panel1' rebuilds only Panel 1;
%        'panel2' rebuilds only Panel 2.

    state = hFig.UserData;
    level = state.levelOpts{get(state.ddLevel, 'Value')};
    comp  = state.compOpts {get(state.ddComp,  'Value')};

    if strcmpi(level, 'event')
        src = state.tblEvent;
    else
        src = state.tblCell;
    end

    % Refresh the Metric dropdown's option list whenever the source
    % schema might have changed (Level switch or first render).
    if any(strcmpi(which, {'all'}))
        populateMetricList(hFig, src);
    end
    state = hFig.UserData;  % reload after populateMetricList may write
    metricOpts = get(state.ddMetric, 'String');
    metric     = metricOpts{get(state.ddMetric, 'Value')};

    % --- Panel 1 ---
    n1 = NaN;
    if any(strcmpi(which, {'all', 'panel1'}))
        snap = readPanelState(state.hPnlLeft);
        tbl1 = prepPanel1(src, comp);
        defXY1 = pickDefaultsPanel1(level);

        excl1 = idColsToHide();
        [x1, y1, g1] = pickXYG(tbl1, snap.x, snap.y, snap.g, defXY1, excl1);

        delete(allchild(state.hPnlLeft));
        tblGUI_scatHist(tbl1, ...
            'Parent', state.hPnlLeft, ...
            'xVar', x1, 'yVar', y1, 'grpVar', g1, ...
            'xScale', snap.xs, 'yScale', snap.ys, 'fitType', snap.ft, ...
            'clr', state.clr, ...
            'varsExclude', excl1);
        n1 = height(tbl1);
    elseif ~isempty(state.hPnlLeft.UserData) && isstruct(state.hPnlLeft.UserData) ...
            && isfield(state.hPnlLeft.UserData, 'tbl')
        n1 = height(state.hPnlLeft.UserData.tbl);
    end

    % --- Panel 2 ---
    n2 = NaN;
    if any(strcmpi(which, {'all', 'panel2'}))
        snap = readPanelState(state.hPnlRight);
        tbl2 = prepPanel2(src, level, comp, metric);
        defXY2 = struct('x', 'cyto', 'y', 'mito');

        excl2 = idColsToHide();
        [x2, y2, g2] = pickXYG(tbl2, snap.x, snap.y, snap.g, defXY2, excl2);

        delete(allchild(state.hPnlRight));
        tblGUI_scatHist(tbl2, ...
            'Parent', state.hPnlRight, ...
            'xVar', x2, 'yVar', y2, 'grpVar', g2, ...
            'xScale', snap.xs, 'yScale', snap.ys, 'fitType', snap.ft, ...
            'clr', state.clr, ...
            'varsExclude', excl2);
        n2 = height(tbl2);
    elseif ~isempty(state.hPnlRight.UserData) && isstruct(state.hPnlRight.UserData) ...
            && isfield(state.hPnlRight.UserData, 'tbl')
        n2 = height(state.hPnlRight.UserData.tbl);
    end

    set(state.hStatus, 'String', sprintf( ...
        'Level: %s | Comp: %s | Metric: %s | n1=%d  n2=%d', ...
        level, comp, metric, n1, n2));

    hFig.UserData = state;
end

function populateMetricList(hFig, src)
% Refresh ddMetric's options from the source table's symmetric metrics.
% Preserves the current selection if still available, else falls back to
% pendingMetric (from input args, used on first render), 'amp', or the
% first available metric.
    state = hFig.UserData;
    metrics = listSymmetricMetrics(src, state.asymPrefix, state.asymExtra);
    if isempty(metrics), metrics = {' '}; end

    current = '';
    items = get(state.ddMetric, 'String');
    if iscell(items) && ~isempty(items) && get(state.ddMetric, 'Value') <= numel(items)
        current = items{get(state.ddMetric, 'Value')};
    end

    pick = '';
    if ismember(current, metrics)
        pick = current;
    elseif isfield(state, 'pendingMetric') && ~isempty(state.pendingMetric) ...
            && ismember(state.pendingMetric, metrics)
        pick = state.pendingMetric;
    elseif ismember('amp', metrics)
        pick = 'amp';
    else
        pick = metrics{1};
    end

    % Atomic set so the intermediate state (new String, old out-of-range
    % Value) never reaches MATLAB's validator.
    set(state.ddMetric, 'String', metrics, ...
        'Value', find(strcmp(metrics, pick), 1));

    % pendingMetric only matters on the very first render.
    state.pendingMetric = '';
    hFig.UserData = state;
end

%% ========================================================================
%  HELPERS: STATE SNAPSHOT, DEFAULTS, EXCLUSIONS
%  ========================================================================

function snap = readPanelState(hPanel)
% Read the embedded tblGUI_scatHist's current selections (X/Y/group, X
% scale, Y scale, fit type). Returns empty fields where unavailable.
    snap = struct('x', '', 'y', '', 'g', '', 'xs', '', 'ys', '', 'ft', '');
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
        snap = struct('x', '', 'y', '', 'g', '', 'xs', '', 'ys', '', 'ft', '');
    end
end

function [x, y, g] = pickXYG(tbl, oldX, oldY, oldG, defXY, excl)
% Validate carried selections against the new schema; fall back to
% per-panel defaults, then let tblGUI_scatHist pick if even those fail.
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

function defXY = pickDefaultsPanel1(level)
% Within-compartment canonical scatter: amp vs the partner metric.
    defXY = struct();
    defXY.x = 'amp';
    if strcmpi(level, 'event')
        defXY.y = 'pairAmp';
    else
        defXY.y = 'pairFlux';
    end
end

function excl = idColsToHide()
% Columns that should never appear in the inner widget's X/Y or group
% dropdowns - identifiers, time indices, and the raw trace blob.
    excl = {'sbjID', 'unitID', 'excluded', 'trace', ...
            'start', 'stop', 'pairIdx'};
end

%% ========================================================================
%  HELPERS: METRICS + TABLE PREPARATION
%  ========================================================================

function metrics = listSymmetricMetrics(src, asymPrefix, asymExtra)
% Numeric columns of the source table whose meaning is symmetric across
% compartments (suitable for cyto / mito pairing in Panel 2).
    cols  = src.Properties.VariableNames;
    isNum = varfun(@isnumeric, src, 'OutputFormat', 'uniform');
    skip  = [idColsToHide(), {'compartment'}];

    metrics = cell(0, 1);
    for iC = 1:numel(cols)
        c = cols{iC};
        if ~isNum(iC),                              continue; end
        if ismember(c, skip),                       continue; end
        if isAsym(c, asymPrefix, asymExtra),        continue; end
        metrics{end+1, 1} = c; %#ok<AGROW>
    end
end

function tbl1 = prepPanel1(src, comp)
% Filter source table to the requested compartment. No reshape; the
% inner widget sees the table as-is.
    if strcmpi(comp, 'cyto')
        keep = src.compartment == 'Cyto';
    else
        keep = src.compartment == 'Mito';
    end
    tbl1 = src(keep, :);
end

function tbl2 = prepPanel2(src, level, comp, metric)
% Wide reshape - one row per paired (cyto, mito) observation. Only two
% metric columns ('cyto', 'mito'), both populated from src.(metric) on
% the appropriate side of each pair.

    if strcmpi(level, 'event')
        % Focal compartment is the selected one; partner via pairIdx.
        if strcmpi(comp, 'cyto')
            focalCmp = 'Cyto';
        else
            focalCmp = 'Mito';
        end

        isFocal     = src.compartment == focalCmp & ~isnan(src.pairIdx);
        focalRows   = find(isFocal);
        partnerRows = src.pairIdx(isFocal);

        if strcmpi(focalCmp, 'Cyto')
            cytoIdx = focalRows;  mitoIdx = partnerRows;
        else
            cytoIdx = partnerRows; mitoIdx = focalRows;
        end

        tbl2 = table();
        tbl2.sbjID    = src.sbjID(focalRows);
        tbl2.genotype = src.genotype(focalRows);
        tbl2.cyto     = src.(metric)(cytoIdx);
        tbl2.mito     = src.(metric)(mitoIdx);

    else
        % Cell level: pair Cyto and Mito rows of tblCell by sbjID.
        isC = src.compartment == 'Cyto';
        isM = src.compartment == 'Mito';
        cytoTbl = src(isC, :);
        mitoTbl = src(isM, :);

        [~, idxM] = ismember(cytoTbl.sbjID, mitoTbl.sbjID);
        keep      = idxM > 0;
        cytoTbl   = cytoTbl(keep, :);
        mitoTbl   = mitoTbl(idxM(keep), :);

        tbl2 = table();
        tbl2.sbjID    = cytoTbl.sbjID;
        tbl2.genotype = cytoTbl.genotype;
        tbl2.cyto     = cytoTbl.(metric);
        tbl2.mito     = mitoTbl.(metric);
    end
end

function flag = isAsym(name, prefix, explicit)
% A column is asymmetric (excluded from Panel 2's metric list) if its
% name starts with the configured prefix OR appears in the explicit
% list. Both are exposed as function parameters so the naming convention
% can evolve without editing this file.
    flag = false;
    if ~isempty(prefix) && startsWith(name, prefix)
        flag = true; return;
    end
    if ~isempty(explicit) && ismember(name, explicit)
        flag = true;
    end
end
