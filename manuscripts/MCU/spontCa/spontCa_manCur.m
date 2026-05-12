function hFig = spontCa_manCur(tbl, fs, varargin)
% SPONTCA_MANCUR Per-cell manual curation GUI for SpontCa events.
%
%   hFig = SPONTCA_MANCUR(TBL, FS, ...) opens an interactive editor for the
%   per-cell event lists in the long-format table TBL (output of
%   SPONTCA_LOAD + the detection loop). The autodetector saturates against
%   plateau noise at fs=3 Hz; this GUI produces a human-curated gold
%   standard by allowing direct drag-edit / add / delete of events on top
%   of the autodetector pre-fill.
%
%   LAYOUT (2x2 + side panel):
%       cyto FULL trace            | cyto ZOOM
%       mito FULL trace            | mito ZOOM
%   Within columns, the two compartments share x. The full trace shows a
%   translucent rectangle marking what the zoom column displays; clicking
%   on a full trace recenters the zoom, dragging the rectangle edges
%   resizes the zoom window.
%
%   INTERACTIONS:
%       FULL trace, left-click           -> recenter zoom on click
%       FULL trace, drag rect edge       -> resize zoom width
%       ZOOM, drag dashed line           -> move event peak (snap to sample)
%       ZOOM, drag dotted line           -> move event stop
%       ZOOM, left-click empty           -> ADD event (stop via flat-d walk)
%       ZOOM, right-click line           -> DELETE event
%   KEYBOARD:
%       LEFT / RIGHT   prev / next event in active compartment
%       UP / DOWN      switch active compartment (cyto <-> mito)
%
%   PERSISTENCE:
%       Curated lists are written to spontCa_curated/ (alongside this
%       file) as <sbjID>_<compartment>.mat. On open, the curated file (if any)
%       overrides the auto-detection pre-fill in TBL. Switching cells or
%       hitting Save persists; switching with unsaved changes auto-saves.
%       Both compartments are written on each save. "Reset to auto"
%       restores tbl.start/stop/... for the current cell and re-saves.
%
%   INPUTS:
%       tbl - long-format table; two rows per sbjID (Cyto and Mito).
%             Required cols: sbjID, compartment, trace (matrix), start,
%             stop, amp, dur, int (cells of doubles).
%       fs  - sampling rate (Hz).
%
%   OPTIONAL (Name-Value):
%       'sbjID' - (char) initial cell to display.
%
%   See also: SPONTCA_GUI, SPONTCA_DETECT, SPONTCA_LOAD

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addRequired(p, 'fs',  @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'sbjID', '', @(x) ischar(x) || isstring(x) || isempty(x));
parse(p, tbl, fs, varargin{:});
initSbj = char(p.Results.sbjID);

cfg     = mcu_cfg;
clrCmp  = cfg.clr.cmp;      % row 1 = cyto, row 2 = mito
clrCyto = clrCmp(1, :);
clrMito = clrCmp(2, :);
clrRect = [0.15, 0.45, 0.85];
alphaLine = 0.65;

nT = size(tbl.trace, 2);
t  = (0:nT-1) / fs;
dt = 1 / fs;

% Local-baseline constants mirrored from spontCa_detect so amp/int
% recomputations on edit match the autodetector's conventions.
detCfg.bslWin   = 30;
detCfg.quantBsl = 20;

curDir = fullfile(fileparts(mfilename('fullpath')), 'spontCa_curated');
if ~exist(curDir, 'dir')
    mkdir(curDir);
end


%% ========================================================================
%  FIGURE LAYOUT
%  ========================================================================

hFig = figure('Name', 'SpontCa ManCur', 'NumberTitle', 'off', 'Color', 'w', ...
    'Units', 'pixels', 'Position', [80, 60, 1500, 900]);

pW = 0.11;
hSide = uipanel('Parent', hFig, 'Units', 'normalized', ...
    'Position', [0, 0, pW, 1], 'BorderType', 'none');
hMain = uipanel('Parent', hFig, 'Units', 'normalized', ...
    'Position', [pW, 0, 1 - pW, 1], 'BorderType', 'none');

margL = 0.06;  margR = 0.02;
margT = 0.04;  margB = 0.08;
gapV  = 0.0;
gapH  = 0.03;

availW = 1 - margL - margR;
availH = 1 - margT - margB;
rowH   = (availH - gapV) / 2;
colW   = (availW - gapH) / 2;

yMito = margB;
yCyto = yMito + rowH + gapV;
xFull = margL;
xZoom = margL + colW + gapH;

axCytoFull = axes('Parent', hMain, 'Position', [xFull, yCyto, colW, rowH]);
axCytoZoom = axes('Parent', hMain, 'Position', [xZoom, yCyto, colW, rowH]);
axMitoFull = axes('Parent', hMain, 'Position', [xFull, yMito, colW, rowH]);
axMitoZoom = axes('Parent', hMain, 'Position', [xZoom, yMito, colW, rowH]);

linkaxes([axCytoFull, axMitoFull], 'x');
linkaxes([axCytoZoom, axMitoZoom], 'x');


%% ========================================================================
%  SIDE PANEL CONTROLS
%  ========================================================================

uicontrol('Parent', hSide, 'Style', 'text', 'String', 'Cell:', ...
    'Units', 'normalized', 'Position', [0.05, 0.94, 0.9, 0.04], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');

cellList = cellstr(string(tbl.sbjID(tbl.compartment == 'Cyto')));
cellList = sort(cellList);

ddCell = uicontrol('Parent', hSide, 'Style', 'popupmenu', ...
    'String', cellList, ...
    'Units', 'normalized', 'Position', [0.05, 0.90, 0.9, 0.04], ...
    'Callback', @onCellChange);

if ~isempty(initSbj)
    idx = find(strcmp(cellList, initSbj), 1);
    if ~isempty(idx), ddCell.Value = idx; end
end

uicontrol('Parent', hSide, 'Style', 'text', 'String', 'Active:', ...
    'Units', 'normalized', 'Position', [0.05, 0.83, 0.9, 0.03], ...
    'HorizontalAlignment', 'left');
lblActive = uicontrol('Parent', hSide, 'Style', 'text', 'String', 'Cyto', ...
    'Units', 'normalized', 'Position', [0.05, 0.79, 0.9, 0.04], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold', ...
    'ForegroundColor', clrCyto);

uicontrol('Parent', hSide, 'Style', 'text', 'String', 'Zoom width (s):', ...
    'Units', 'normalized', 'Position', [0.05, 0.72, 0.9, 0.03], ...
    'HorizontalAlignment', 'left');
edZoomW = uicontrol('Parent', hSide, 'Style', 'edit', 'String', '20', ...
    'Units', 'normalized', 'Position', [0.05, 0.68, 0.9, 0.04], ...
    'Callback', @onZoomWidthEdit);

uicontrol('Parent', hSide, 'Style', 'text', 'String', 'Zoom center (s):', ...
    'Units', 'normalized', 'Position', [0.05, 0.63, 0.9, 0.03], ...
    'HorizontalAlignment', 'left');
edZoomCenter = uicontrol('Parent', hSide, 'Style', 'edit', 'String', '0', ...
    'Units', 'normalized', 'Position', [0.05, 0.59, 0.9, 0.04], ...
    'Callback', @onZoomCenterEdit);

uicontrol('Parent', hSide, 'Style', 'pushbutton', 'String', 'Save', ...
    'Units', 'normalized', 'Position', [0.05, 0.49, 0.9, 0.05], ...
    'Callback', @(~,~) saveCurrent(true));

lblDirty = uicontrol('Parent', hSide, 'Style', 'text', 'String', 'saved', ...
    'Units', 'normalized', 'Position', [0.05, 0.44, 0.9, 0.04], ...
    'HorizontalAlignment', 'center', 'ForegroundColor', [0.2, 0.6, 0.2]);

uicontrol('Parent', hSide, 'Style', 'pushbutton', 'String', 'Load auto', ...
    'Units', 'normalized', 'Position', [0.05, 0.37, 0.9, 0.05], ...
    'Callback', @onLoadAuto);

uicontrol('Parent', hSide, 'Style', 'pushbutton', 'String', 'Load man', ...
    'Units', 'normalized', 'Position', [0.05, 0.31, 0.9, 0.05], ...
    'Callback', @onLoadMan);

uicontrol('Parent', hSide, 'Style', 'text', 'String', ...
    sprintf(['Drag dashed: peak\n' ...
             'Drag dotted: stop\n' ...
             'L-click empty: add\n' ...
             'R-click line: delete\n' ...
             'Arrows L/R: window\n' ...
             'Arrows U/D: comp']), ...
    'Units', 'normalized', 'Position', [0.05, 0.04, 0.9, 0.24], ...
    'HorizontalAlignment', 'left', 'FontAngle', 'italic', ...
    'ForegroundColor', [0.4, 0.4, 0.4]);


%% ========================================================================
%  STATE
%  ========================================================================
% Held in figure UserData so callbacks can mutate it. S.events.(Cyto|Mito)
% is a struct of column vectors. S.bsl is the rolling 20th-pct baseline
% used for amp/int recomputation on edits, computed once per cell load.

S.cellList    = cellList;
S.currentCell = '';
S.activeCmp   = 'Cyto';
S.events      = struct('Cyto', emptyEv(), 'Mito', emptyEv());
S.traces      = struct('Cyto', [], 'Mito', []);
S.bsl         = struct('Cyto', [], 'Mito', []);
S.zoomCenter  = 0;
S.zoomWidth   = 20;
S.dirty       = false;
S.drag        = [];
hFig.UserData = S;

set(hFig, 'WindowButtonDownFcn',   @onMouseDown);
set(hFig, 'WindowButtonUpFcn',     @onMouseUp);
set(hFig, 'WindowButtonMotionFcn', @onMouseMove);
set(hFig, 'KeyPressFcn',           @onKeyPress);

onCellChange();


%% ========================================================================
%  CELL LOAD / SAVE
%  ========================================================================

    function onCellChange(~, ~)
        S = hFig.UserData;
        items = get(ddCell, 'String');
        if isempty(items), return; end
        sName = items{get(ddCell, 'Value')};

        if S.dirty && ~isempty(S.currentCell)
            saveCurrent(false);
            S = hFig.UserData;
        end

        iC = find(tbl.sbjID == sName & tbl.compartment == 'Cyto');
        iM = find(tbl.sbjID == sName & tbl.compartment == 'Mito');
        if isempty(iC) || isempty(iM), return; end

        S.currentCell = char(sName);
        S.traces.Cyto = tbl.trace(iC, :);
        S.traces.Mito = tbl.trace(iM, :);
        S.bsl.Cyto = rollingPercentileLocal(S.traces.Cyto, ...
            round(detCfg.bslWin * fs), detCfg.quantBsl);
        S.bsl.Mito = rollingPercentileLocal(S.traces.Mito, ...
            round(detCfg.bslWin * fs), detCfg.quantBsl);

        S.events.Cyto = loadEvForRow(S.currentCell, 'Cyto', iC);
        S.events.Mito = loadEvForRow(S.currentCell, 'Mito', iM);

        evA = S.events.(S.activeCmp);
        if ~isempty(evA.start)
            S.zoomCenter = evA.start(1);
        else
            S.zoomCenter = t(round(end / 2));
        end

        S.dirty = false;
        hFig.UserData = S;
        setDirty(false);
        redrawAll();
    end


    function ev = loadEvForRow(sName, cmp, iRow)
        fpath = fullfile(curDir, sprintf('%s_%s.mat', sName, cmp));
        if exist(fpath, 'file')
            load(fpath, 'cur');
            ev = struct('start', cur.start(:), 'stop', cur.stop(:), ...
                'amp', cur.amp(:), 'dur', cur.dur(:), 'int', cur.int(:));
        else
            ev = struct( ...
                'start', tbl.start{iRow}(:), ...
                'stop',  tbl.stop{iRow}(:), ...
                'amp',   tbl.amp{iRow}(:), ...
                'dur',   tbl.dur{iRow}(:), ...
                'int',   tbl.int{iRow}(:));
        end
        ev = sortEv(ev);
    end


    function saveCurrent(verbose)
        S = hFig.UserData;
        if isempty(S.currentCell), return; end
        bkupDir = fullfile(curDir, 'bkup');
        stamp   = datestr(now, 'yymmdd_HHMMSS'); %#ok<TNOW1,DATST>
        cmps = {'Cyto', 'Mito'};
        for k = 1:2
            cmp = cmps{k};
            ev = S.events.(cmp);
            fpath = fullfile(curDir, sprintf('%s_%s.mat', S.currentCell, cmp));
            % Back up any prior version before overwrite. Each save makes a
            % stamped copy so curation is recoverable across sessions.
            if exist(fpath, 'file')
                if ~exist(bkupDir, 'dir'), mkdir(bkupDir); end
                bkupName = sprintf('%s_%s_%s.mat', ...
                    S.currentCell, cmp, stamp);
                copyfile(fpath, fullfile(bkupDir, bkupName));
            end
            cur = struct( ...
                'sbjID',       S.currentCell, ...
                'compartment', cmp, ...
                'fs',          fs, ...
                'start',       ev.start(:), ...
                'stop',        ev.stop(:), ...
                'amp',         ev.amp(:), ...
                'dur',         ev.dur(:), ...
                'int',         ev.int(:), ...
                'savedAt',     datestr(now, 'yyyy-mm-dd HH:MM:SS')); %#ok<TNOW1,DATST>
            save(fpath, 'cur');
        end
        S.dirty = false;
        hFig.UserData = S;
        setDirty(false);
        if verbose
            fprintf('Saved %s (cyto + mito) to %s\n', S.currentCell, curDir);
        end
    end


    function onLoadAuto(~, ~)
        % Restore the autodetector pre-fill for the current cell.
        % Marks dirty (caller can Save to persist).
        S = hFig.UserData;
        if isempty(S.currentCell), return; end
        sName = S.currentCell;
        iC = find(tbl.sbjID == sName & tbl.compartment == 'Cyto');
        iM = find(tbl.sbjID == sName & tbl.compartment == 'Mito');
        S.events.Cyto = sortEv(struct( ...
            'start', tbl.start{iC}(:), 'stop', tbl.stop{iC}(:), ...
            'amp',   tbl.amp{iC}(:),   'dur',  tbl.dur{iC}(:), ...
            'int',   tbl.int{iC}(:)));
        S.events.Mito = sortEv(struct( ...
            'start', tbl.start{iM}(:), 'stop', tbl.stop{iM}(:), ...
            'amp',   tbl.amp{iM}(:),   'dur',  tbl.dur{iM}(:), ...
            'int',   tbl.int{iM}(:)));
        S.dirty = true;
        hFig.UserData = S;
        setDirty(true);
        redrawAll();
    end


    function onLoadMan(~, ~)
        % Load the most recent manual save for the current cell. Per
        % compartment, prefer spontCa_curated/<sbjID>_<cmp>.mat; fall back
        % to the newest spontCa_curated/bkup/<sbjID>_<cmp>_*.mat by file
        % mtime. If neither exists for a compartment, leave it untouched.
        S = hFig.UserData;
        if isempty(S.currentCell), return; end
        sName = S.currentCell;
        bkupDir = fullfile(curDir, 'bkup');
        cmps = {'Cyto', 'Mito'};
        nLoaded = 0;
        for k = 1:2
            cmp = cmps{k};
            fpath = fullfile(curDir, sprintf('%s_%s.mat', sName, cmp));
            src = '';
            if exist(fpath, 'file')
                src = fpath;
            elseif exist(bkupDir, 'dir')
                files = dir(fullfile(bkupDir, ...
                    sprintf('%s_%s_*.mat', sName, cmp)));
                if ~isempty(files)
                    [~, idx] = max([files.datenum]);
                    src = fullfile(bkupDir, files(idx).name);
                end
            end
            if isempty(src)
                fprintf('No manual save found for %s %s\n', sName, cmp);
                continue;
            end
            L = load(src, 'cur');
            if ~isfield(L, 'cur')
                fprintf('Skipped %s: no ''cur'' struct\n', src);
                continue;
            end
            cur = L.cur;
            S.events.(cmp) = sortEv(struct( ...
                'start', cur.start(:), 'stop', cur.stop(:), ...
                'amp',   cur.amp(:),   'dur',  cur.dur(:), ...
                'int',   cur.int(:)));
            nLoaded = nLoaded + 1;
            fprintf('Loaded %s %s from %s\n', sName, cmp, src);
        end
        if nLoaded > 0
            S.dirty = true;
            hFig.UserData = S;
            setDirty(true);
            redrawAll();
        end
    end


%% ========================================================================
%  REDRAW
%  ========================================================================

    function redrawAll()
        S = hFig.UserData;
        if isempty(S.currentCell), return; end

        set(hFig, 'Name', sprintf( ...
            'SpontCa ManCur - %s | cyto n=%d | mito n=%d', ...
            S.currentCell, ...
            numel(S.events.Cyto.start), numel(S.events.Mito.start)));

        drawFull(axCytoFull, S.traces.Cyto, S.events.Cyto, clrCyto);
        drawFull(axMitoFull, S.traces.Mito, S.events.Mito, clrMito);
        drawZoom(axCytoZoom, S.traces.Cyto, S.events.Cyto, clrCyto);
        drawZoom(axMitoZoom, S.traces.Mito, S.events.Mito, clrMito);

        applyZoomLims();
        drawZoomIndicator();
        updateZoomCenterDisplay();
        highlightActive();
    end


    function drawFull(ax, trace, ev, clr)
        cla(ax, 'reset');
        hold(ax, 'on');
        plot(ax, t, trace, 'Color', clr, 'LineWidth', 0.7, ...
            'HandleVisibility', 'off', 'PickableParts', 'none');
        axis(ax, 'tight');
        yL = ylim(ax);
        plotMarks(ax, ev.start, yL, clr, '--', alphaLine);
        if isMitoAx(ax)
            plotMarks(ax, ev.stop, yL, clr, ':',  alphaLine);
        end
        applyAxisLabels(ax);
        hold(ax, 'off');
    end


    function drawZoom(ax, trace, ev, clr)
        cla(ax, 'reset');
        hold(ax, 'on');
        plot(ax, t, trace, 'Color', clr, 'LineWidth', 0.9, ...
            'HandleVisibility', 'off', 'PickableParts', 'none', ...
            'Tag', 'zoomTrace');
        axis(ax, 'tight');
        yL = ylim(ax);
        % Each event line is a distinct Line object tagged with role+idx so
        % the mouse handlers can hit-test and drag it directly. Cyto stop
        % lines are suppressed (cyto stops are no longer biologically
        % meaningful at fs=3); start lines remain editable.
        drawStops = isMitoAx(ax);
        for iE = 1:numel(ev.start)
            line(ax, [ev.start(iE), ev.start(iE)], yL, ...
                'Color', [clr, alphaLine], 'LineStyle', '--', ...
                'LineWidth', 1.2, 'Tag', sprintf('startLine_%d', iE), ...
                'UserData', struct('role', 'start', 'idx', iE), ...
                'PickableParts', 'visible');
            if drawStops
                line(ax, [ev.stop(iE), ev.stop(iE)], yL, ...
                    'Color', [clr, alphaLine], 'LineStyle', ':', ...
                    'LineWidth', 1.2, 'Tag', sprintf('stopLine_%d', iE), ...
                    'UserData', struct('role', 'stop', 'idx', iE), ...
                    'PickableParts', 'visible');
            end
        end
        applyAxisLabels(ax);
        hold(ax, 'off');
    end


    function tf = isMitoAx(ax)
        tf = (ax == axMitoFull) || (ax == axMitoZoom);
    end


    function applyAxisLabels(ax)
        % Only the mito (bottom) row carries the x-label and tick labels;
        % the cyto (top) row sits flush above it with x-axis linked.
        if ax == axMitoFull || ax == axMitoZoom
            xlabel(ax, 'Time (s)');
        else
            set(ax, 'XTickLabel', []);
        end
    end


    function plotMarks(ax, x, yL, clr, ls, a)
        if isempty(x), return; end
        x  = x(:);
        xx = reshape([x, x, nan(length(x), 1)]', [], 1);
        yy = repmat([yL(1); yL(2); NaN], length(x), 1);
        plot(ax, xx, yy, 'Color', [clr, a], 'LineStyle', ls, ...
            'LineWidth', 0.7, 'HandleVisibility', 'off', ...
            'PickableParts', 'none');
    end


    function applyZoomLims()
        S = hFig.UserData;
        hw = S.zoomWidth / 2;
        xl = [max(0, S.zoomCenter - hw), min(t(end), S.zoomCenter + hw)];
        if diff(xl) <= 0
            xl = [0, min(S.zoomWidth, t(end))];
        end
        xlim(axCytoZoom, xl);
    end


    function drawZoomIndicator()
        % findall (not findobj) so HandleVisibility=off patches get deleted.
        % Without this, every redraw leaves a stale patch behind and the
        % rectangle accumulates into a gradient.
        xl = xlim(axCytoZoom);
        for ax = [axCytoFull, axMitoFull]
            delete(findall(ax, 'Tag', 'zoomRect'));
            delete(findall(ax, 'Tag', 'zoomEdgeL'));
            delete(findall(ax, 'Tag', 'zoomEdgeR'));
            yL = ylim(ax);
            patch(ax, [xl(1), xl(2), xl(2), xl(1)], ...
                  [yL(1), yL(1), yL(2), yL(2)], clrRect, ...
                  'FaceAlpha', 0.06, 'EdgeColor', 'none', ...
                  'Tag', 'zoomRect', 'PickableParts', 'none');
            line(ax, [xl(1), xl(1)], yL, 'Color', [clrRect, 0.3], ...
                'LineWidth', 1.0, 'Tag', 'zoomEdgeL', ...
                'PickableParts', 'visible');
            line(ax, [xl(2), xl(2)], yL, 'Color', [clrRect, 0.3], ...
                'LineWidth', 1.0, 'Tag', 'zoomEdgeR', ...
                'PickableParts', 'visible');
        end
    end


    function highlightActive()
        S = hFig.UserData;
        if strcmp(S.activeCmp, 'Cyto')
            set(lblActive, 'String', 'Cyto', 'ForegroundColor', clrCyto);
            set([axCytoFull, axCytoZoom], 'LineWidth', 1.8, ...
                'XColor', [0 0 0], 'YColor', [0 0 0]);
            set([axMitoFull, axMitoZoom], 'LineWidth', 0.5, ...
                'XColor', [0.4 0.4 0.4], 'YColor', [0.4 0.4 0.4]);
        else
            set(lblActive, 'String', 'Mito', 'ForegroundColor', clrMito);
            set([axMitoFull, axMitoZoom], 'LineWidth', 1.8, ...
                'XColor', [0 0 0], 'YColor', [0 0 0]);
            set([axCytoFull, axCytoZoom], 'LineWidth', 0.5, ...
                'XColor', [0.4 0.4 0.4], 'YColor', [0.4 0.4 0.4]);
        end
    end


    function setDirty(flag)
        S = hFig.UserData;
        S.dirty = flag;
        hFig.UserData = S;
        if flag
            set(lblDirty, 'String', 'modified', ...
                'ForegroundColor', [0.85, 0.35, 0.15]);
        else
            set(lblDirty, 'String', 'saved', ...
                'ForegroundColor', [0.2, 0.6, 0.2]);
        end
    end


%% ========================================================================
%  MOUSE HANDLERS
%  ========================================================================

    function onMouseDown(~, ~)
        S = hFig.UserData;
        ax = gca;
        if ~ismember(ax, [axCytoFull, axCytoZoom, axMitoFull, axMitoZoom])
            return;
        end
        cp = ax.CurrentPoint;
        xClick = cp(1, 1);
        clickType = get(hFig, 'SelectionType');

        isFull = (ax == axCytoFull) || (ax == axMitoFull);
        cmpAx  = compartmentOfAx(ax);

        if isFull
            edge = hitTestZoomEdge(ax, xClick);
            if ~isempty(edge)
                S.drag = struct('mode', 'zoomEdge', 'edge', edge, 'ax', ax);
                hFig.UserData = S;
                return;
            end
            S.zoomCenter = xClick;
            hFig.UserData = S;
            applyZoomLims();
            drawZoomIndicator();
            updateZoomCenterDisplay();
            return;
        end

        % ZOOM panel: try event-line hit first.
        [hit, role, iE] = hitTestEventLine(ax, xClick);

        if ~isempty(hit) && strcmp(clickType, 'alt')
            S.events.(cmpAx) = deleteEvent(S.events.(cmpAx), iE);
            S.activeCmp = cmpAx;
            % Clear any in-flight drag - redrawAll about to wipe handles.
            S.drag = [];
            hFig.UserData = S;
            setDirty(true);
            redrawAll();
            return;
        end

        if ~isempty(hit) && strcmp(clickType, 'normal')
            S.drag = struct('mode', 'evLine', 'role', role, 'idx', iE, ...
                'ax', ax, 'cmp', cmpAx, 'hLine', hit);
            S.activeCmp = cmpAx;
            hFig.UserData = S;
            highlightActive();
            return;
        end

        if strcmp(clickType, 'normal')
            tNew = snapTime(xClick, fs);
            S.events.(cmpAx) = addEventLocal(S.events.(cmpAx), tNew, ...
                S.traces.(cmpAx), S.bsl.(cmpAx));
            S.activeCmp = cmpAx;
            hFig.UserData = S;
            setDirty(true);
            redrawAll();
        end
    end


    function onMouseMove(~, ~)
        S = hFig.UserData;
        if isempty(S.drag), return; end
        ax = S.drag.ax;
        cp = ax.CurrentPoint;
        xNow = cp(1, 1);

        switch S.drag.mode
            case 'zoomEdge'
                xl = xlim(axCytoZoom);
                if strcmp(S.drag.edge, 'L')
                    xl(1) = min(xNow, xl(2) - dt);
                else
                    xl(2) = max(xNow, xl(1) + dt);
                end
                xl(1) = max(0, xl(1));
                xl(2) = min(t(end), xl(2));
                S.zoomCenter = mean(xl);
                S.zoomWidth  = diff(xl);
                hFig.UserData = S;
                set(edZoomW, 'String', sprintf('%.2f', S.zoomWidth));
                applyZoomLims();
                drawZoomIndicator();
                updateZoomCenterDisplay();

            case 'evLine'
                % Preview: move the line only; commit on mouse-up.
                % isvalid guard: handle can go stale if a redrawAll
                % fired between mouse-down and mouse-move (e.g. a
                % delete from a right-click landed in between).
                if ~isvalid(S.drag.hLine)
                    S.drag = [];
                    hFig.UserData = S;
                    return;
                end
                xSnap = snapTime(xNow, fs);
                set(S.drag.hLine, 'XData', [xSnap, xSnap]);
        end
    end


    function onMouseUp(~, ~)
        S = hFig.UserData;
        if isempty(S.drag), return; end
        switch S.drag.mode
            case 'evLine'
                xData = get(S.drag.hLine, 'XData');
                tNew  = snapTime(xData(1), fs);
                cmp   = S.drag.cmp;
                ev    = S.events.(cmp);
                tr    = S.traces.(cmp);
                bs    = S.bsl.(cmp);
                iE    = S.drag.idx;
                if strcmp(S.drag.role, 'start')
                    tNew = min(tNew, ev.stop(iE) - dt);
                    tNew = max(tNew, 0);
                    ev.start(iE) = tNew;
                else
                    tNew = max(tNew, ev.start(iE) + dt);
                    tNew = min(tNew, t(end));
                    ev.stop(iE) = tNew;
                end
                ev = recomputeOneLocal(ev, iE, tr, bs);
                ev = sortEv(ev);
                S.events.(cmp) = ev;
                S.drag = [];
                hFig.UserData = S;
                setDirty(true);
                redrawAll();
            otherwise
                S.drag = [];
                hFig.UserData = S;
        end
    end


    function edge = hitTestZoomEdge(ax, xClick)
        % 1.5% of visible x-range counts as an edge grab.
        xl = xlim(axCytoZoom);
        ax_xl = xlim(ax);
        tol = 0.015 * diff(ax_xl);
        if abs(xClick - xl(1)) < tol
            edge = 'L';
        elseif abs(xClick - xl(2)) < tol
            edge = 'R';
        else
            edge = '';
        end
    end


    function [hit, role, iE] = hitTestEventLine(ax, xClick)
        ax_xl = xlim(ax);
        tol = 0.012 * diff(ax_xl);
        lines = findobj(ax, '-regexp', 'Tag', '^(start|stop)Line_');
        hit = []; role = ''; iE = NaN;
        bestD = inf;
        for k = 1:numel(lines)
            xd = get(lines(k), 'XData');
            d  = abs(xd(1) - xClick);
            if d < tol && d < bestD
                bestD = d;
                hit = lines(k);
                ud  = get(lines(k), 'UserData');
                role = ud.role;
                iE   = ud.idx;
            end
        end
    end


    function cmp = compartmentOfAx(ax)
        if ax == axCytoFull || ax == axCytoZoom
            cmp = 'Cyto';
        else
            cmp = 'Mito';
        end
    end


    function onZoomWidthEdit(src, ~)
        v = str2double(get(src, 'String'));
        S = hFig.UserData;
        if ~isfinite(v) || v <= 0
            set(src, 'String', sprintf('%.2f', S.zoomWidth));
            return;
        end
        S.zoomWidth = min(v, t(end));
        hFig.UserData = S;
        applyZoomLims();
        drawZoomIndicator();
    end


    function onZoomCenterEdit(src, ~)
        v = str2double(get(src, 'String'));
        S = hFig.UserData;
        if ~isfinite(v)
            set(src, 'String', sprintf('%.2f', S.zoomCenter));
            return;
        end
        S.zoomCenter = min(max(0, v), t(end));
        hFig.UserData = S;
        applyZoomLims();
        drawZoomIndicator();
        updateZoomCenterDisplay();
    end


    function updateZoomCenterDisplay()
        S = hFig.UserData;
        set(edZoomCenter, 'String', sprintf('%.2f', S.zoomCenter));
    end


%% ========================================================================
%  KEYBOARD
%  ========================================================================

    function onKeyPress(~, evt)
        S = hFig.UserData;
        if isempty(S.currentCell), return; end
        % Step keeps a 20% overlap with the previous window so events
        % straddling a boundary aren't missed on a fast scan.
        stepFrac = 0.8;
        switch evt.Key
            case 'rightarrow'
                S.zoomCenter = min(t(end), ...
                    S.zoomCenter + stepFrac * S.zoomWidth);
                hFig.UserData = S;
                applyZoomLims();
                drawZoomIndicator();
                updateZoomCenterDisplay();
            case 'leftarrow'
                S.zoomCenter = max(0, ...
                    S.zoomCenter - stepFrac * S.zoomWidth);
                hFig.UserData = S;
                applyZoomLims();
                drawZoomIndicator();
                updateZoomCenterDisplay();
            case {'uparrow', 'downarrow'}
                if strcmp(S.activeCmp, 'Cyto')
                    S.activeCmp = 'Mito';
                else
                    S.activeCmp = 'Cyto';
                end
                hFig.UserData = S;
                highlightActive();
        end
    end


%% ========================================================================
%  EVENT MUTATION (nested: needs fs, detCfg, dt in scope)
%  ========================================================================

    function ev = addEventLocal(ev, tStart, trace, bsl)
        % Default stop is two samples past the start: just enough to keep
        % the dotted line visible. The walk-forward heuristic was wrong
        % for clicks in non-flat regions (ran off the trace end). User
        % drags the stop afterwards if they want a wider event.
        pkSmp   = max(1, min(numel(trace), round(tStart * fs) + 1));
        stopSmp = min(numel(trace), pkSmp + 2);
        tStartS = (pkSmp  - 1) * dt;
        tStopS  = (stopSmp - 1) * dt;
        ampV = trace(pkSmp) - bsl(pkSmp);
        seg  = trace(pkSmp:stopSmp) - bsl(pkSmp:stopSmp);
        seg(isnan(seg)) = 0;
        intV = trapz(seg) * dt;
        durV = tStopS - tStartS;
        ev.start(end + 1, 1) = tStartS;
        ev.stop(end + 1, 1)  = tStopS;
        ev.amp(end + 1, 1)   = ampV;
        ev.dur(end + 1, 1)   = durV;
        ev.int(end + 1, 1)   = intV;
        ev = sortEv(ev);
    end


    function ev = recomputeOneLocal(ev, iE, trace, bsl)
        nTrc    = numel(trace);
        pkSmp   = max(1, min(nTrc, round(ev.start(iE) * fs) + 1));
        stopSmp = max(1, min(nTrc, round(ev.stop(iE)  * fs) + 1));
        stopSmp = max(stopSmp, pkSmp);
        ev.amp(iE) = trace(pkSmp) - bsl(pkSmp);
        seg = trace(pkSmp:stopSmp) - bsl(pkSmp:stopSmp);
        seg(isnan(seg)) = 0;
        ev.int(iE) = trapz(seg) * dt;
        ev.dur(iE) = ev.stop(iE) - ev.start(iE);
    end


end     % spontCa_manCur


%% ========================================================================
%  FILE-LOCAL HELPERS
%  ========================================================================

function tSnap = snapTime(t, fs)
tSnap = round(t * fs) / fs;
end


function ev = emptyEv()
ev = struct('start', zeros(0, 1), 'stop', zeros(0, 1), ...
            'amp',   zeros(0, 1), 'dur',  zeros(0, 1), ...
            'int',   zeros(0, 1));
end


function ev = sortEv(ev)
if isempty(ev.start), return; end
[~, ord] = sort(ev.start);
ev.start = ev.start(ord);
ev.stop  = ev.stop(ord);
ev.amp   = ev.amp(ord);
ev.dur   = ev.dur(ord);
ev.int   = ev.int(ord);
end


function ev = deleteEvent(ev, iE)
ev.start(iE) = [];
ev.stop(iE)  = [];
ev.amp(iE)   = [];
ev.dur(iE)   = [];
ev.int(iE)   = [];
end


function out = rollingPercentileLocal(x, winSmp, q)
% Centered, NaN-tolerant rolling percentile (matches spontCa_detect's
% internal rollingPercentile so baselines line up).
nT = length(x);
out = nan(1, nT);
halfWin = floor(winSmp / 2);
for i = 1:nT
    lo = max(1, i - halfWin);
    hi = min(nT, i + halfWin);
    seg = x(lo:hi);
    seg = seg(~isnan(seg));
    if ~isempty(seg)
        out(i) = prctile(seg, q);
    end
end
end
