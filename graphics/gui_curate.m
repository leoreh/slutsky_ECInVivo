function hFig = gui_curate(basepath, varargin)
% GUI_CURATE General-purpose event-curation viewer (EDs, ripples, ...).
%
%   hFig = GUI_CURATE(basepath, varargin)
%
%   SUMMARY:
%       Single-session signal viewer on uifigure, built on the +tblgui helper
%       layer. It is modality-agnostic: an "event" is just a peak time (with
%       optional start/stop and an accepted flag), and any 1-D signal can be a
%       panel. What differs between modalities (which file to load, which
%       signals to show, the default layout, where to save) is captured by a
%       PRESET. A Preset dropdown (EDs, Ripples, ...) loads the relevant file
%       and applies that preset's default view; switching presets reloads.
%       Presets live in curate_presets.m (one entry per modality).
%
%       Plot area: a stack of PANELS in two regions separated by a thin divider:
%           Top (wide)   - full-session overview (drawn once; x-range zoomable).
%           Bottom (narrow) - a window around the cursor t0 (redrawn as t0 moves).
%       All panels live in one tiledlayout, so every plot box shares one left
%       gutter and width. The Top prints its x-axis at the bottom of its lowest
%       panel (just above the divider) and the Bottom at the very bottom.
%
%       Panel TYPES (the draw function is chosen by the input type):
%           trace      - a 1-D signal (.data vector, .fs, optional .ylim, .clr)
%           spec       - a spectrogram adapter (.data = struct .s/.freq/.tstamps)
%           hypnogram  - sleep-state strip (.data = cell of [start end] in hours)
%           raster     - spike raster (.data = cell of spike-time vectors [s])
%           eventTicks - event marks + the event module (.data = struct with
%                        .peakTime, optional .times [N x 2], optional .accepted)
%
%       Navigation: a single cursor t0 drives the Bottom window. Clicking ANY
%       panel moves t0 to the clicked time; prev/next/accept/reject and the
%       event index set t0 to an event. The Top overview is fixed; a cursor line
%       and a shaded band on it mark t0 and the Bottom window.
%
%       Interaction:
%           Left / Right          previous / next event
%           Up / Down             accept / reject (auto-advances)
%           Ctrl+S                save (via the preset's save function)
%           + / -                 zoom the active region in / out
%           0                     reset the active region
%           click any panel       move t0 (the Bottom) to the clicked time
%           scroll over a region  zoom that region (by pointer position)
%
%   INPUTS:
%       basepath    - (Char) Session directory.
%       varargin    - Parameter/Value pairs:
%           'preset'   - (Char) Initial preset name (see curate_presets). If
%                        empty, auto-detected from the files present.
%           'inputs'   - (Struct array) Bypass presets with explicit inputs
%                        (each: .name .type [.data .fs .ylim .clr .label .height
%                        .defRegion .defOrder]). Pairs with 'panels'/'saveFcn'.
%           'panels'   - (Struct array) Each: .source .region (custom path).
%           'saveFcn'  - (Fcn) @(accepted) persistence for the custom path.
%           'winPlot'  - (Num) Bottom-window full width [s] (custom path). {1.0}
%           'basename' - (Char) Override (defaults to the folder name).
%           'Visible'  - (Char) 'on' (default) | 'off' (headless).
%
%   OUTPUT:
%       hFig        - (uifigure) Handle to the viewer window.
%
%   DEPENDENCIES:
%       curate_presets, tblgui.layout / labeledControl / eventPanel / notify /
%       chooseDialog, plot_spec, plot_hypnogram (style 'strip'), plot_raster.
%
%   HISTORY:
%       Created:  22 Jun 2026 (as ed_gui).
%       Renamed:  23 Jun 2026 -> gui_curate; ed/sSig dependency removed; Preset
%                               dropdown + curate_presets registry; events are
%                               generic (peak + optional start/stop + accepted).

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'basepath', @ischar);
addParameter(p, 'preset', '', @ischar);
addParameter(p, 'inputs', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'panels', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'events', [], @(x) isempty(x) || isnumeric(x) || isstruct(x));
addParameter(p, 'saveFcn', [], @(x) isempty(x) || isa(x, 'function_handle'));
addParameter(p, 'winPlot', 1.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'basename', '', @ischar);
addParameter(p, 'Visible', 'on', @(x) any(strcmpi(char(x), {'on', 'off'})));

parse(p, basepath, varargin{:});
basepath    = p.Results.basepath;
presetArg   = p.Results.preset;
inputsParam = p.Results.inputs;
panelsParam = p.Results.panels;
eventsParam = p.Results.events;
saveFcn     = p.Results.saveFcn;
winPlot     = p.Results.winPlot;
vis         = char(p.Results.Visible);

basename = p.Results.basename;
if isempty(basename)
    [~, basename] = fileparts(basepath);
end

%% ========================================================================
%  RESOLVE INITIAL CONFIG (preset or explicit inputs)
%  ========================================================================
presets = curate_presets();

if ~isempty(inputsParam) || ~isempty(eventsParam)
    presetName = 'Custom';
    inp = normalizeInputs(inputsParam);                 % [] -> empty struct array
    if ~isempty(eventsParam), inp = [inp, makeEventInputStruct(eventsParam, 'Events')]; end
    config = struct('inputs', {inp}, 'panels', {panelsParam}, 'saveFcn', saveFcn, 'winPlot', winPlot);
else
    presetName = resolvePreset(presetArg, presets, basepath, basename);
    config = loadPresetConfig(presets, presetName, basepath, basename, winPlot, []);
end

%% ========================================================================
%  STATE
%  ========================================================================
d = struct();
d.basepath   = basepath;
d.basename   = basename;
d.presets    = presets;
d.presetName = presetName;
d.clrAccept  = [0.10 0.55 0.10];
d.clrReject  = [0.65 0.15 0.15];
d.unit       = struct('wide', 'hr', 'narrow', 's');   % per-region x-axis units
d.hInd       = gobjects(0);
d.hCenter    = gobjects(0);
d.dirty      = false;
d.win        = winPlot;
d.winDefault = winPlot;
d.saveFcn    = [];
d = applyConfig(d, config);    % sets inputs, events, panels, saveFcn, t0, win

%% ========================================================================
%  FIGURE + LAYOUT
%  ========================================================================

% Build hidden, then reveal once everything is drawn (avoids incremental
% repaints while the control column and panels are populated).
hFig = uifigure('Name', sprintf('%s - Curate: %s', basename, presetName), ...
    'Position', [50, 50, 1320, 820], 'Visible', 'off', ...
    'WindowKeyPressFcn', @onKey, 'WindowScrollWheelFcn', @onScroll);

[~, gPlot, gCtrl, gActions] = tblgui.layout(hFig, 'CtrlWidth', 240, 'CtrlSide', 'left');

% --- Controls: preset ---
tblgui.labeledControl(gCtrl, 'label', 'PRESET', 'FontWeight', 'bold');
d.hPresetDD = tblgui.labeledControl(gCtrl, 'dropdown', '', ...
    'Items', presetItems(presets, presetName), 'Value', presetName, ...
    'ValueChangedFcn', @(s, ~) loadPreset(s.Value));

% --- Controls: load events / signals from the base workspace ---
% Events: an Nx1 [peak], Nx2 [start stop], or Nx3 [start peak stop] matrix.
% Signal: a numeric vector + the fs you type. Both are grabbed by variable name.
tblgui.labeledControl(gCtrl, 'label', 'LOAD (base workspace)', 'FontWeight', 'bold');
hL = tblgui.labeledControl(gCtrl, 'panel', '', 'RowHeight', 86);
gL = uigridlayout(hL, [3, 2], 'RowHeight', {'fit', 'fit', 'fit'}, ...
    'ColumnWidth', {'1x', '1x'}, 'Padding', 2, 'RowSpacing', 2, 'ColumnSpacing', 4);
d.hWsVar = uieditfield(gL, 'text', 'Placeholder', 'variable name');
d.hWsVar.Layout.Row = 1; d.hWsVar.Layout.Column = [1, 2];
lblFs = uilabel(gL, 'Text', 'fs (Hz)');
lblFs.Layout.Row = 2; lblFs.Layout.Column = 1;
d.hWsFs = uieditfield(gL, 'numeric', 'Value', 1250, 'Limits', [eps, Inf]);
d.hWsFs.Layout.Row = 2; d.hWsFs.Layout.Column = 2;
d.hWsEventsBtn = uibutton(gL, 'Text', 'Events', 'ButtonPushedFcn', @(~,~) onLoadWsEvents());
d.hWsEventsBtn.Layout.Row = 3; d.hWsEventsBtn.Layout.Column = 1;
d.hWsSignalBtn = uibutton(gL, 'Text', 'Signal', 'ButtonPushedFcn', @(~,~) onLoadWsSignal());
d.hWsSignalBtn.Layout.Row = 3; d.hWsSignalBtn.Layout.Column = 2;

% --- Controls: layout configuration (per region) ---
if strcmp(d.cfgRegion, 'wide'), nInitCfg = numel(d.wideP); else, nInitCfg = numel(d.narrowP); end
tblgui.labeledControl(gCtrl, 'label', 'LAYOUT', 'FontWeight', 'bold');

% Region dropdown and # Panels side by side.
hRP = tblgui.labeledControl(gCtrl, 'panel', '', 'RowHeight', 56);
gRP = uigridlayout(hRP, [2, 2], 'RowHeight', {'fit', 'fit'}, ...
    'ColumnWidth', {'1x', '1x'}, 'Padding', 2, 'RowSpacing', 2, 'ColumnSpacing', 4);
lblR = uilabel(gRP, 'Text', 'Region', 'FontWeight', 'bold');
lblR.Layout.Row = 1; lblR.Layout.Column = 1;
lblP = uilabel(gRP, 'Text', '# Panels', 'FontWeight', 'bold');
lblP.Layout.Row = 1; lblP.Layout.Column = 2;
d.hRegionDD = uidropdown(gRP, 'Items', {'Top', 'Bottom'}, 'Value', regDisp(d.cfgRegion), ...
    'ValueChangedFcn', @onRegionSel);
d.hRegionDD.Layout.Row = 2; d.hRegionDD.Layout.Column = 1;
d.hPanelsN = uieditfield(gRP, 'numeric', 'Limits', [1, 6], 'RoundFractionalValues', 'on', ...
    'Value', max(1, nInitCfg), 'ValueChangedFcn', @onPanelsN);
d.hPanelsN.Layout.Row = 2; d.hPanelsN.Layout.Column = 2;

hSrcPanel  = tblgui.labeledControl(gCtrl, 'panel', 'Panels', 'RowHeight', 160);
d.hSrcGrid = uigridlayout(hSrcPanel, [1, 1], 'Padding', 2, 'RowSpacing', 3, ...
    'ColumnWidth', {'1x'}, 'Scrollable', 'on');

% --- Controls: view (t0 / window in seconds; per-region x-unit) ---
tblgui.labeledControl(gCtrl, 'label', 'VIEW', 'FontWeight', 'bold');
hVF = tblgui.labeledControl(gCtrl, 'panel', '', 'RowHeight', 84);
gVF = uigridlayout(hVF, [3, 2], 'RowHeight', {'fit', 'fit', 'fit'}, ...
    'ColumnWidth', {'1x', '1x'}, 'Padding', 2, 'RowSpacing', 2, 'ColumnSpacing', 4);
lblT0 = uilabel(gVF, 'Text', 't0 (s)', 'FontWeight', 'bold');
lblT0.Layout.Row = 1; lblT0.Layout.Column = 1;
lblWin = uilabel(gVF, 'Text', 'Window (s)', 'FontWeight', 'bold');
lblWin.Layout.Row = 1; lblWin.Layout.Column = 2;
d.hT0 = uieditfield(gVF, 'numeric', 'Limits', [0, Inf], 'ValueDisplayFormat', '%.2f', ...
    'Value', d.t0, 'ValueChangedFcn', @onEditT0);
d.hT0.Layout.Row = 2; d.hT0.Layout.Column = 1;
d.hWin = uieditfield(gVF, 'numeric', 'Limits', [0.02, Inf], 'Value', d.win, ...
    'ValueChangedFcn', @onEditWin);
d.hWin.Layout.Row = 2; d.hWin.Layout.Column = 2;
lblU = uilabel(gVF, 'Text', 'X units', 'FontWeight', 'bold');
lblU.Layout.Row = 3; lblU.Layout.Column = 1;
d.hUnit = uidropdown(gVF, 'Items', {'ms', 's', 'min', 'hr'}, 'Value', d.unit.(d.cfgRegion), ...
    'ValueChangedFcn', @onUnit);
d.hUnit.Layout.Row = 3; d.hUnit.Layout.Column = 2;

hZ = tblgui.labeledControl(gCtrl, 'panel', '', 'RowHeight', 30);
gZ = uigridlayout(hZ, [1, 2], 'ColumnWidth', {'1x', '1x'}, 'Padding', 2, 'ColumnSpacing', 4);
bZi = uibutton(gZ, 'Text', 'Zoom In (+)',  'ButtonPushedFcn', @(~,~) zoomActive(1/1.5));
bZi.Layout.Column = 1;
bZo = uibutton(gZ, 'Text', 'Zoom Out (-)', 'ButtonPushedFcn', @(~,~) zoomActive(1.5));
bZo.Layout.Column = 2;
tblgui.labeledControl(gCtrl, 'button', '', 'Text', 'Reset (0)', 'ButtonPushedFcn', @(~,~) resetActive());
tblgui.labeledControl(gCtrl, 'spacer', '');

% --- Event module (always present; shows "/ 0" when a view has no events) ---
api = struct('prev', @() navStep(-1), 'next', @() navStep(1), ...
    'accept', @() setAccept(true), 'reject', @() setAccept(false), ...
    'save', @() onSave(), 'setIdx', @(v) setIdx(v));
d.ev = tblgui.eventPanel(gActions, api);

% --- Plot area: one tiledlayout holding every panel ---
d.hPanel = uipanel(gPlot, 'BorderType', 'none');
d = buildPlot(d);

hFig.UserData = d;

%% ========================================================================
%  INITIAL RENDER (Bottom first for faster reveal, then Top overview)
%  ========================================================================

populateSrc(hFig);
renderNarrow(hFig);
renderWideStatic(hFig);
finalizeInteractions(hFig);

drawnow;
di = hFig.UserData;
if numel(di.axWide)   > 1, linkaxes(di.axWide, 'x');   end
if numel(di.axNarrow) > 1, linkaxes(di.axNarrow, 'x'); end

updateMarker(hFig);
refreshEvent(hFig);

hFig.CloseRequestFcn = @(~,~) onClose();
hFig.Visible = vis;

%% ========================================================================
%  PRESET SWITCHING
%  ========================================================================

    function loadPreset(name)
        data = hFig.UserData;
        if strcmp(name, data.presetName), return; end
        if ~any(strcmp(name, {data.presets.name}))   % 'Custom'/'ws:...'/unknown
            data.hPresetDD.Value = data.presetName; return;
        end
        % offer to save unsaved curation before switching
        if data.dirty
            sel = tblgui.chooseDialog(hFig, 'Unsaved curation. Save before switching?', ...
                {'Save', 'Discard'});
            if isempty(sel), data.hPresetDD.Value = data.presetName; return; end
            if strcmp(sel, 'Save'), onSave(); data = hFig.UserData; end
        end
        dlg = busyOn(sprintf('Loading %s...', name));   % already-loaded signals are reused
        try
            cfg = loadPresetConfig(data.presets, name, data.basepath, data.basename, ...
                data.winDefault, data.pool);
        catch ME
            busyOff(dlg);
            tblgui.notify(hFig, sprintf('Preset "%s" failed: %s', name, ME.message), 'error');
            data.hPresetDD.Value = data.presetName; hFig.UserData = data; return;
        end
        data = applyConfig(data, cfg);
        data.presetName = name;
        data.dirty = false;
        hFig.Name = sprintf('%s - Curate: %s', data.basename, name);
        hFig.UserData = data;
        data.hRegionDD.Value = regDisp(data.cfgRegion);
        data.hPanelsN.Value  = max(1, numel(data.(regField(data.cfgRegion))));
        data.hUnit.Value     = data.unit.(data.cfgRegion);
        data.hWin.Value      = data.win;
        rebuildPlot();
        populateSrc(hFig);
        refreshEvent(hFig);
        busyOff(dlg);
    end

    % busy indicator around slow work (preset load / signal prep)
    function dlg = busyOn(msg)
        dlg = [];
        try
            dlg = uiprogressdlg(hFig, 'Indeterminate', 'on', 'Message', msg, 'Title', 'Working');
            drawnow;
        catch
        end
    end
    function busyOff(dlg)
        try, if ~isempty(dlg) && isvalid(dlg), close(dlg); end, catch, end
    end

    % set the active event set (from the workspace) while keeping the signal pool
    function setEvents(evInput, sv, label)
        data = hFig.UserData;
        ev = normalizeInputs(evInput);
        data.inputs   = [data.pool, ev(1)];
        data.saveFcn  = sv;
        data.presetName = label;
        data.dirty    = false;
        data.Tend_s   = max(eps, computeTend(data.inputs));
        [data.ed, data.accepted, data.nEvents, data.hasEvents] = extractEvents(data.inputs);
        data.currIdx  = 1;
        if data.hasEvents && data.nEvents > 0, data.t0 = data.ed.peakTime(1); else, data.t0 = data.Tend_s / 2; end
        if ~any(strcmp({data.wideP.source}, 'eventTicks')) && ~any(strcmp({data.narrowP.source}, 'eventTicks'))
            data.wideP = [data.wideP, makePanel('eventTicks', 'wide', data.inputs)];
        end
        items = data.hPresetDD.Items;
        if ~any(strcmp(label, items)), data.hPresetDD.Items = [items, {label}]; end
        data.hPresetDD.Value = label;
        hFig.Name = sprintf('%s - Curate: %s', data.basename, label);
        hFig.UserData = data;
        rebuildPlot();
        populateSrc(hFig);
        refreshEvent(hFig);
    end

    function onLoadWsEvents()
        data = hFig.UserData;
        name = strtrim(data.hWsVar.Value);
        if isempty(name), tblgui.notify(hFig, 'Enter a base workspace variable name.', 'warning'); return; end
        try, M = evalin('base', name); catch, tblgui.notify(hFig, sprintf('No base variable "%s".', name), 'error'); return; end
        try, evd = eventsFromMatrix(M); catch ME, tblgui.notify(hFig, ME.message, 'error'); return; end
        evInput = struct('name', 'eventTicks', 'type', 'eventTicks', 'data', evd, 'label', name);
        sv = @(acc) assignin('base', [name, '_accepted'], logical(acc(:)));
        setEvents(evInput, sv, ['ws:' name]);
        tblgui.notify(hFig, sprintf('Loaded %d events from "%s"; Save writes %s_accepted to base.', ...
            numel(evd.peakTime), name, name), 'success');
    end

    function onLoadWsSignal()
        data = hFig.UserData;
        name = strtrim(data.hWsVar.Value);
        if isempty(name), tblgui.notify(hFig, 'Enter a base workspace variable name.', 'warning'); return; end
        try, v = evalin('base', name); catch, tblgui.notify(hFig, sprintf('No base variable "%s".', name), 'error'); return; end
        if ~isnumeric(v) || ~isvector(v) || numel(v) < 2
            tblgui.notify(hFig, 'Signal must be a numeric vector.', 'error'); return;
        end
        fsv = data.hWsFs.Value;
        tr = normalizeInputs(struct('name', name, 'type', 'trace', 'data', double(v(:)), ...
            'fs', fsv, 'ylim', prctile(double(v(:)), [0.1, 99.9]), 'label', name));
        data.pool = mergePool(data.pool, tr);
        evIx = find(strcmp({data.inputs.type}, 'eventTicks'), 1);
        if ~isempty(evIx), data.inputs = [data.pool, data.inputs(evIx)]; else, data.inputs = data.pool; end
        data.Tend_s = max(eps, computeTend(data.inputs));
        hFig.UserData = data;
        populateSrc(hFig);
        tblgui.notify(hFig, sprintf('Loaded "%s" (%d samp @ %g Hz). Pick it in a panel dropdown.', ...
            name, numel(v), fsv), 'success');
    end

%% ========================================================================
%  PLOT CONSTRUCTION (single tiledlayout)
%  ========================================================================

    function dd = buildPlot(dd)
        % (re)build the tiledlayout from dd.wideP / dd.narrowP. Per-panel
        % heights become integer row spans; a 1-row black tile divides the
        % regions; a blank gap above it holds the Top x-axis. A larger K makes
        % the divider proportionally thinner while preserving panel ratios.
        if isfield(dd, 'tl') && ~isempty(dd.tl) && isvalid(dd.tl), delete(dd.tl); end
        K = 20;
        nW = numel(dd.wideP); nN = numel(dd.narrowP);
        wSpans = ones(1, nW); for i = 1:nW, wSpans(i) = max(1, round(dd.wideP(i).height * K)); end
        nSpans = ones(1, nN); for i = 1:nN, nSpans(i) = max(1, round(dd.narrowP(i).height * K)); end
        hasDiv = nW > 0 && nN > 0;
        base = sum(wSpans) + sum(nSpans) + hasDiv;
        xgap = 0;
        if hasDiv, xgap = max(4, round(0.06 * base)); end
        totalRows = base + xgap;

        tl = tiledlayout(dd.hPanel, max(1, totalRows), 1, ...
            'TileSpacing', 'none', 'Padding', 'tight');
        dd.tl = tl;

        r = 1;
        dd.axWide = gobjects(1, nW);
        for i = 1:nW
            ax = nexttile(tl, r, [wSpans(i), 1]); r = r + wSpans(i);
            dd.wideP(i).ax = ax; dd.axWide(i) = ax;
        end

        r = r + xgap;          % blank gap holds the wide x-tick labels

        dd.divider = gobjects(0);
        if hasDiv
            dax = nexttile(tl, r, [1, 1]); r = r + 1;
            set(dax, 'Color', [0 0 0], 'XTick', [], 'YTick', [], ...
                'XColor', 'none', 'YColor', 'none', 'Box', 'off');
            try, disableDefaultInteractivity(dax); catch, end
            dax.Toolbar = []; dax.HitTest = 'off';
            dd.divider = dax;
        end

        dd.axNarrow = gobjects(1, nN);
        for i = 1:nN
            ax = nexttile(tl, r, [nSpans(i), 1]); r = r + nSpans(i);
            dd.narrowP(i).ax = ax; dd.axNarrow(i) = ax;
        end
    end

%% ========================================================================
%  RENDER LAYER
%  ========================================================================

    function renderWideStatic(fig)
        % Top overview: drawn once (full session). x in the Top's unit. Each
        % panel gets a cursor line (t0) + a band (the Bottom window).
        data = fig.UserData;
        xf = unitSec(data.unit.wide);
        a = 0; b = data.Tend_s;
        data.hInd    = gobjects(1, numel(data.wideP));
        data.hCenter = gobjects(1, numel(data.wideP));
        for i = 1:numel(data.wideP)
            pn = data.wideP(i); ax = pn.ax;
            cla(ax); hold(ax, 'on');
            drawPanel(data, ax, pn.source, a, b, xf);
            title(ax, ''); ylabel(ax, pn.label);
            set(get(ax, 'YLabel'), 'Color', [0.15 0.15 0.15]);
            ax.XLim = [a, b] / xf;
            data.hInd(i) = xregion(ax, a/xf, a/xf + eps, 'FaceColor', [0 0.2 0.8], ...
                'FaceAlpha', 0.15, 'HandleVisibility', 'off');
            data.hCenter(i) = xline(ax, 0, 'Color', [0 0.25 0.9], ...
                'LineWidth', 1.25, 'HandleVisibility', 'off');
            mute(ax);
        end
        fig.UserData = data;
        hideOuterXTicks(data.axWide, ['Time (' data.unit.wide ')'], 'bottom');
        updateMarker(fig);
    end

    function renderTicksPanel(fig)
        % recolour the Top eventTicks panel(s) after an acceptance change
        data = fig.UserData;
        xf = unitSec(data.unit.wide);
        for i = 1:numel(data.wideP)
            if strcmp(data.wideP(i).source, 'eventTicks')
                ax = data.wideP(i).ax;
                cla(ax); hold(ax, 'on');
                drawTicks(ax, data, xf);
                ylabel(ax, data.wideP(i).label);
                set(get(ax, 'YLabel'), 'Color', [0.15 0.15 0.15]);
                ax.XLim = [0, data.Tend_s] / xf;
                data.hInd(i) = xregion(ax, 0, eps, 'FaceColor', [0 0.2 0.8], ...
                    'FaceAlpha', 0.15, 'HandleVisibility', 'off');
                data.hCenter(i) = xline(ax, 0, 'Color', [0 0.25 0.9], ...
                    'LineWidth', 1.25, 'HandleVisibility', 'off');
                mute(ax);
            end
        end
        fig.UserData = data;
        hideOuterXTicks(data.axWide, ['Time (' data.unit.wide ')'], 'bottom');
        updateMarker(fig);
    end

    function renderNarrow(fig)
        % Bottom window [t0 +/- win/2], x in the Bottom's unit; redrawn as t0
        % or win changes (the Top is untouched, only its markers move).
        data = fig.UserData;
        if isempty(data.narrowP), return; end
        xf = unitSec(data.unit.narrow);
        ws = data.t0 - data.win / 2;
        we = data.t0 + data.win / 2;
        if ws < 0,            we = we - ws;              ws = 0; end
        if we > data.Tend_s,  ws = ws - (we - data.Tend_s); we = data.Tend_s; end
        ws = max(0, ws);
        ev = data.currIdx;
        showEv = data.hasEvents && data.nEvents > 0 && ...
            data.ed.peakTime(ev) >= ws && data.ed.peakTime(ev) <= we;
        for i = 1:numel(data.narrowP)
            pn = data.narrowP(i); ax = pn.ax;
            cla(ax); hold(ax, 'on');
            drawPanel(data, ax, pn.source, ws, we, xf);
            if showEv
                xline(ax, eventMarks(data.ed, ev) / xf, '--b', 'HandleVisibility', 'off');
            end
            ax.XLim = [ws, we] / xf;
            ylabel(ax, pn.label);
            set(get(ax, 'YLabel'), 'Color', [0.15 0.15 0.15]);
            mute(ax);
        end
        hideOuterXTicks(data.axNarrow, ['Time (' data.unit.narrow ')'], 'bottom');
        fig.UserData = data;
    end

    function drawPanel(data, ax, source, a, b, xf)
        % dispatch one panel by the source input's TYPE; x in display units
        % (= seconds / xf), covering the seconds range [a, b]
        inp = getInput(data.inputs, source);
        if isempty(inp), return; end
        switch inp.type
            case 'trace',      drawTrace(ax, inp, a, b, xf);
            case 'spec',       drawSpec(ax, inp, xf);
            case 'hypnogram',  drawHypno(ax, inp, xf);
            case 'eventTicks', if data.hasEvents, drawTicks(ax, data, xf); end
            case 'raster',     drawRaster(ax, inp, a, b, xf);
        end
    end

    function updateMarker(fig)
        % move the Top cursor line (t0) + window band (the Bottom window)
        data = fig.UserData;
        xf = unitSec(data.unit.wide);
        c  = data.t0 / xf;
        lo = max(0, data.t0 - data.win / 2) / xf;
        hi = max(lo + eps, (data.t0 + data.win / 2) / xf);
        for k = 1:numel(data.hInd)
            if isvalid(data.hInd(k)),    data.hInd(k).Value = [lo, hi]; end
            if isvalid(data.hCenter(k)), data.hCenter(k).Value = c;     end
        end
        if isfield(data, 'hT0') && isvalid(data.hT0), data.hT0.Value = data.t0; end
    end

    function refreshEvent(fig)
        % event module shows index / total only; accept-reject is read off the
        % plot panels, so no status text is set here
        data = fig.UserData;
        if isempty(data.ev), return; end
        if ~data.hasEvents || data.nEvents == 0
            data.ev.refresh(0, 0, false, '');
            return;
        end
        data.ev.refresh(data.currIdx, data.nEvents, data.accepted(data.currIdx), '');
    end

%% ========================================================================
%  NAVIGATION / STATE (t0 drives the Bottom; the Top is a fixed overview)
%  ========================================================================

    function navStep(step)
        data = hFig.UserData;
        if ~data.hasEvents || data.nEvents == 0, return; end
        data.currIdx = min(max(1, data.currIdx + step), data.nEvents);
        data.t0 = data.ed.peakTime(data.currIdx);
        hFig.UserData = data;
        renderNarrow(hFig); updateMarker(hFig); refreshEvent(hFig);
    end

    function setIdx(v)
        data = hFig.UserData;
        if ~data.hasEvents || data.nEvents == 0, return; end
        data.currIdx = min(max(1, round(v)), data.nEvents);
        data.t0 = data.ed.peakTime(data.currIdx);
        hFig.UserData = data;
        renderNarrow(hFig); updateMarker(hFig); refreshEvent(hFig);
    end

    function jumpToTime(tSec)
        % free move: set t0, re-centre the Bottom, leave currIdx untouched
        data = hFig.UserData;
        data.t0 = min(max(0, tSec), data.Tend_s);
        hFig.UserData = data;
        renderNarrow(hFig); updateMarker(hFig);
    end

    function setAccept(tf)
        data = hFig.UserData;
        if ~data.hasEvents || data.nEvents == 0, return; end
        data.accepted(data.currIdx) = tf;
        data.dirty = true;
        hFig.UserData = data;
        renderTicksPanel(hFig);           % recolour Top eventTicks
        if data.currIdx < data.nEvents
            navStep(1);                   % auto-advance (re-renders Bottom + markers)
        else
            renderNarrow(hFig); refreshEvent(hFig);
        end
    end

    function onEditT0(src, ~)
        jumpToTime(src.Value);
    end

    function onEditWin(src, ~)
        data = hFig.UserData;
        data.win = max(0.02, src.Value);
        hFig.UserData = data;
        renderNarrow(hFig); updateMarker(hFig);
    end

    function onUnit(src, ~)
        data = hFig.UserData;
        data.unit.(data.cfgRegion) = src.Value;
        hFig.UserData = data;
        if strcmp(data.cfgRegion, 'wide')
            renderWideStatic(hFig);       % redraws + updateMarker
        else
            renderNarrow(hFig);
        end
    end

%% ========================================================================
%  ZOOM (active region: Top scales its x-range, Bottom scales the window)
%  ========================================================================

    function zoomActive(f)
        data = hFig.UserData;
        if strcmp(data.cfgRegion, 'wide'), zoomTop(f, []); else, zoomWin(f); end
    end

    function resetActive()
        data = hFig.UserData;
        if strcmp(data.cfgRegion, 'wide'), resetOverview(); else, resetWin(); end
    end

    function zoomTop(f, centerDisp)
        % scale the Top x-range by f (f<1 zooms in) about centerDisp (display
        % units); default centre is t0
        data = hFig.UserData;
        if isempty(data.axWide), return; end
        ax = data.axWide(1);
        xl = xlim(ax);
        xf = unitSec(data.unit.wide);
        if isempty(centerDisp), c = data.t0 / xf; else, c = centerDisp; end
        if ~isfinite(c) || c < xl(1) || c > xl(2), c = mean(xl); end
        lo = max(0, c - (c - xl(1)) * f);
        hi = min(data.Tend_s / xf, c + (xl(2) - c) * f);
        if hi > lo, xlim(ax, [lo, hi]); end
    end

    function zoomWin(f)
        data = hFig.UserData;
        data.win = min(max(0.02, data.win * f), data.Tend_s);
        hFig.UserData = data;
        data.hWin.Value = data.win;
        renderNarrow(hFig); updateMarker(hFig);
    end

    function resetOverview()
        data = hFig.UserData;
        if isempty(data.axWide), return; end
        xlim(data.axWide(1), [0, data.Tend_s / unitSec(data.unit.wide)]);
    end

    function resetWin()
        data = hFig.UserData;
        data.win = data.winDefault;
        hFig.UserData = data;
        data.hWin.Value = data.win;
        renderNarrow(hFig); updateMarker(hFig);
    end

%% ========================================================================
%  INPUT (keyboard / scroll / click)
%  ========================================================================

    function onKey(~, evt)
        if any(strcmpi(evt.Modifier, 'control')) && strcmpi(evt.Key, 's')
            onSave(); return;
        end
        switch evt.Key
            case 'leftarrow',  navStep(-1);
            case 'rightarrow', navStep(1);
            case 'uparrow',    setAccept(true);
            case 'downarrow',  setAccept(false);
            case {'equal', 'add'},       zoomActive(1 / 1.5);
            case {'hyphen', 'subtract'}, zoomActive(1.5);
            case {'0', 'numpad0'},       resetActive();
        end
    end

    function onScroll(~, evt)
        data = hFig.UserData;
        cp = hFig.CurrentPoint;
        if evt.VerticalScrollCount > 0, f = 1.4; else, f = 1 / 1.4; end
        if ~isempty(data.axNarrow) && pointerOver(data.axNarrow, cp)
            zoomWin(f);
        elseif ~isempty(data.axWide) && pointerOver(data.axWide, cp)
            zoomTop(f, data.axWide(1).CurrentPoint(1, 1));
        end
    end

    function tf = pointerOver(axs, cp)
        tf = false;
        for i = 1:numel(axs)
            if ~isvalid(axs(i)), continue; end
            r = getpixelposition(axs(i), true);
            if cp(1) >= r(1) && cp(1) <= r(1) + r(3) && ...
                    cp(2) >= r(2) && cp(2) <= r(2) + r(4)
                tf = true; return;
            end
        end
    end

    function onClickAxis(region, ax)
        % click any panel -> move t0 (in that panel's region units) -> Bottom
        data = hFig.UserData;
        jumpToTime(ax.CurrentPoint(1, 1) * unitSec(data.unit.(region)));
    end

    function finalizeInteractions(fig)
        data = fig.UserData;
        for k = 1:numel(data.axWide),   armAxis(data.axWide(k),   'wide');   end
        for k = 1:numel(data.axNarrow), armAxis(data.axNarrow(k), 'narrow'); end
        data.jumpToTimeFcn = @jumpToTime;       % exposed for hosts / tests
        fig.UserData = data;
    end

    function armAxis(ax, region)
        if ~isvalid(ax), return; end
        try, disableDefaultInteractivity(ax); catch, end
        ax.Toolbar = [];
        ax.PickableParts = 'all';
        ax.HitTest = 'on';
        ax.ButtonDownFcn = @(s, ~) onClickAxis(region, s);
    end

%% ========================================================================
%  LAYOUT RECONFIGURATION (region-switch)
%  ========================================================================

    function onRegionSel(src, ~)
        data = hFig.UserData;
        data.cfgRegion = regFromDisp(src.Value);
        hFig.UserData = data;
        data.hPanelsN.Value = max(1, numel(data.(regField(data.cfgRegion))));
        populateSrc(hFig);
        data.hUnit.Value = data.unit.(data.cfgRegion);
    end

    function onPanelsN(src, ~)
        data = hFig.UserData;
        region = data.cfgRegion;
        f = regField(region);
        pArr = data.(f);
        n = round(src.Value);
        cur = numel(pArr);
        if n > cur
            opts = availSources(data.inputs);
            for i = cur + 1:n
                pArr(i) = makePanel(opts{1}, region, data.inputs);
            end
        elseif n < cur
            pArr = pArr(1:max(1, n));
        end
        data.(f) = pArr;
        hFig.UserData = data;
        rebuildPlot();
        populateSrc(hFig);
    end

    function onSrcChange(k, src)
        data = hFig.UserData;
        region = data.cfgRegion;
        f = regField(region);
        pArr = data.(f);
        if k > numel(pArr), return; end
        inp = getInput(data.inputs, src.Value);
        pArr(k).source = src.Value;
        pArr(k).label  = inp.label;
        pArr(k).height = inp.height;
        data.(f) = pArr;
        hFig.UserData = data;
        rebuildPlot();
    end

    function populateSrc(fig)
        data = fig.UserData;
        region = data.cfgRegion;
        pArr = data.(regField(region));
        g = data.hSrcGrid;
        delete(g.Children);
        nP = max(1, numel(pArr));
        g.RowHeight = repmat({'fit'}, 1, nP);
        opts = availSources(data.inputs);
        for k = 1:numel(pArr)
            val = pArr(k).source;
            if ~any(strcmp(val, opts)), val = opts{1}; end
            dd = uidropdown(g, 'Items', opts, 'Value', val, ...
                'ValueChangedFcn', @(s, ~) onSrcChange(k, s));
            dd.Layout.Row = k; dd.Layout.Column = 1;
        end
    end

    function rebuildPlot()
        dd = hFig.UserData;
        dd = buildPlot(dd);
        hFig.UserData = dd;
        renderNarrow(hFig);
        renderWideStatic(hFig);
        finalizeInteractions(hFig);
        drawnow;
        dd = hFig.UserData;
        if numel(dd.axWide) > 1, linkaxes(dd.axWide, 'x'); end
        if numel(dd.axNarrow) > 1, linkaxes(dd.axNarrow, 'x'); end
        updateMarker(hFig);
        refreshEvent(hFig);
    end

%% ========================================================================
%  SAVE / CLOSE
%  ========================================================================

    function onSave()
        data = hFig.UserData;
        if ~data.hasEvents, return; end
        if isempty(data.saveFcn)
            tblgui.notify(hFig, 'No save target for this view.', 'warning');
            return;
        end
        data.saveFcn(data.accepted(:));
        data.dirty = false;
        hFig.UserData = data;
        tblgui.notify(hFig, sprintf('Saved %d accepted / %d events.', ...
            sum(data.accepted), data.nEvents), 'success');
    end

    function onClose()
        data = hFig.UserData;
        if isfield(data, 'dirty') && data.dirty
            sel = tblgui.chooseDialog(hFig, 'Unsaved curation. Save before closing?', ...
                {'Save and close', 'Close without saving'});
            if isempty(sel), return; end
            if strcmp(sel, 'Save and close'), onSave(); end
        end
        delete(hFig);
    end

end     % MAIN

%% ========================================================================
%  CONFIG / PRESET HELPERS (pure)
%  ========================================================================

function name = resolvePreset(arg, presets, basepath, basename)
% chosen preset, else the first preset whose file exists, else the first
names = {presets.name};
if ~isempty(arg)
    ix = find(strcmpi(arg, names), 1);
    if ~isempty(ix), name = names{ix}; return; end
end
for i = 1:numel(presets)
    if isfile(fullfile(basepath, [basename, '.', presets(i).var, '.mat']))
        name = presets(i).name; return;
    end
end
name = presets(1).name;
end

function config = loadPresetConfig(presets, name, basepath, basename, winDefault, pool)
ix = find(strcmp({presets.name}, name), 1);
if isempty(ix), error('gui_curate:preset', 'unknown preset "%s"', name); end
config = presets(ix).load(basepath, basename, pool);
if ~isfield(config, 'inputs') || isempty(config.inputs)
    error('gui_curate:preset', 'preset "%s" returned no inputs', name);
end
if ~isfield(config, 'panels'),  config.panels  = []; end
if ~isfield(config, 'saveFcn'), config.saveFcn = []; end
if ~isfield(config, 'winPlot') || isempty(config.winPlot), config.winPlot = winDefault; end
end

function items = presetItems(presets, presetName)
items = {presets.name};
if strcmp(presetName, 'Custom'), items = [{'Custom'}, items]; end
end

function d = applyConfig(d, config)
% apply a resolved config onto state struct d. SIGNAL inputs are merged into a
% persistent pool (loaded once, reused across presets); the eventTicks input is
% the active event set. d.inputs = pool + active events, so every loaded signal
% stays available regardless of which preset is active.
allIn = normalizeInputs(config.inputs);
isEvt = strcmp({allIn.type}, 'eventTicks');
sigIn = allIn(~isEvt);
evtIn = allIn(isEvt);
if ~isfield(d, 'pool') || isempty(d.pool)
    d.pool = sigIn;
else
    d.pool = mergePool(d.pool, sigIn);
end
if ~isempty(evtIn), d.inputs = [d.pool, evtIn(1)]; else, d.inputs = d.pool; end

d.saveFcn = config.saveFcn;
if isfield(config, 'winPlot') && ~isempty(config.winPlot)
    d.win = config.winPlot; d.winDefault = config.winPlot;
end
d.Tend_s = max(eps, computeTend(d.inputs));

[d.ed, d.accepted, d.nEvents, d.hasEvents] = extractEvents(d.inputs);
d.currIdx = 1;
if d.hasEvents && d.nEvents > 0, d.t0 = d.ed.peakTime(1); else, d.t0 = d.Tend_s / 2; end

panels = normalizePanels(config, d.inputs);
d.wideP   = panels(strcmp({panels.region}, 'wide'));
d.narrowP = panels(strcmp({panels.region}, 'narrow'));
if isempty(d.wideP) && ~isempty(d.narrowP), d.cfgRegion = 'narrow'; else, d.cfgRegion = 'wide'; end
end

function [ed, accepted, nEv, has] = extractEvents(inputs)
% pull the eventTicks input's data into the live event state
ed = []; accepted = logical([]); nEv = 0; has = false;
ix = find(strcmp({inputs.type}, 'eventTicks'), 1);
if isempty(ix), return; end
data = inputs(ix).data;
if ~isstruct(data) || ~isfield(data, 'peakTime') || isempty(data.peakTime), return; end
ed = data; ed.peakTime = data.peakTime(:);
nEv = numel(ed.peakTime); has = true;
if isfield(data, 'accepted') && numel(data.accepted) == nEv
    accepted = logical(data.accepted(:));
else
    accepted = true(nEv, 1);
end
end

function x = eventMarks(ed, ev)
% marks to draw for event ev: [start peak stop] if start/stop exist, else peak
x = ed.peakTime(ev);
if isfield(ed, 'times') && ~isempty(ed.times) && size(ed.times, 1) >= ev
    x = [ed.times(ev, 1), ed.peakTime(ev), ed.times(ev, 2)];
end
end

function panels = normalizePanels(config, inputs)
% build full panel structs from config.panels (.source/.region) or, if none,
% from the inputs' default regions; then drop panels with no matching input
if isempty(config.panels)
    panels = defaultPanels(inputs);
else
    cp = config.panels;
    P = cell(1, numel(cp));
    for i = 1:numel(cp)
        reg = 'narrow';
        if isfield(cp, 'region') && ~isempty(cp(i).region), reg = cp(i).region; end
        P{i} = makePanel(cp(i).source, reg, inputs);
    end
    panels = [P{:}];
end
panels = filterAvail(panels, inputs);
panels = addAxField(panels);
end

function ev = eventsFromMatrix(M)
% Nx1 [peak], Nx2 [start stop], or Nx3 [start peak stop] -> events struct
if ~isnumeric(M) || isempty(M)
    error('gui_curate:events', 'events must be a non-empty numeric matrix');
end
M = double(M);
if isvector(M)
    ev = struct('peakTime', M(:));
elseif size(M, 2) == 3
    ev = struct('peakTime', M(:, 2), 'times', [M(:, 1), M(:, 3)]);
elseif size(M, 2) == 2
    ev = struct('peakTime', mean(M, 2), 'times', M);
else
    error('gui_curate:events', 'events matrix must be Nx1, Nx2, or Nx3');
end
end

function ev = makeEventInputStruct(eventsData, label)
% a normalized eventTicks input from a matrix (Nx1/2/3) or an events struct
if isnumeric(eventsData), evd = eventsFromMatrix(eventsData); else, evd = eventsData; end
ev = normalizeInputs(struct('name', 'eventTicks', 'type', 'eventTicks', 'data', evd, 'label', label));
end

function pool = mergePool(pool, sigIn)
% append signal inputs not already present (cache; dedup by name)
for i = 1:numel(sigIn)
    if isempty(pool) || ~any(strcmp(sigIn(i).name, {pool.name}))
        pool(end + 1) = sigIn(i); %#ok<AGROW>
    end
end
end

%% ========================================================================
%  PANEL DRAW FUNCTIONS (selected by input type)
%  ========================================================================

function drawTrace(ax, inp, a, b, xf)
% windowed raw trace, decimated for display, x in display units
sig = inp.data; fs = inp.fs;
s1 = max(1, floor(a * fs) + 1);
s2 = min(numel(sig), ceil(b * fs) + 1);
if s2 < s1, s2 = s1; end
rng = s1:s2;
t = ((rng - 1) / fs) / xf;
np = numel(rng); maxPts = 20000;
if np > maxPts
    st = ceil(np / maxPts);
    plot(ax, t(1:st:end), sig(rng(1:st:end)), 'Color', inp.clr);
else
    plot(ax, t, sig(rng), 'Color', inp.clr);
end
if ~isempty(inp.ylim), ax.YLim = inp.ylim; end
end

function drawSpec(ax, inp, xf)
plot_spec(inp.data, 'axh', ax, 'saveFig', false, 'xtime', xf);
end

function drawHypno(ax, inp, xf)
% bout times arrive in hours; convert to display units. Pin sstates to the
% number of bout-cells provided so the strip is independent of cfg.nstates.
bt = cellfun(@(x) x * 3600 / xf, inp.data, 'uni', false);
plot_hypnogram('boutTimes', bt, 'sstates', 1:numel(bt), 'style', 'strip', 'hAx', ax);
end

function drawRaster(ax, inp, a, b, xf)
% inward ticks so the unit numbers stay but no tick marks protrude left
spk = cellfun(@(s) s(s >= a & s <= b) / xf, inp.data, 'uni', false);
if all(cellfun(@isempty, spk)), return; end   % skip (avoids plot_raster's empty warning)
plot_raster(spk, 'hAx', ax, 'xLim', [a, b] / xf, 'flgLbls', false, 'tickDir', 'in');
end

function drawTicks(ax, data, xf)
% draw only the ACCEPTED events; rejected are intentionally not shown (load a
% different event list from the workspace to inspect a different set)
drawTickLine(ax, data.ed.peakTime(data.accepted) / xf, [0.15 0.15 0.15]);
ylim(ax, [0, 1]);
end

function drawTickLine(ax, x, clr)
if isempty(x), return; end
x = x(:)';
X = [x; x; nan(1, numel(x))];
Y = repmat([0; 1; NaN], 1, numel(x));
line(ax, X(:), Y(:), 'Color', clr, 'HandleVisibility', 'off');
end

%% ========================================================================
%  INPUTS MODEL (pure)
%  ========================================================================

function b = blankInput()
b = struct('name', '', 'type', '', 'data', [], 'fs', NaN, 'ylim', [], ...
    'clr', 'k', 'label', '', 'height', 1, 'defRegion', 'narrow', 'defOrder', 99);
end

function out = normalizeInputs(in)
% fill defaults on an inputs struct array (requires .name .type)
n = numel(in);
out = repmat(blankInput(), 1, n);
for i = 1:n
    s = in(i);
    if ~isfield(s, 'name') || isempty(s.name), error('gui_curate:inputs', 'input %d needs a name', i); end
    if ~isfield(s, 'type') || isempty(s.type), error('gui_curate:inputs', 'input %d needs a type', i); end
    td = typeDefaults(s.type);
    out(i).name      = s.name;
    out(i).type      = s.type;
    out(i).data      = pick(s, 'data', []);
    out(i).fs        = pick(s, 'fs', NaN);
    out(i).ylim      = pick(s, 'ylim', []);
    out(i).clr       = pick(s, 'clr', 'k');
    out(i).label     = pick(s, 'label', td.label);
    if isempty(out(i).label), out(i).label = s.name; end
    out(i).height    = pick(s, 'height', td.height);
    out(i).defRegion = pick(s, 'defRegion', td.defRegion);
    out(i).defOrder  = pick(s, 'defOrder', td.defOrder);
end
end

function v = pick(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end

function td = typeDefaults(type)
% per-type fallback height / label / default region / stacking order
switch type
    case 'spec',       td = struct('height', 1.4,  'label', 'Freq (Hz)', 'defRegion', 'wide',   'defOrder', 20);
    case 'hypnogram',  td = struct('height', 0.28, 'label', 'State',     'defRegion', 'wide',   'defOrder', 10);
    case 'eventTicks', td = struct('height', 0.28, 'label', 'Events',    'defRegion', 'wide',   'defOrder', 40);
    case 'raster',     td = struct('height', 1.2,  'label', 'Units',     'defRegion', 'narrow', 'defOrder', 60);
    case 'trace',      td = struct('height', 1.0,  'label', '',          'defRegion', 'narrow', 'defOrder', 50);
    otherwise,         td = struct('height', 1.0,  'label', '',          'defRegion', 'narrow', 'defOrder', 99);
end
end

function T = computeTend(inputs)
% session length [s] = the largest time spanned by any input
T = eps;
for i = 1:numel(inputs)
    inp = inputs(i);
    switch inp.type
        case 'trace'
            if ~isempty(inp.data) && isfinite(inp.fs) && inp.fs > 0
                T = max(T, (numel(inp.data) - 1) / inp.fs);
            end
        case 'spec'
            if isstruct(inp.data) && isfield(inp.data, 'tstamps') && ~isempty(inp.data.tstamps)
                T = max(T, max(inp.data.tstamps(:)));
            end
        case 'hypnogram'
            for s = 1:numel(inp.data)
                if ~isempty(inp.data{s}), T = max(T, max(inp.data{s}(:)) * 3600); end
            end
        case 'raster'
            for u = 1:numel(inp.data)
                if ~isempty(inp.data{u}), T = max(T, max(inp.data{u})); end
            end
        case 'eventTicks'
            if isstruct(inp.data)
                if isfield(inp.data, 'times') && ~isempty(inp.data.times)
                    T = max(T, max(inp.data.times(:)));
                elseif isfield(inp.data, 'peakTime') && ~isempty(inp.data.peakTime)
                    T = max(T, max(inp.data.peakTime(:)));
                end
            end
    end
end
end

function inp = getInput(inputs, name)
ix = find(strcmp({inputs.name}, name), 1);
if isempty(ix), inp = []; else, inp = inputs(ix); end
end

%% ========================================================================
%  PANEL / LAYOUT HELPERS (pure)
%  ========================================================================

function panels = defaultPanels(inputs)
% fallback layout (custom path with no panels): each input in its default
% region, ordered by defOrder; a hypnogram is also added to the narrow region
P = {};
for region = {'wide', 'narrow'}
    reg = region{1};
    sel = find(strcmp({inputs.defRegion}, reg));
    if isempty(sel), continue; end
    [~, ord] = sort([inputs(sel).defOrder]);
    sel = sel(ord);
    for k = 1:numel(sel)
        P{end + 1} = makePanel(inputs(sel(k)).name, reg, inputs); %#ok<AGROW>
    end
end
h = find(strcmp({inputs.type}, 'hypnogram'), 1);
if ~isempty(h) && ~strcmp(inputs(h).defRegion, 'narrow')
    P{end + 1} = makePanel(inputs(h).name, 'narrow', inputs);
end
panels = [P{:}];
end

function panels = filterAvail(panels, inputs)
% keep only panels whose source names an existing input
if isempty(panels), return; end
keep = ismember({panels.source}, {inputs.name});
panels = panels(keep);
end

function panels = addAxField(panels)
% guarantee an .ax field on every entry (and on an empty array)
if isempty(panels)
    panels = struct('source', {}, 'region', {}, 'height', {}, 'label', {}, 'ax', {});
    return;
end
for i = 1:numel(panels)
    panels(i).ax = gobjects(1);
end
end

function pn = makePanel(name, region, inputs)
inp = getInput(inputs, name);
if isempty(inp)
    pn = struct('source', name, 'region', region, 'height', 1, 'label', name, 'ax', gobjects(1));
else
    pn = struct('source', name, 'region', region, 'height', inp.height, ...
        'label', inp.label, 'ax', gobjects(1));
end
end

function s = availSources(inputs)
% every input is offered as a source in both regions (stable defOrder order)
[~, ord] = sort([inputs.defOrder]);
s = {inputs(ord).name};
end

function s = regDisp(region)
if strcmp(region, 'wide'), s = 'Top'; else, s = 'Bottom'; end
end

function r = regFromDisp(disp)
if strcmpi(disp, 'Top'), r = 'wide'; else, r = 'narrow'; end
end

function f = regField(region)
if strcmp(region, 'wide'), f = 'wideP'; else, f = 'narrowP'; end
end

function xf = unitSec(unit)
% seconds per display unit
switch unit
    case 'ms',  xf = 1e-3;
    case 's',   xf = 1;
    case 'min', xf = 60;
    case 'hr',  xf = 3600;
    otherwise,  xf = 1;
end
end

function hideOuterXTicks(axs, label, location)
% show x-tick labels + the axis label only on the region's outer axis. For the
% Top (wide) region this is called with 'bottom' so the overview's x-axis prints
% on the bottom of its lowest panel, just above the divider.
n = numel(axs);
for i = 1:n
    if ~isvalid(axs(i)), continue; end
    isOuter = (strcmp(location, 'top') && i == 1) || (strcmp(location, 'bottom') && i == n);
    if isOuter
        axs(i).XAxisLocation = location;
        xlabel(axs(i), label);
    else
        axs(i).XTickLabel = [];
        xlabel(axs(i), '');
    end
end
end

function mute(ax)
% make plotted children non-pickable so the axis ButtonDownFcn fires on click
ch = allchild(ax);
if ~isempty(ch), set(ch, 'HitTest', 'off'); end
end

% EOF
