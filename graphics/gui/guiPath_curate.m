function [hFig, cfgData] = guiPath_curate(basepath, varargin)
% GUIPATH_CURATE General-purpose event-curation viewer (EDs, ripples, ...).
%
%   hFig = GUIPATH_CURATE(basepath, varargin)
%
%   SUMMARY:
%       Single-session signal viewer on uifigure, built on the shared gui_
%       helpers. It is modality-agnostic: an "event" is just a peak time (with
%       optional start/stop and an accepted flag), and any 1-D signal can be a
%       panel. What differs between modalities (which file to load, which
%       signals to show, the default layout, where to save) is captured by a
%       PRESET. A Preset dropdown (EDs, Ripples, ...) loads the relevant file
%       and applies that preset's default view; switching presets reloads.
%       Presets live in guiPath_presets.m (one entry per modality).
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
%       basepath    - (Char) Session directory. Defaults to pwd.
%       varargin    - Parameter/Value pairs:
%           'preset'   - (Char) Initial preset name (see guiPath_presets). If
%                        empty, auto-detected from the files present.
%           'cfgData'  - (Struct) Panels, flat: one field per panel (see
%                        guiPath_panel). Bypasses presets. See guiPath_doc.
%           'cfgGui'   - (Struct) Behaviour for a cfgData: .mode .win .save.
%                        Anything absent is derived from cfgData.
%           'basename' - (Char) Override (defaults to the folder name).
%           'Visible'  - (Char) 'on' (default) | 'off' (headless).
%
%   OUTPUT:
%       hFig        - (uifigure) Handle to the viewer window.
%       cfgData     - (Struct) The full panels (data loaded), for a fast reopen.
%                     The live version is also in hFig.UserData.cfgData.
%
%   DEPENDENCIES:
%       guiPath_presets, guiPath_load, guiPath_panel; gui_layout,
%       gui_labeledControl, gui_eventPanel, gui_notify, gui_chooseDialog;
%       plot_spec, plot_hypnogram, plot_raster.
%
%   HISTORY:
%       Created:  22 Jun 2026 (as ed_gui).
%       Renamed:  23 Jun 2026; ed/sSig dependency removed; Preset dropdown +
%                 preset registry; events are generic (peak + start/stop).
%       Renamed:  06 Jul 2026 -> guiPath_curate; moved to graphics/gui; the
%                 shared GUI package flattened to gui_* helpers.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
if nargin < 1 || isempty(basepath), basepath = pwd; end   % default to the current folder
p = inputParser;
addRequired(p, 'basepath', @ischar);
addParameter(p, 'preset', '', @ischar);
addParameter(p, 'cfgData', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'cfgGui', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'basename', '', @ischar);
addParameter(p, 'Visible', 'on', @(x) any(strcmpi(char(x), {'on', 'off'})));

parse(p, basepath, varargin{:});
basepath   = p.Results.basepath;
presetArg  = p.Results.preset;
cfgDataArg = p.Results.cfgData;
cfgGuiArg  = p.Results.cfgGui;
vis        = char(p.Results.Visible);

basename = p.Results.basename;
if isempty(basename)
    [~, basename] = fileparts(basepath);
end

%% ========================================================================
%  RESOLVE INITIAL CONFIG (preset or explicit inputs)
%  ========================================================================
presets = guiPath_presets();    % {name, file} list for dropdown + auto-detect

% pick the panels (cfgData) + behaviour (cfgGui): an explicit cfgData, else the
% named / auto-detected preset, else an empty template
if ~isempty(cfgDataArg)
    presetName = 'Custom';
    cfgData = cfgDataArg;
    if ~isempty(cfgGuiArg), cfgGui = cfgGuiArg; else, cfgGui = struct(); end
else
    presetName = resolvePreset(presetArg, presets, basepath, basename);
    if isempty(presetName), presetName = 'template'; end
    [cfgData, cfgGui] = guiPath_presets(presetName);
    presetName = cfgGui.name;
end

% load the data (skips panels that already carry it), then finalize behaviour and
% adapt to the internal render config. ctx carries basename + a within-session cache
ctx = guiPath_ctx(basepath, basename);
try
    cfgData = guiPath_load(cfgData, basepath, ctx);
catch ME
    warning('guiPath_curate:load', 'load failed (%s); opening empty.', ME.message);
    cfgData = struct();
end
cfgGui = finalizeGui(cfgGui, cfgData);
config = cfgDataToConfig(cfgData, cfgGui, basepath, basename);

%% ========================================================================
%  STATE
%  ========================================================================
d = struct();
d.basepath   = basepath;
d.basename   = basename;
d.presets    = presets;
d.presetName = presetName;
d.ctx        = ctx;            % basename + within-session file cache
d.cfgData    = cfgData;        % the live, full panels (returned for a fast reopen)
d.cfgGui     = cfgGui;
d.clrAccept  = [0.10 0.55 0.10];
d.clrReject  = [0.65 0.15 0.15];
d.unit       = struct('wide', 'hr', 'narrow', 's');   % per-region x-axis units
d.hInd       = gobjects(0);
d.hCenter    = gobjects(0);
d.dirty      = false;
d.win        = cfgGui.win;
d.winDefault = cfgGui.win;
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

[~, gPlot, gCtrl, gActions] = gui_layout(hFig, 'CtrlWidth', 240, 'CtrlSide', 'left');

% --- Controls: preset ---
gui_labeledControl(gCtrl, 'label', 'PRESET', 'FontWeight', 'bold');
d.hPresetDD = gui_labeledControl(gCtrl, 'dropdown', '', ...
    'Items', presetItems(presets, presetName), 'Value', presetName, ...
    'ValueChangedFcn', @(s, ~) loadPreset(s.Value));

% --- Controls: unified Load (opens a progressive dialog; nothing else lives
% here permanently). The dialog reveals inputs as you choose type + source. ---
gui_labeledControl(gCtrl, 'label', 'LOAD', 'FontWeight', 'bold');
gui_labeledControl(gCtrl, 'button', '', 'Text', 'Load...', 'ButtonPushedFcn', @(~, ~) onLoadUnified());

% --- Controls: layout configuration (per region) ---
if strcmp(d.cfgRegion, 'wide'), nInitCfg = numel(d.wideP); else, nInitCfg = numel(d.narrowP); end
gui_labeledControl(gCtrl, 'label', 'LAYOUT', 'FontWeight', 'bold');

% Region dropdown and # Panels side by side.
hRP = gui_labeledControl(gCtrl, 'panel', '', 'RowHeight', 56);
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

hSrcPanel  = gui_labeledControl(gCtrl, 'panel', 'Panels', 'RowHeight', 160);
d.hSrcGrid = uigridlayout(hSrcPanel, [1, 1], 'Padding', 2, 'RowSpacing', 3, ...
    'ColumnWidth', {'1x'}, 'Scrollable', 'on');

% --- Controls: view (t0 / window in seconds; per-region x-unit) ---
gui_labeledControl(gCtrl, 'label', 'VIEW', 'FontWeight', 'bold');
hVF = gui_labeledControl(gCtrl, 'panel', '', 'RowHeight', 84);
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

hZ = gui_labeledControl(gCtrl, 'panel', '', 'RowHeight', 30);
gZ = uigridlayout(hZ, [1, 2], 'ColumnWidth', {'1x', '1x'}, 'Padding', 2, 'ColumnSpacing', 4);
bZi = uibutton(gZ, 'Text', 'Zoom In (+)',  'ButtonPushedFcn', @(~,~) zoomActive(1/1.5));
bZi.Layout.Column = 1;
bZo = uibutton(gZ, 'Text', 'Zoom Out (-)', 'ButtonPushedFcn', @(~,~) zoomActive(1.5));
bZo.Layout.Column = 2;
gui_labeledControl(gCtrl, 'button', '', 'Text', 'Reset (0)', 'ButtonPushedFcn', @(~,~) resetActive());
gui_labeledControl(gCtrl, 'spacer', '');

% --- Action widget: gui_eventPanel (events) or gui_statePanel (states).
% Rebuilt by buildActionWidget when a preset switch changes the curation mode.
d.gActions = gActions;
d = buildActionWidget(d);

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
        prevMode = data.mode;
        if ~any(strcmp(name, {data.presets.name}))   % 'Custom'/'ws:...'/unknown
            data.hPresetDD.Value = data.presetName; return;
        end
        % offer to save unsaved curation before switching
        if data.dirty
            sel = gui_chooseDialog(hFig, 'Unsaved curation. Save before switching?', ...
                {'Save', 'Discard'});
            if isempty(sel), data.hPresetDD.Value = data.presetName; return; end
            if strcmp(sel, 'Save'), onSave(); data = hFig.UserData; end
        end
        dlg = busyOn(sprintf('Loading %s...', name));   % panels already loaded are reused
        try
            % reuse already-loaded panels from the live cfgData, load only the rest
            [newData, newGui] = guiPath_presets(name, data.cfgData);
            newData = guiPath_load(newData, data.basepath, data.ctx);
            newGui  = finalizeGui(newGui, newData);
            config  = cfgDataToConfig(newData, newGui, data.basepath, data.basename);
        catch ME
            busyOff(dlg);
            gui_notify(hFig, sprintf('Preset "%s" failed: %s', name, ME.message), 'error');
            data.hPresetDD.Value = data.presetName; hFig.UserData = data; return;
        end
        data.cfgData = newData; data.cfgGui = newGui;
        data = applyConfig(data, config);
        if ~strcmp(prevMode, data.mode)
            data = buildActionWidget(data);   % swap event <-> state widget
        end
        data.presetName = newGui.name;
        data.dirty = false;
        hFig.Name = sprintf('%s - Curate: %s', data.basename, data.presetName);
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

    function data = buildActionWidget(data)
        % (re)build the pinned action widget for the current mode:
        % gui_statePanel for states, gui_eventPanel for events. Removes prior.
        if isfield(data, 'ev') && ~isempty(data.ev) && ...
                isfield(data.ev, 'grid') && isvalid(data.ev.grid)
            delete(data.ev.grid);
        end
        if strcmp(data.mode, 'states')
            apiS = struct('prev', @() navStep(-1), 'next', @() navStep(1), ...
                'assign', @(k) assignState(k), 'save', @() onSave(), 'setIdx', @(v) setIdx(v));
            data.ev = gui_statePanel(data.gActions, apiS, data.stateNames, data.stateColors);
        else
            api = struct('prev', @() navStep(-1), 'next', @() navStep(1), ...
                'accept', @() setAccept(true), 'reject', @() setAccept(false), ...
                'save', @() onSave(), 'setIdx', @(v) setIdx(v));
            data.ev = gui_eventPanel(data.gActions, api);
        end
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
        prevMode = data.mode;
        data.mode = 'events';
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
        if ~strcmp(prevMode, data.mode)
            data = buildActionWidget(data);   % swap state -> event widget
        end
        hFig.UserData = data;
        rebuildPlot();
        populateSrc(hFig);
        refreshEvent(hFig);
    end

    function onLoadUnified()
        % open the progressive Load dialog, then integrate the chosen source via
        % loadCore (the same address machinery as the presets). The dialog owns
        % the type/source/placement choices; here we only wire the save target.
        data = hFig.UserData;
        sel = gui_loadDialog(hFig, data.basepath);
        if isempty(sel), return; end
        switch sel.from
            case 'ws',  src = ['ws:', sel.value];  label = sel.value;
            case 'bin', src = ['bin:', sel.value]; label = ['ch ', sel.value];
            otherwise,  src = sel.value; label = sel.var;      % file: inline value
        end
        saveFcn = [];
        if strcmp(sel.type, 'eventTicks')
            if strcmp(sel.from, 'ws')
                saveFcn = @(acc) assignin('base', matlab.lang.makeValidName([sel.var, '_accepted']), logical(acc(:)));
            elseif strcmp(sel.from, 'file')
                saveFcn = @(acc) saveVarToFile(sel.file, [sel.var, '_accepted'], logical(acc(:)));
            end
        elseif strcmp(sel.type, 'stateStrip')
            saveFcn = @(lb) saveLabelsFlow(fullfile(data.basepath, [data.basename, '.sleep_labelsMan.mat']), lb);
        end
        switch sel.type
            case 'eventTicks', loadNm = 'eventTicks';
            case 'stateStrip', loadNm = 'states';
            otherwise,         loadNm = matlab.lang.makeValidName(label);
        end
        % single-wrap src so an inline cell/array value (a File source) lands in
        % ONE struct, not a struct array (the classic struct(...,cell,...) trap)
        loadCore(struct('src', {src}, 'type', sel.type, 'region', sel.region, 'fs', sel.fs, ...
            'name', loadNm, 'label', label, 'saveFcn', saveFcn));
    end

    function loadCore(p)
        % load the one declared source (same path as presets) and integrate it: a
        % signal joins the pool + a new panel; events / states become the active
        % curation target. The loaded panel is also recorded in data.cfgData.
        data = hFig.UserData;
        ps = guiPath_panel(p.type, p.region, p.src, 'name', p.name, 'label', p.label, 'fs', p.fs);
        one = struct('item', ps);                    % a one-panel flat cfgData
        try
            one = guiPath_load(one, data.basepath, data.ctx);
        catch ME
            gui_notify(hFig, sprintf('Load failed: %s', ME.message), 'error'); return;
        end
        if isempty(fieldnames(one))                  % the panel was dropped (load failed)
            gui_notify(hFig, 'Could not load that source (check the value / type).', 'error'); return;
        end
        inp = recFromPanel(one.item);                % the loaded panel -> an input
        reg = regAlias(one.item.region);             % 'wide' / 'narrow'
        % capture what the message needs BEFORE the integrate call: 'inp' is a
        % shared nested-scope variable that rebuildPlot -> drawPanel reassigns
        itype = inp.type; iname = inp.name;
        % integrate + render inside a guard: a bad source must notify, never crash
        try
            switch itype
                case 'eventTicks'
                    nEv = numel(inp.data.peakTime);
                    setEvents(inp, p.saveFcn, p.label);
                    gui_notify(hFig, sprintf('Loaded %d events from "%s".', nEv, p.label), 'success');
                case 'stateStrip'
                    nLab = numel(inp.data.labels);
                    setStatesTarget(inp, p.saveFcn, p.label, data.winDefault, reg);
                    gui_notify(hFig, sprintf('Loaded %d state epochs from "%s".', nLab, p.label), 'success');
                otherwise
                    placeSignal(inp, reg);
                    gui_notify(hFig, sprintf('Loaded "%s" into the %s.', iname, regDisp(reg)), 'success');
            end
            d2 = hFig.UserData;                      % record the panel for the round-trip
            d2.cfgData.(matlab.lang.makeValidName(iname)) = one.item;
            hFig.UserData = d2;
        catch ME
            gui_notify(hFig, sprintf('Could not display "%s": %s', iname, ME.message), 'error');
        end
    end

    function placeSignal(inp, reg)
        % add a signal input to the persistent pool (unique name) and open a panel
        % for it in the chosen region, preserving any active curation target
        data = hFig.UserData;
        base = inp.name; nm = base; k = 1;
        while ~isempty(data.pool) && any(strcmp(nm, {data.pool.name}))
            k = k + 1; nm = sprintf('%s_%d', base, k);
        end
        inp.name = nm; inp.defRegion = reg;
        data.pool = mergePool(data.pool, inp);
        evIx = find(ismember({data.inputs.type}, {'eventTicks', 'stateStrip'}), 1);
        if ~isempty(evIx), data.inputs = [data.pool, data.inputs(evIx)]; else, data.inputs = data.pool; end
        data.Tend_s = max(eps, computeTend(data.inputs));
        pn = makePanel(nm, reg, data.inputs);
        if strcmp(reg, 'wide'), data.wideP = [data.wideP, pn]; else, data.narrowP = [data.narrowP, pn]; end
        hFig.UserData = data;
        rebuildPlot();
        populateSrc(hFig);
    end

    function setStatesTarget(inp, sv, label, winPlot, reg)
        % make a stateStrip input the active curation target (states mode); the
        % sibling of setEvents. Adds one strip panel in the chosen region.
        data = hFig.UserData;
        prevMode = data.mode;
        data.mode = 'states';
        strip = normalizeInputs(inp);
        data.inputs = [data.pool, strip(1)];
        data.saveFcn = sv;
        data.presetName = label;
        data.dirty = false;
        if ~isempty(winPlot), data.win = winPlot; data.winDefault = winPlot; end
        data.Tend_s = max(eps, computeTend(data.inputs));
        [data.labels, epochT, data.epochLen, data.nstates, data.stateNames, data.stateColors] = extractStates(data.inputs);
        nEp = numel(data.labels);
        data.ed = struct('peakTime', epochT(:));
        data.accepted = true(max(0, nEp), 1);
        data.nEvents = nEp; data.hasEvents = nEp > 0;
        data.currIdx = 1;
        if data.hasEvents && data.nEvents > 0, data.t0 = data.ed.peakTime(1); else, data.t0 = data.Tend_s / 2; end
        if ~any(strcmp({data.wideP.source}, strip(1).name)) && ~any(strcmp({data.narrowP.source}, strip(1).name))
            pn = makePanel(strip(1).name, reg, data.inputs);
            if strcmp(reg, 'wide'), data.wideP = [data.wideP, pn]; else, data.narrowP = [data.narrowP, pn]; end
        end
        items = data.hPresetDD.Items;
        if ~any(strcmp(label, items)), data.hPresetDD.Items = [items, {label}]; end
        data.hPresetDD.Value = label;
        hFig.Name = sprintf('%s - Curate: %s', data.basename, label);
        if ~strcmp(prevMode, data.mode)
            data = buildActionWidget(data);          % swap event -> state widget
        end
        hFig.UserData = data;
        rebuildPlot();
        populateSrc(hFig);
        refreshEvent(hFig);
    end

    function saveLabelsFlow(file, labels)
        % Load-flow states save: AccuSleep-compatible labels vector, backed up
        backup_file(file);
        labels = labels(:);
        save(file, 'labels');
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
            [data.hInd(i), data.hCenter(i)] = addCursor(ax);
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
                [data.hInd(i), data.hCenter(i)] = addCursor(ax);
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
                xline(ax, eventMarks(data.ed, ev) / xf, '--b');
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
            case 'stateStrip', drawStateStrip(ax, data, xf);
            case 'raster',     drawRaster(ax, inp, a, b, xf);
        end
    end

    function [hI, hC] = addCursor(ax)
        % overview cursor line (t0) + window band. Left with default
        % HandleVisibility so the next cla clears them (no stale-overlay pileup);
        % updateMarker repositions them in place between redraws.
        hI = xregion(ax, 0, eps, 'FaceColor', [0 0.2 0.8], 'FaceAlpha', 0.15);
        hC = xline(ax, 0, 'Color', [0 0.25 0.9], 'LineWidth', 1.25);
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
        % sync the action widget: index / total (+ current state in states mode).
        % accept-reject is read off the plot panels, so no status text there
        data = fig.UserData;
        if isempty(data.ev), return; end
        if strcmp(data.mode, 'states')
            if data.nEvents == 0, data.ev.refresh(0, 0, []); return; end
            data.ev.refresh(data.currIdx, data.nEvents, data.labels(data.currIdx));
            return;
        end
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
        % step to the next (step>0) or previous (step<0) event RELATIVE TO the
        % current view t0. When t0 sits on an event this is simply +/-1; after a
        % free click it continues from wherever the view now is. Events are
        % chronological (peakTime ascending), so index order tracks time order.
        data = hFig.UserData;
        if ~data.hasEvents || data.nEvents == 0, return; end
        pk = data.ed.peakTime;
        if step > 0
            idx = find(pk > data.t0, 1, 'first');
            if isempty(idx), idx = data.nEvents; end        % already past the last
        else
            idx = find(pk < data.t0, 1, 'last');
            if isempty(idx), idx = 1; end                   % already before the first
        end
        data.currIdx = idx;
        data.t0 = pk(idx);
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
        % click / t0 edit: move the Bottom to tSec and re-anchor the event cursor
        % so arrow-stepping continues from here. States snap the view onto the
        % epoch (the epoch IS the position); events keep t0 where you clicked but
        % point the cursor at the nearest event.
        data = hFig.UserData;
        tSec = min(max(0, tSec), data.Tend_s);
        if data.hasEvents && data.nEvents > 0
            [~, idx] = min(abs(data.ed.peakTime - tSec));
            data.currIdx = idx;
            if strcmp(data.mode, 'states'), tSec = data.ed.peakTime(idx); end
        end
        data.t0 = tSec;
        hFig.UserData = data;
        renderNarrow(hFig); updateMarker(hFig); refreshEvent(hFig);
    end

    function setAccept(tf)
        data = hFig.UserData;
        if ~strcmp(data.mode, 'events'), return; end
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

    function assignState(stateVal)
        % states mode: set the current epoch's label, recolour, auto-advance
        data = hFig.UserData;
        if ~strcmp(data.mode, 'states') || data.nEvents == 0, return; end
        stateVal = round(stateVal);
        if stateVal < 1 || stateVal > data.nstates + 1, return; end   % N+1 = undefined
        data.labels(data.currIdx) = stateVal;
        data.dirty = true;
        hFig.UserData = data;
        renderStateStrip(hFig);           % recolour Top state strip
        if data.currIdx < data.nEvents
            navStep(1);                   % auto-advance (re-renders Bottom + markers)
        else
            renderNarrow(hFig); refreshEvent(hFig);
        end
    end

    function jumpNextUndefined()
        % states mode: jump to the next epoch whose label is undefined (> N)
        data = hFig.UserData;
        if ~strcmp(data.mode, 'states') || data.nEvents == 0, return; end
        rel = find(data.labels((data.currIdx + 1):end) > data.nstates, 1);
        if isempty(rel)
            nxt = find(data.labels > data.nstates, 1);      % wrap to first
            if isempty(nxt), gui_notify(hFig, 'No undefined epochs.', 'info'); return; end
        else
            nxt = data.currIdx + rel;
        end
        setIdx(nxt);
    end

    function renderStateStrip(fig)
        % recolour the Top state strip after a label change (mirror of
        % renderTicksPanel): redraw the strip, restore its cursor + window
        % markers; the Bottom is refreshed by the auto-advance / caller
        data = fig.UserData;
        xf = unitSec(data.unit.wide);
        for i = 1:numel(data.wideP)
            if isStrip(data, data.wideP(i).source)
                ax = data.wideP(i).ax;
                cla(ax); hold(ax, 'on');
                drawStateStrip(ax, data, xf);
                ylabel(ax, data.wideP(i).label);
                set(get(ax, 'YLabel'), 'Color', [0.15 0.15 0.15]);
                ax.XLim = [0, data.Tend_s] / xf;
                [data.hInd(i), data.hCenter(i)] = addCursor(ax);
                mute(ax);
            end
        end
        fig.UserData = data;
        hideOuterXTicks(data.axWide, ['Time (' data.unit.wide ')'], 'bottom');
        updateMarker(fig);
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

    function tf = isEditingField()
        % true while a numeric/text edit field holds focus, so the figure-level
        % keyboard shortcuts stand down and the field receives the keystrokes
        co = hFig.CurrentObject;
        tf = ~isempty(co) && ...
            (isa(co, 'matlab.ui.control.NumericEditField') || ...
             isa(co, 'matlab.ui.control.EditField'));
    end

    function onKey(~, evt)
        if isEditingField(), return; end   % typing in a field: let it keep the key
        if any(strcmpi(evt.Modifier, 'control')) && strcmpi(evt.Key, 's')
            onSave(); return;
        end
        data = hFig.UserData;
        if strcmp(data.mode, 'states')
            switch evt.Key
                case {'rightarrow', 'uparrow'},   navStep(1);
                case {'leftarrow', 'downarrow'},  navStep(-1);
                case {'1', '2', '3', '4', '5', '6', '7', '8', '9'}
                    assignState(str2double(evt.Key));
                case {'numpad1', 'numpad2', 'numpad3', 'numpad4', 'numpad5', ...
                        'numpad6', 'numpad7', 'numpad8', 'numpad9'}
                    assignState(str2double(evt.Key(end)));
                case 'x',                    assignState(data.nstates + 1);
                case 'n',                    jumpNextUndefined();
                case {'equal', 'add'},       zoomActive(1 / 1.5);
                case {'hyphen', 'subtract'}, zoomActive(1.5);
                case {'0', 'numpad0'},       resetActive();
            end
            return;
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
        data.loadCoreFcn   = @loadCore;         % exposed for tests (drives the Load flow)
        data.loadPresetFcn = @loadPreset;       % exposed for tests (drives a preset switch)
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
            if isempty(opts)
                gui_notify(hFig, 'Load a signal or events first (no sources to show).', 'warning');
                src.Value = max(1, cur); return;
            end
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
        if isempty(data.saveFcn)
            gui_notify(hFig, 'No save target for this view.', 'warning');
            return;
        end
        if strcmp(data.mode, 'states')
            if data.nEvents == 0, return; end
            data.saveFcn(data.labels(:));
            data.dirty = false;
            hFig.UserData = data;
            gui_notify(hFig, sprintf('Saved %d epoch labels.', data.nEvents), 'success');
            return;
        end
        if ~data.hasEvents, return; end
        data.saveFcn(data.accepted(:));
        data.dirty = false;
        hFig.UserData = data;
        gui_notify(hFig, sprintf('Saved %d accepted / %d events.', ...
            sum(data.accepted), data.nEvents), 'success');
    end

    function onClose()
        data = hFig.UserData;
        if isfield(data, 'dirty') && data.dirty
            sel = gui_chooseDialog(hFig, 'Unsaved curation. Save before closing?', ...
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
% chosen preset, else the first preset whose file exists, else '' (empty ->
% the caller opens an empty Custom view; guiPath_curate never requires a file)
names = {presets.name};
if ~isempty(arg)
    ix = find(strcmpi(arg, names), 1);
    if ~isempty(ix), name = names{ix}; return; end
end
for i = 1:numel(presets)
    if isfile(fullfile(basepath, [basename, '.', presets(i).file, '.mat']))
        name = presets(i).name; return;
    end
end
name = '';
end

function g = finalizeGui(g, cfgData)
% fill any cfgGui field the caller left out: mode from the target panel, window
% from the mode, name/file/save to safe defaults
if ~isfield(g, 'name') || isempty(g.name), g.name = 'Custom'; end
if ~isfield(g, 'file'), g.file = ''; end
if ~isfield(g, 'mode') || isempty(g.mode)
    g.mode = 'events';
    fns = fieldnames(cfgData);
    for i = 1:numel(fns)
        if strcmp(cfgData.(fns{i}).type, 'stateStrip'), g.mode = 'states'; break; end
    end
end
if ~isfield(g, 'win') || isempty(g.win)
    if strcmp(g.mode, 'states'), g.win = 10; else, g.win = 1; end
end
if ~isfield(g, 'save'), g.save = ''; end
end

function fn = resolveSave(save, mode, basepath, basename, cfgData)
% turn a cfgGui.save spec into a save handle: a function handle is used as is;
% states write labels to sleep_labelsMan; an events token writes the accepted
% mask to <basename>.<token>.mat (or a 'ws:VAR' target to a base variable)
if isa(save, 'function_handle'), fn = save; return; end
if strcmp(mode, 'states')
    fn = @(labels) saveLabels(fullfile(basepath, [basename, '.sleep_labelsMan.mat']), labels);
    return;
end
if ischar(save) && ~isempty(save) && ~strcmp(save, 'labelsMan'), tok = save; else, tok = targetToken(cfgData); end
fn = [];
if isempty(tok), return; end
if any(tok == ':')
    ci = find(tok == ':', 1);
    if strcmp(tok(1:ci - 1), 'ws')
        vn = tok(ci + 1:end);
        fn = @(acc) assignin('base', matlab.lang.makeValidName([vn, '_accepted']), logical(acc(:)));
    end
    return;
end
file = fullfile(basepath, [basename, '.', tok, '.mat']);
if isfile(file), fn = @(acc) saveAccepted(file, tok, acc); end
end

function tok = targetToken(cfgData)
% the src token of the eventTicks target ('ed' / 'ripp' / a 'ws:VAR' address)
tok = '';
fns = fieldnames(cfgData);
for i = 1:numel(fns)
    p = cfgData.(fns{i});
    if strcmp(p.type, 'eventTicks') && ischar(p.src) && ~isempty(p.src), tok = p.src; return; end
end
end

function saveAccepted(file, varName, accepted)
% back up the file, then set <var>.accepted and write all variables back
backup_file(file);
S = load(file);
S.(varName).accepted = logical(accepted(:));
save(file, '-struct', 'S', '-v7.3');
end

function saveLabels(file, labels)
% AccuSleep-compatible manual labels: back up any existing file, write the vector
backup_file(file);
labels = labels(:);
save(file, 'labels');
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
% the curation TARGET (eventTicks or stateStrip) is the active set; everything
% else is a signal that goes to the persistent pool
isEvt = strcmp({allIn.type}, 'eventTicks') | strcmp({allIn.type}, 'stateStrip');
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

d.mode = 'events';
if isfield(config, 'mode') && ~isempty(config.mode), d.mode = config.mode; end
if strcmp(d.mode, 'states')
    % epochs navigate as "events" (t0 steps by epoch); labels are the target
    [d.labels, epochT, d.epochLen, d.nstates, d.stateNames, d.stateColors] = extractStates(d.inputs);
    nEp = numel(d.labels);
    d.ed = struct('peakTime', epochT(:));
    d.accepted = true(max(0, nEp), 1);      % unused in states mode
    d.nEvents = nEp; d.hasEvents = nEp > 0;
else
    [d.ed, d.accepted, d.nEvents, d.hasEvents] = extractEvents(d.inputs);
    d.labels = []; d.epochLen = 1; d.nstates = 0; d.stateNames = {}; d.stateColors = {};
end
d.currIdx = 1;
if d.hasEvents && d.nEvents > 0, d.t0 = d.ed.peakTime(1); else, d.t0 = d.Tend_s / 2; end

panels = normalizePanels(config, d.inputs);
d.wideP   = panels(strcmp({panels.region}, 'wide'));
d.narrowP = panels(strcmp({panels.region}, 'narrow'));
if isempty(d.wideP) && ~isempty(d.narrowP), d.cfgRegion = 'narrow'; else, d.cfgRegion = 'wide'; end
end

function config = cfgDataToConfig(cfgData, cfgGui, basepath, basename)
% adapt the flat, full cfgData into the internal render config: inputs deduped by
% panel name (a top+bottom pair shares one loaded input) + a {source,region} panel
% per field in stacking order. Behaviour (mode / win) comes from cfgGui; the save
% handle is resolved from cfgGui.save + the session location.
fns = fieldnames(cfgData);
recs = {}; names = {}; P = {};
for i = 1:numel(fns)
    pc = cfgData.(fns{i});
    nm = pick(pc, 'name', fns{i});
    if ~any(strcmp(nm, names))
        recs{end + 1} = panelToRec(pc, nm, i);                        %#ok<AGROW>
        names{end + 1} = nm;                                          %#ok<AGROW>
    end
    P{end + 1} = struct('source', nm, 'region', regAlias(pc.region)); %#ok<AGROW>
end
if isempty(recs), inputs = normalizeInputs([]); else, inputs = normalizeInputs([recs{:}]); end
if isempty(P), panels = []; else, panels = [P{:}]; end
saveFcn = resolveSave(cfgGui.save, cfgGui.mode, basepath, basename, cfgData);
config = struct('inputs', {inputs}, 'panels', {panels}, ...
    'saveFcn', saveFcn, 'winPlot', cfgGui.win, 'mode', cfgGui.mode);
end

function inp = recFromPanel(pc)
% one normalized input from a single enriched panel (the Load-dialog path)
inp = normalizeInputs(panelToRec(pc, pick(pc, 'name', 'item'), 99));
end

function rec = panelToRec(pc, nm, order)
% a pre-normalization input record from an enriched panel. data single-wrapped so
% a cell payload (raster / hypnogram) stays in one field.
rec = struct('name', nm, 'type', pc.type, 'data', {pick(pc, 'data', [])}, ...
    'fs', pick(pc, 'fs', NaN), 'ylim', pick(pc, 'ylim', []), 'clr', pick(pc, 'clr', 'k'), ...
    'label', pick(pc, 'label', ''), 'height', pick(pc, 'height', 1), ...
    'defRegion', regAlias(pc.region), 'defOrder', order);
end

function r = regAlias(region)
% user-facing top/bottom -> internal wide/narrow (both accepted)
switch lower(region)
    case {'top', 'wide'},      r = 'wide';
    case {'bottom', 'narrow'}, r = 'narrow';
    otherwise, r = region;
end
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

function [labels, epochT, epochLen, nstates, names, colors] = extractStates(inputs)
% pull the stateStrip input's data into the live state-scoring state
labels = []; epochT = []; epochLen = 1; nstates = 0; names = {}; colors = {};
ix = find(strcmp({inputs.type}, 'stateStrip'), 1);
if isempty(ix), return; end
D = inputs(ix).data;
if ~isstruct(D) || ~isfield(D, 'labels') || isempty(D.labels), return; end
labels = double(D.labels(:));
nEp = numel(labels);
if isfield(D, 'epochT') && numel(D.epochT) == nEp
    epochT = double(D.epochT(:));
else
    epochT = (0:nEp - 1)';
end
if nEp > 1, epochLen = median(diff(epochT)); end
if isfield(D, 'names')  && ~isempty(D.names),  names  = D.names;  end
if isfield(D, 'colors') && ~isempty(D.colors), colors = D.colors; end
if isfield(D, 'nstates') && ~isempty(D.nstates)
    nstates = D.nstates;
else
    nstates = max(1, numel(names));
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

function saveVarToFile(file, varName, value)
% back up the file (if present) then write value under varName, preserving any
% other variables already in the file
backup_file(file);
tmp = struct(varName, value);
if isfile(file)
    save(file, '-struct', 'tmp', '-append');
else
    save(file, '-struct', 'tmp');
end
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
if s2 < s1, return; end          % window outside this signal's extent -> nothing to draw
rng = s1:s2;
t = ((rng - 1) / fs) / xf;
np = numel(rng); maxPts = 20000;
if np > maxPts
    st = ceil(np / maxPts);
    plot(ax, t(1:st:end), sig(rng(1:st:end)), 'Color', inp.clr);
else
    plot(ax, t, sig(rng), 'Color', inp.clr);
end
% only a valid, increasing, finite range (a flat / NaN signal gives lo==hi)
if numel(inp.ylim) == 2 && all(isfinite(inp.ylim)) && inp.ylim(2) > inp.ylim(1)
    ax.YLim = inp.ylim;
end
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
line(ax, X(:), Y(:), 'Color', clr);
end

function drawStateStrip(ax, data, xf)
% per-epoch coloured label strip drawn as one truecolor image (fast to redraw).
% epoch centres come from the navigation vector (ed.peakTime), labels from
% data.labels; undefined (> nstates) render gray.
if isempty(data.labels) || ~isfield(data.ed, 'peakTime'), return; end
T = data.ed.peakTime(:)';
L = data.labels(:)';
n = min(numel(T), numel(L));
if n == 0, return; end
T = T(1:n); L = L(1:n);
cmap = stateCmap(data);
Lc = min(max(round(L), 1), size(cmap, 1));
cdata = reshape(cmap(Lc, :), [1, n, 3]);
if n == 1, xl = [T(1) - 0.5, T(1) + 0.5]; else, xl = [T(1), T(end)]; end
image(ax, 'XData', xl / xf, 'YData', [0, 1], 'CData', cdata);
ax.YLim = [0, 1]; ax.YTick = [];
end

function cmap = stateCmap(data)
% (nstates+1) x 3 state colormap; the extra row (undefined) is gray
ns = max(1, data.nstates);
cmap = repmat([0.6 0.6 0.6], ns + 1, 1);
C = data.stateColors;
for i = 1:min(ns, numel(C))
    c = C{i};
    if numel(c) >= 3, cmap(i, :) = c(1:3); end
end
end

function tf = isStrip(data, source)
% true if the panel source names a stateStrip input
inp = getInput(data.inputs, source);
tf = ~isempty(inp) && strcmp(inp.type, 'stateStrip');
end

%% ========================================================================
%  INPUTS MODEL (pure)
%  ========================================================================

function b = blankInput()
b = struct('name', '', 'type', '', 'data', [], 'fs', NaN, 'ylim', [], ...
    'clr', 'k', 'label', '', 'height', 1, 'defRegion', 'narrow', 'defOrder', 99);
end

function out = normalizeInputs(in)
% fill defaults on an inputs struct array (requires .name .type). Per-type
% appearance is set upstream by guiPath_panel, so these are generic fallbacks only.
n = numel(in);
out = repmat(blankInput(), 1, n);
for i = 1:n
    s = in(i);
    if ~isfield(s, 'name') || isempty(s.name), error('guiPath_curate:inputs', 'input %d needs a name', i); end
    if ~isfield(s, 'type') || isempty(s.type), error('guiPath_curate:inputs', 'input %d needs a type', i); end
    out(i).name      = s.name;
    out(i).type      = s.type;
    out(i).data      = pick(s, 'data', []);
    out(i).fs        = pick(s, 'fs', NaN);
    out(i).ylim      = pick(s, 'ylim', []);
    out(i).clr       = pick(s, 'clr', 'k');
    out(i).label     = pick(s, 'label', '');
    if isempty(out(i).label), out(i).label = s.name; end
    out(i).height    = pick(s, 'height', 1);
    out(i).defRegion = pick(s, 'defRegion', 'narrow');
    out(i).defOrder  = pick(s, 'defOrder', 99);
end
end

function v = pick(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
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
        case 'stateStrip'
            if isstruct(inp.data) && isfield(inp.data, 'epochT') && ~isempty(inp.data.epochT)
                T = max(T, max(inp.data.epochT(:)));
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
