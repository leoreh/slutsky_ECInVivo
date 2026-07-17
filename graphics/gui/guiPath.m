function [hFig, cfgData] = guiPath(basepath, varargin)

% Opens a single-session signal viewer for curating events or scoring states.
%
% Modality-agnostic: an "event" is just a peak time (with optional start / stop
% and an accepted flag), and any 1-D signal can be a panel. Two things are
% chosen independently, by two dropdowns:
%   PRESET  the arrangement - which file to load, which signals to show, the
%           default layout - one per modality in guiPath_presets.
%   CURATE  which loaded event / state set is the editable target, or None (view
%           only: Prev / Next then step the window itself). Several event sets
%           can be loaded at once (Load...); everything loaded stays available
%           across preset switches. Every set is drawn - Top ticks and Bottom
%           window marks, each its own colour - and CURATE picks the one you
%           step through and accept / reject. The viewer auto-curates the
%           starting preset's set on open; switching preset later changes the
%           arrangement only, not the target.
%
% Panels stack in two regions separated by a thin divider:
%       Top (wide)      full-session overview; drawn once, x-range zoomable.
%       Bottom (narrow) a window around the cursor t0; redrawn as t0 moves.
% All panels share one tiledlayout, so every plot box shares a left gutter and
% width. The Top prints its x-axis just above the divider, the Bottom at the
% very bottom. Panel types, addresses and the cfgData contract are documented
% in guiPath_doc; the draw function per type lives in guiPath_draw.
%
% INTERACTION
% - Left / Right            previous / next event (or the window, in view mode)
% - Up / Down               accept / reject, then auto-advance
% - 1-9 / x                 assign a state / undefined  (states mode)
% - n                       jump to the next undefined epoch  (states mode)
% - Ctrl+S                  save, via the preset's save function
% - + / - / 0               zoom the active region in / out / reset (time)
% - Shift + / - / 0         amplitude of the active panel: bigger / smaller /
%                           reset. Per type: a trace's y-limits, a stack's
%                           per-channel gain, a spectrogram's brightness.
% - Shift + scroll          the same, for the panel under the pointer
% - click any panel         activate its region and move t0 (so the Bottom)
%                           to the clicked time; the panel becomes the target
%                           for Shift +/-/0
% - scroll over a region    zoom that region's time axis, about the pointer
%
% EXAMPLES
% - guiPath(basepath)
%   auto-detects the preset from the files present in basepath.
%
% - guiPath(basepath, 'preset', 'Ripples')
%   forces a preset by name.
%
% - [hFig, cfgData] = guiPath(basepath, 'preset', 'EDs');
%   guiPath(basepath, 'cfgData', cfgData)
%   reopens fast: the returned cfgData already carries the loaded data.
%
% INPUTS
% - basepath        <char>(opt) session directory. {pwd}
% - preset          <char>(opt) initial preset name, see guiPath_presets.
%                   Auto-detected from the files present when empty.
% - cfgData         <struct>(opt) panels, flat: one field per panel, see
%                   guiPath_panel. Bypasses presets. See guiPath_doc.
% - cfgGui          <struct>(opt) behaviour for a cfgData: .mode .win .save.
%                   Anything absent is derived from cfgData.
% - basename        <char>(opt) override. {the folder name}
% - Visible         <char>(opt) 'on' | 'off' for headless. {'on'}
%
% OUTPUTS
% - hFig            <uifigure> handle to the viewer window.
% - cfgData         <struct> the full panels, data loaded, for a fast reopen.
%                   The live version is also in hFig.UserData.cfgData.
%
% SEE ALSO
% - guiPath_doc
% - guiPath_presets
% - guiPath_panel
% - guiPath_load
% - guiPath_src
% - guiPath_draw
%
% HISTORY
% - 260622          created as ed_gui.
% - 260623          renamed; ed / sSig dependency dropped; Preset dropdown and
%                   registry; events are generic (peak + start / stop).
% - 260706          renamed to guiPath_curate; moved to graphics/gui; the
%                   shared GUI package flattened to gui_* helpers.
% - 260716          rejected events are drawn, not hidden; a drifted .accepted
%                   warns rather than silently resetting; the event stepper
%                   shows the current event's vigilance state.
% - 260716          renamed to guiPath, the family's entry point; the draw
%                   functions extracted to guiPath_draw.
% - 260716          a click activates its region (Region dropdown follows);
%                   shift+scroll and shift+/-/0 adjust a panel's amplitude
%                   (per-panel yAdjust, applied in guiPath_draw).
% - 260716          one blue for every event mark (top ticks + Bottom lines),
%                   the Bottom peak solid and start / stop dashed; marker strips
%                   get a horizontal y-label and no y-ticks, but still carry the
%                   region's Time axis when they are its bottom panel (the event
%                   strip is last, so it shows the hours); the # Panels field is
%                   replaced by an add / delete / reorder panel list.
% - 260717          curation split from arrangement: event / state sets are
%                   plural inputs, a CURATE dropdown elects the editable target
%                   (or None); every set draws at once, each its own colour;
%                   Bottom window shows all sets' in-window marks. The Region
%                   dropdown is gone (a click sets the region; the Panels header
%                   shows it). PRESET now means the arrangement only.
% - 260717          None is a real "view" mode (gui_viewPanel): the accept /
%                   reject widget is replaced by a Prev / Next that steps the
%                   window itself (20% overlap). Tick strips draw from their
%                   own set, not the curated target, so they stay filled at
%                   None. d.inputs ACCUMULATES across preset switches: a preset
%                   changes the arrangement and loads its data but never drops
%                   what is loaded, so the panel list and CURATE keep every
%                   loaded signal / set; the target is set on first open only.
%                   PRESET (arrangement) and CURATE (target) are independent.

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
    warning('guiPath:load', 'load failed (%s); opening empty.', ME.message);
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
d.clrEvt     = [0 0.2 0.8];    % one blue for every event mark (top ticks + Bottom lines)
d.unit       = struct('wide', 'hr', 'narrow', 's');   % per-region x-axis units
d.hInd       = gobjects(0);
d.hCenter    = gobjects(0);
d.dirty      = false;
d.modShift   = false;         % Shift held? (tracked for shift+scroll amplitude)
d.activeSrc  = '';            % last-clicked panel source (shift+/-/0 target)
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
    'WindowKeyPressFcn', @onKey, 'WindowKeyReleaseFcn', @onKeyRelease, ...
    'WindowScrollWheelFcn', @onScroll);

[~, gPlot, gCtrl, gActions] = gui_layout(hFig, 'CtrlWidth', 240, 'CtrlSide', 'left');

% --- Controls: PRESET (arrangement) + CURATE (target). Two label + dropdown
% rows. PRESET loads a modality's signal set + panel layout; CURATE picks which
% loaded event / state set is the editable target (None = view only). The two
% are independent: you can curate one set while other sets stay visible. ---
hPC = gui_labeledControl(gCtrl, 'panel', '', 'RowHeight', 64);
gPC = uigridlayout(hPC, [2, 2], 'RowHeight', {'fit', 'fit'}, ...
    'ColumnWidth', {'fit', '1x'}, 'Padding', 2, 'RowSpacing', 4, 'ColumnSpacing', 6);
lblP = uilabel(gPC, 'Text', 'PRESET', 'FontWeight', 'bold');
lblP.Layout.Row = 1; lblP.Layout.Column = 1;
d.hPresetDD = uidropdown(gPC, 'Items', presetItems(presets, presetName), ...
    'Value', presetName, 'ValueChangedFcn', @(s, ~) loadPreset(s.Value));
d.hPresetDD.Layout.Row = 1; d.hPresetDD.Layout.Column = 2;
lblC = uilabel(gPC, 'Text', 'CURATE', 'FontWeight', 'bold');
lblC.Layout.Row = 2; lblC.Layout.Column = 1;
d.hCurateDD = uidropdown(gPC, 'Items', curateItems(d), ...
    'Value', curateDisp(d), 'ValueChangedFcn', @(s, ~) onCurateSel(s.Value));
d.hCurateDD.Layout.Row = 2; d.hCurateDD.Layout.Column = 2;

% --- Controls: unified Load (opens a progressive dialog; nothing else lives
% here permanently). The dialog reveals inputs as you choose type + source. ---
gui_labeledControl(gCtrl, 'label', 'LOAD', 'FontWeight', 'bold');
gui_labeledControl(gCtrl, 'button', '', 'Text', 'Load...', 'ButtonPushedFcn', @(~, ~) onLoadUnified());

% --- Controls: panel list for the active region. The region is set by clicking
% a panel (no dropdown); the header shows which one is being edited. Each row is
% a source dropdown with up / down (reorder) and X (delete); an Add row appends
% one. There is no panel-count field - the list IS the layout. ---
d.hPanelsLbl = gui_labeledControl(gCtrl, 'label', ...
    ['Panels (' regDisp(d.cfgRegion) ')'], 'FontWeight', 'bold');
hSrcPanel  = gui_labeledControl(gCtrl, 'panel', '', 'RowHeight', 180);
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
updateTitle(hFig);

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
        data = applyConfig(data, config);     % merges inputs, keeps the curate target
        if ~strcmp(prevMode, data.mode)
            data = buildActionWidget(data);   % swap event <-> state <-> view widget
        end
        data.presetName = newGui.name;
        data.dirty = false;
        hFig.UserData = data;
        updateTitle(hFig);
        refreshCurateDD(hFig);
        data.hPanelsLbl.Text = ['Panels (' regDisp(data.cfgRegion) ')'];
        data.hUnit.Value     = data.unit.(data.cfgRegion);
        data.hWin.Value      = data.win;
        rebuildPlot();
        populateSrc(hFig);
        refreshEvent(hFig);
        busyOff(dlg);
    end

    function data = buildActionWidget(data)
        % (re)build the pinned action widget for the current mode: gui_statePanel
        % (states), gui_viewPanel (view / None), else gui_eventPanel. Removes prior.
        if isfield(data, 'ev') && ~isempty(data.ev) && ...
                isfield(data.ev, 'grid') && isvalid(data.ev.grid)
            delete(data.ev.grid);
        end
        if strcmp(data.mode, 'states')
            apiS = struct('prev', @() navStep(-1), 'next', @() navStep(1), ...
                'assign', @(k) assignState(k), 'save', @() onSave(), 'setIdx', @(v) setIdx(v));
            data.ev = gui_statePanel(data.gActions, apiS, data.stateNames, data.stateColors);
        elseif strcmp(data.mode, 'view')
            apiV = struct('prev', @() navWindow(-1), 'next', @() navWindow(1));
            data.ev = gui_viewPanel(data.gActions, apiV);
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

    % --- CURATE: elect which loaded set is the editable target (or None).
    % refreshCurateDD keeps the selector's items in sync with the loaded sets;
    % onCurateSel guards an unsaved switch, setCurate does the work. ---
    function onCurateSel(val)
        data = hFig.UserData;
        target = ''; if ~strcmp(val, 'None'), target = val; end
        if strcmp(target, data.curate), return; end
        if data.dirty
            sel = gui_chooseDialog(hFig, 'Unsaved curation. Save before switching?', ...
                {'Save', 'Discard'});
            if isempty(sel), data.hCurateDD.Value = curateDisp(data); return; end
            if strcmp(sel, 'Save'), onSave(); end
        end
        setCurate(target);
    end

    function setCurate(target)
        % switch the editable target to the set named TARGET ('' = None): write
        % the outgoing set's edits back into its input, mirror the incoming set
        % into the live working state, swap the action widget if the mode flips,
        % then recolour the strips and redraw the window.
        data = hFig.UserData;
        data = writeBackCurate(data);
        prevMode = data.mode;
        data.curate  = target;
        data.currIdx = 1;
        data = syncCurate(data);
        if data.hasEvents && data.nEvents > 0, data.t0 = data.ed.peakTime(1); end
        data.dirty = false;
        if ~strcmp(prevMode, data.mode)
            data = buildActionWidget(data);   % events <-> states <-> view
        end
        hFig.UserData = data;
        updateTitle(hFig);
        refreshCurateDD(hFig);
        renderCurationStrips(hFig);           % recolour outgoing + incoming strips
        renderNarrow(hFig);                   % window marks + emphasis
        updateMarker(hFig);
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
        % name the set from its source so it reads clearly in CURATE (a collision
        % with an existing set is made unique in addCurationSet)
        loadNm = matlab.lang.makeValidName(label);
        % single-wrap src so an inline cell/array value (a File source) lands in
        % ONE struct, not a struct array (the classic struct(...,cell,...) trap)
        loadCore(struct('src', {src}, 'type', sel.type, 'region', sel.region, 'fs', sel.fs, ...
            'name', loadNm, 'label', label, 'saveFcn', saveFcn));
    end

    function loadCore(p)
        % load the one declared source (same path as presets) and integrate it: a
        % signal joins d.inputs + a new panel; an event / state set joins as a Top
        % strip and is listed in CURATE (curated only if none is yet). The loaded
        % panel is also recorded in data.cfgData for a fast reopen.
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
                    addCurationSet(inp, p.saveFcn);
                    gui_notify(hFig, sprintf('Loaded %d events from "%s".', nEv, p.label), 'success');
                case 'stateStrip'
                    nLab = numel(inp.data.labels);
                    addCurationSet(inp, p.saveFcn);
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
        % add a signal input (unique name) and open a panel for it in the chosen
        % region. The curated target and every other loaded set are untouched.
        data = hFig.UserData;
        inp.name = uniqueName(inp.name, data.inputs);
        inp.defRegion = reg;
        data.inputs = addInput(data.inputs, inp);
        data.Tend_s = max(eps, computeTend(data.inputs));
        f = regField(reg);
        data.(f) = [data.(f), makePanel(inp.name, reg, data.inputs)];
        hFig.UserData = data;
        rebuildPlot();
        populateSrc(hFig);
    end

    function addCurationSet(inp, sv)
        % add an event / state set (unique name) as a Top strip and list it in
        % CURATE. It is visible at once; it becomes the editable target only if
        % nothing is curated yet (otherwise it stays context and you elect it
        % from CURATE). The strip carries the set's own save handle.
        data = hFig.UserData;
        inp.name = uniqueName(inp.name, data.inputs);
        inp.defRegion = 'wide';
        inp.saveFcn = sv;
        data.inputs = addInput(data.inputs, inp);
        data.inputs = assignEvtColors(data.inputs);
        data.Tend_s = max(eps, computeTend(data.inputs));
        if ~any(strcmp({data.wideP.source}, inp.name)) && ...
                ~any(strcmp({data.narrowP.source}, inp.name))
            data.wideP = [data.wideP, makePanel(inp.name, 'wide', data.inputs)];
        end
        autoCur = isempty(data.curate);
        prevMode = data.mode;
        if autoCur, data.curate = inp.name; data = syncCurate(data); end
        if autoCur && ~strcmp(prevMode, data.mode)
            data = buildActionWidget(data);          % first target sets the widget
        end
        hFig.UserData = data;
        rebuildPlot();
        populateSrc(hFig);
        refreshCurateDD(hFig);
        updateTitle(hFig);
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
        % (re)build the tiledlayout from dd.wideP / dd.narrowP. Per-panel heights
        % become integer row spans; a 1-row black tile divides the regions; a
        % blank gap above the divider holds the Top region's x-axis, which prints
        % on its bottom panel. A larger K makes the divider proportionally
        % thinner while preserving the panel ratios.
        if isfield(dd, 'tl') && ~isempty(dd.tl) && isvalid(dd.tl), delete(dd.tl); end
        K = 20;
        nW = numel(dd.wideP); nN = numel(dd.narrowP);
        wSpans = ones(1, nW); for i = 1:nW, wSpans(i) = max(1, round(dd.wideP(i).height * K)); end
        nSpans = ones(1, nN); for i = 1:nN, nSpans(i) = max(1, round(dd.narrowP(i).height * K)); end
        hasDiv = nW > 0 && nN > 0;
        base = sum(wSpans) + sum(nSpans) + hasDiv;
        xgap = 0;                          % blank rows for the Top x-tick labels
        if hasDiv, xgap = max(3, round(0.04 * base)); end
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
            title(ax, ''); styleYLabel(ax, pn.label, pn.source, data.inputs);
            ax.XLim = [a, b] / xf;
            [data.hInd(i), data.hCenter(i)] = addCursor(ax);
            mute(ax);
        end
        fig.UserData = data;
        hideOuterXTicks(data.axWide, ['Time (' data.unit.wide ')']);
        updateMarker(fig);
    end

    function renderCurationStrips(fig)
        % redraw every Top curation strip (eventTicks / stateStrip) - used on a
        % curate switch, where the outgoing and incoming sets both change look.
        data = fig.UserData;
        for i = 1:numel(data.wideP)
            inp = getInput(data.inputs, data.wideP(i).source);
            if ~isempty(inp) && any(strcmp(inp.type, {'eventTicks', 'stateStrip'}))
                renderWideSource(fig, data.wideP(i).source);
            end
        end
    end

    function renderNarrow(fig)
        % Bottom window [t0 +/- win/2], x in the Bottom's unit; redrawn as t0 or
        % win changes (the Top is untouched, only its markers move). Every event
        % set draws its in-window marks (its own colour); the curated set's
        % current event is emphasised (peak solid, start / stop dashed).
        data = fig.UserData;
        if isempty(data.narrowP), return; end
        xf = unitSec(data.unit.narrow);
        ws = data.t0 - data.win / 2;
        we = data.t0 + data.win / 2;
        if ws < 0,            we = we - ws;              ws = 0; end
        if we > data.Tend_s,  ws = ws - (we - data.Tend_s); we = data.Tend_s; end
        ws = max(0, ws);
        for i = 1:numel(data.narrowP)
            pn = data.narrowP(i); ax = pn.ax;
            cla(ax); hold(ax, 'on');
            drawPanel(data, ax, pn.source, ws, we, xf);
            drawWindowMarks(ax, data, ws, we, xf);
            ax.XLim = [ws, we] / xf;
            styleYLabel(ax, pn.label, pn.source, data.inputs);
            mute(ax);
        end
        hideOuterXTicks(data.axNarrow, ['Time (' data.unit.narrow ')']);
        fig.UserData = data;
    end

    function drawPanel(data, ax, source, a, b, xf)
        % resolve the panel's input and hand it to guiPath_draw; x in display
        % units (= seconds / xf), covering the seconds range [a, b]. Each event
        % set draws from its OWN input, so a tick strip shows regardless of which
        % set (if any) is the curated target - drawTicks no-ops on empty data.
        inp = getInput(data.inputs, source);
        if isempty(inp), return; end
        guiPath_draw(ax, inp, data, a, b, xf);
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
        % accept-reject is read off the plot panels, so it carries no status
        % text; the status line shows the event's vigilance state instead.
        data = fig.UserData;
        if isempty(data.ev), return; end
        if strcmp(data.mode, 'view')
            data.ev.refresh(data.t0, data.Tend_s);
            return;
        end
        if strcmp(data.mode, 'states')
            if data.nEvents == 0, data.ev.refresh(0, 0, []); return; end
            data.ev.refresh(data.currIdx, data.nEvents, data.labels(data.currIdx));
            return;
        end
        if ~data.hasEvents || data.nEvents == 0
            data.ev.refresh(0, 0, false, '');
            return;
        end
        data.ev.refresh(data.currIdx, data.nEvents, data.accepted(data.currIdx), ...
            evtStatus(data));
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

    function navWindow(dirn)
        % view mode (CURATE = None): step the window by ~its own width, keeping a
        % 20% overlap so nothing at the seam is missed. dirn = +1 next, -1 prev.
        data = hFig.UserData;
        step = data.win * 0.8;                     % 20% overlap between windows
        data.t0 = min(max(0, data.t0 + dirn * step), data.Tend_s);
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
        renderWideSource(hFig, data.curate);   % recolour the curated Top strip
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
        renderWideSource(hFig, data.curate);   % recolour the curated Top strip
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
        setModShift(any(strcmpi(evt.Modifier, 'shift')));   % track for shift+scroll
        if isEditingField(), return; end   % typing in a field: let it keep the key
        if any(strcmpi(evt.Modifier, 'control')) && strcmpi(evt.Key, 's')
            onSave(); return;
        end
        data = hFig.UserData;
        % shift + / - / 0 : amplitude of the active (last-clicked) panel. Before
        % the mode switch, so it works while curating events or scoring states.
        if any(strcmpi(evt.Modifier, 'shift'))
            switch evt.Key
                case {'equal', 'add'},       adjustAmp(data.activeSrc, 'in');    return;
                case {'hyphen', 'subtract'}, adjustAmp(data.activeSrc, 'out');   return;
                case {'0', 'numpad0'},       adjustAmp(data.activeSrc, 'reset'); return;
            end
        end
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
        if strcmp(data.mode, 'view')
            switch evt.Key
                case {'leftarrow', 'downarrow'},  navWindow(-1);
                case {'rightarrow', 'uparrow'},   navWindow(1);
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
        up = evt.VerticalScrollCount < 0;
        % Shift+scroll adjusts the amplitude of the panel under the pointer, in
        % either region (up = bigger signal / brighter spec); plain scroll zooms
        % the pointed region's time axis. Over a raster / tick strip it no-ops.
        if data.modShift
            src = srcAtPointer(data, cp);
            if ~isempty(src), adjustAmp(src, ternary(up, 'in', 'out')); end
            return;
        end
        if evt.VerticalScrollCount > 0, f = 1.4; else, f = 1 / 1.4; end
        if ~isempty(data.axNarrow) && pointerOver(data.axNarrow, cp)
            zoomWin(f);
        elseif ~isempty(data.axWide) && pointerOver(data.axWide, cp)
            zoomTop(f, data.axWide(1).CurrentPoint(1, 1));
        end
    end

    function onKeyRelease(~, evt)
        % clear the Shift flag the moment Shift comes up, so a later plain scroll
        % is not mistaken for shift+scroll
        if strcmpi(evt.Key, 'shift')
            setModShift(false);
        else
            setModShift(any(strcmpi(evt.Modifier, 'shift')));
        end
    end

    function setModShift(tf)
        data = hFig.UserData;
        if data.modShift ~= tf, data.modShift = tf; hFig.UserData = data; end
    end

    function adjustAmp(src, action)
        % nudge one panel's yAdjust factor (in = bigger signal / brighter spec,
        % reset = back to the loaded look) and redraw it. Works in either region.
        % Limited to trace / traces / spec - a raster or a tick strip has no
        % amplitude to scale, so pointing at one does nothing.
        data = hFig.UserData;
        if isempty(src) || ~knownSrc(data, src)
            src = firstAmpSrc(data);   % nothing clicked yet -> first amp panel
        end
        ix = find(strcmp({data.inputs.name}, src), 1);
        if isempty(ix) || ~any(strcmp(data.inputs(ix).type, {'trace', 'traces', 'spec'}))
            return;
        end
        ya = data.inputs(ix).yAdjust;
        if isempty(ya) || ya <= 0, ya = 1; end
        switch action
            case 'in',    ya = min(ya * 1.25, 50);
            case 'out',   ya = max(ya / 1.25, 0.05);
            case 'reset', ya = 1;
        end
        data.inputs(ix).yAdjust = ya;
        hFig.UserData = data;
        if ~isempty(data.narrowP) && any(strcmp({data.narrowP.source}, src))
            renderNarrow(hFig);
        end
        if ~isempty(data.wideP) && any(strcmp({data.wideP.source}, src))
            renderWideSource(hFig, src);
        end
    end

    function renderWideSource(fig, src)
        % redraw one Top panel in place, by source name: clear, draw over the
        % whole session, restore its cursor + window band. Used to recolour the
        % curated strip after an accept / assign and on a curate switch.
        data = fig.UserData;
        xf = unitSec(data.unit.wide);
        for i = 1:numel(data.wideP)
            if strcmp(data.wideP(i).source, src)
                ax = data.wideP(i).ax;
                cla(ax); hold(ax, 'on');
                drawPanel(data, ax, src, 0, data.Tend_s, xf);
                styleYLabel(ax, data.wideP(i).label, data.wideP(i).source, data.inputs);
                ax.XLim = [0, data.Tend_s] / xf;
                [data.hInd(i), data.hCenter(i)] = addCursor(ax);
                mute(ax);
            end
        end
        fig.UserData = data;
        hideOuterXTicks(data.axWide, ['Time (' data.unit.wide ')']);
        updateMarker(fig);
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
        % click any panel -> that region becomes active (so +/-/0 and the Region
        % dropdown follow the click), the panel becomes the amplitude target for
        % shift+/-/0, and t0 moves to the clicked time (driving the Bottom).
        setActiveRegion(region);
        data = hFig.UserData;                          % re-read: setActiveRegion wrote it
        data.activeSrc = panelSrcFromAx(data, ax);
        hFig.UserData = data;
        jumpToTime(ax.CurrentPoint(1, 1) * unitSec(data.unit.(region)));
    end

    function finalizeInteractions(fig)
        data = fig.UserData;
        for k = 1:numel(data.axWide),   armAxis(data.axWide(k),   'wide');   end
        for k = 1:numel(data.axNarrow), armAxis(data.axNarrow(k), 'narrow'); end
        data.jumpToTimeFcn = @jumpToTime;       % exposed for hosts / tests
        data.loadCoreFcn   = @loadCore;         % exposed for tests (drives the Load flow)
        data.loadPresetFcn = @loadPreset;       % exposed for tests (drives a preset switch)
        data.setCurateFcn  = @setCurate;        % exposed for tests (drives a curate switch)
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

    function setActiveRegion(region)
        % make one region active (set by clicking a panel): update the Panels
        % header + the controls that read the region (the panel list, the
        % x-unit). Guarded so a click within the already-active region does
        % nothing (the rebuild is not free).
        data = hFig.UserData;
        if strcmp(data.cfgRegion, region), return; end
        data.cfgRegion = region;
        hFig.UserData = data;
        data.hPanelsLbl.Text = ['Panels (' regDisp(region) ')'];
        populateSrc(hFig);
        data.hUnit.Value = data.unit.(region);
    end

    % --- panel-list edits: each mutates the active region's panel array, then
    % rebuilds the plot and the list. addPanel / removePanel / movePanel are the
    % list's Add / X / up-down; onSrcChange is a row's dropdown.
    function addPanel()
        data = hFig.UserData; region = data.cfgRegion; f = regField(region);
        opts = availSources(data.inputs);
        if isempty(opts)
            gui_notify(hFig, 'Load a signal or events first (no sources to show).', 'warning');
            return;
        end
        data.(f)(end + 1) = makePanel(opts{1}, region, data.inputs);
        hFig.UserData = data;
        rebuildPlot(); populateSrc(hFig);
    end

    function removePanel(k)
        data = hFig.UserData; f = regField(data.cfgRegion); pArr = data.(f);
        if k < 1 || k > numel(pArr) || numel(pArr) <= 1, return; end   % keep >= 1
        pArr(k) = []; data.(f) = pArr;
        hFig.UserData = data;
        rebuildPlot(); populateSrc(hFig);
    end

    function movePanel(k, step)
        data = hFig.UserData; f = regField(data.cfgRegion); pArr = data.(f);
        j = k + step;
        if k < 1 || k > numel(pArr) || j < 1 || j > numel(pArr), return; end
        pArr([k, j]) = pArr([j, k]); data.(f) = pArr;
        hFig.UserData = data;
        rebuildPlot(); populateSrc(hFig);
    end

    function onSrcChange(k, src)
        data = hFig.UserData;
        f = regField(data.cfgRegion);
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
        % rebuild the active region's panel list: a row per panel (source
        % dropdown + up / down / X), then an Add row.
        data = fig.UserData;
        pArr = data.(regField(data.cfgRegion));
        g = data.hSrcGrid;
        delete(g.Children);
        nP = numel(pArr);
        g.RowHeight = repmat({'fit'}, 1, nP + 1);
        g.ColumnWidth = {'1x', 22, 22, 22};
        opts = availSources(data.inputs);
        for k = 1:nP
            val = pArr(k).source;
            if ~any(strcmp(val, opts)), val = opts{1}; end
            dd = uidropdown(g, 'Items', opts, 'Value', val, ...
                'ValueChangedFcn', @(s, ~) onSrcChange(k, s));
            dd.Layout.Row = k; dd.Layout.Column = 1;
            bU = uibutton(g, 'Text', char(9650), 'ButtonPushedFcn', @(~, ~) movePanel(k, -1));
            bU.Layout.Row = k; bU.Layout.Column = 2; bU.Enable = tf2e(k > 1);
            bD = uibutton(g, 'Text', char(9660), 'ButtonPushedFcn', @(~, ~) movePanel(k, +1));
            bD.Layout.Row = k; bD.Layout.Column = 3; bD.Enable = tf2e(k < nP);
            bX = uibutton(g, 'Text', 'X', 'ButtonPushedFcn', @(~, ~) removePanel(k));
            bX.Layout.Row = k; bX.Layout.Column = 4; bX.Enable = tf2e(nP > 1);
        end
        bA = uibutton(g, 'Text', 'Add', 'ButtonPushedFcn', @(~, ~) addPanel());
        bA.Layout.Row = nP + 1; bA.Layout.Column = [1, 4];
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
            data = writeBackCurate(data);
            hFig.UserData = data;
            gui_notify(hFig, sprintf('Saved %d epoch labels.', data.nEvents), 'success');
            return;
        end
        if ~data.hasEvents, return; end
        data.saveFcn(data.accepted(:));
        data.dirty = false;
        data = writeBackCurate(data);
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
% the caller opens an empty Custom view; guiPath never requires a file)
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

function items = curateItems(data)
% the CURATE selector's items: None + every loaded event / state set (by name)
names = {};
for k = 1:numel(data.inputs)
    if any(strcmp(data.inputs(k).type, {'eventTicks', 'stateStrip'}))
        names{end + 1} = data.inputs(k).name; %#ok<AGROW>
    end
end
items = [{'None'}, names];
end

function v = curateDisp(data)
% the CURATE selector's value for the current target ('None' when nothing is)
if ~isfield(data, 'curate') || isempty(data.curate), v = 'None'; else, v = data.curate; end
end

function refreshCurateDD(fig)
% resync the CURATE selector's items + value with the loaded sets
data = fig.UserData;
items = curateItems(data);
v = curateDisp(data);
if ~any(strcmp(v, items)), v = 'None'; end
data.hCurateDD.Items = items;
data.hCurateDD.Value = v;
end

function updateTitle(fig)
% window title: session - preset [curating: SET | view]
data = fig.UserData;
if isfield(data, 'curate') && ~isempty(data.curate)
    tag = ['curating: ' data.curate];
else
    tag = 'view';
end
fig.Name = sprintf('%s - %s [%s]', data.basename, data.presetName, tag);
end

function d = applyConfig(d, config)
% apply a resolved config onto state struct d. All inputs - signals AND every
% event / state set - live in one flat list d.inputs, which ACCUMULATES across
% preset switches: a preset changes the arrangement (which panels are shown) and
% loads its own data, but never drops what is already loaded, so every signal /
% set stays available in the panel list and in CURATE. The editable target
% (d.curate) is chosen only on the first open; a later preset switch leaves it
% alone - PRESET (arrangement) and CURATE (target) are independent.
newInputs = normalizeInputs(config.inputs);
% pin the preset's resolved save handle onto its own target set
if isfield(config, 'saveFcn') && ~isempty(config.saveFcn)
    tix = find(ismember({newInputs.type}, {'eventTicks', 'stateStrip'}), 1);
    if ~isempty(tix), newInputs(tix).saveFcn = config.saveFcn; end
end
presetTarget = '';
tix = find(ismember({newInputs.type}, {'eventTicks', 'stateStrip'}), 1);
if ~isempty(tix), presetTarget = newInputs(tix).name; end

firstOpen = ~isfield(d, 'curate');
if ~isfield(d, 'inputs') || isempty(d.inputs)
    d.inputs = newInputs;
else
    d.inputs = mergeInputs(d.inputs, newInputs);   % accumulate; keep what is loaded
end
d.inputs = assignEvtColors(d.inputs);
if isfield(config, 'winPlot') && ~isempty(config.winPlot)
    d.win = config.winPlot; d.winDefault = config.winPlot;
end
d.Tend_s = max(eps, computeTend(d.inputs));

% target auto-selects the preset's set ONLY on the first open; a preset switch
% keeps whatever is being curated (independent of the arrangement)
if firstOpen, d.curate = presetTarget; end
if ~isfield(d, 'currIdx') || isempty(d.currIdx), d.currIdx = 1; end
d = syncCurate(d);
if d.hasEvents && d.nEvents > 0, d.t0 = d.ed.peakTime(1); else, d.t0 = d.Tend_s / 2; end

% the arrangement (which panels are shown, and where) follows the preset
panels = normalizePanels(config, d.inputs);
d.wideP   = panels(strcmp({panels.region}, 'wide'));
d.narrowP = panels(strcmp({panels.region}, 'narrow'));
if isempty(d.wideP) && ~isempty(d.narrowP), d.cfgRegion = 'narrow'; else, d.cfgRegion = 'wide'; end
end

function d = syncCurate(d)
% mirror the curated input (named d.curate; '' = None) into the live working
% state: the event list + accept mask (events) or the labels + epoch centres
% (states), the mode, and the save handle. Navigation reads this working copy;
% edits are written back into the input by writeBackCurate on a switch / save.
ix = curateIx(d);
if isempty(ix)
    d.mode = 'view'; d.saveFcn = [];   % None: browse the window, no editable target
    d.ed = struct('peakTime', []); d.accepted = logical([]);
    d.nEvents = 0; d.hasEvents = false;
    d.labels = []; d.epochLen = 1; d.nstates = 0; d.stateNames = {}; d.stateColors = {};
    d.currIdx = 1;
    return;
end
inp = d.inputs(ix);
d.saveFcn = inp.saveFcn;
if strcmp(inp.type, 'stateStrip')
    d.mode = 'states';
    [d.labels, epochT, d.epochLen, d.nstates, d.stateNames, d.stateColors] = ...
        stateFromData(inp.data);
    nEp = numel(d.labels);
    d.ed = struct('peakTime', epochT(:));
    d.accepted = true(max(0, nEp), 1);      % unused in states mode
    d.nEvents = nEp; d.hasEvents = nEp > 0;
else
    d.mode = 'events';
    [d.ed, d.accepted, d.nEvents, d.hasEvents] = eventFromData(inp.data);
    d.labels = []; d.epochLen = 1; d.nstates = 0; d.stateNames = {}; d.stateColors = {};
end
if ~isfield(d, 'currIdx') || isempty(d.currIdx) || d.currIdx < 1, d.currIdx = 1; end
d.currIdx = min(d.currIdx, max(1, d.nEvents));
end

function ix = curateIx(d)
% index of the curated input in d.inputs; [] for None / not found
ix = [];
if ~isfield(d, 'curate') || isempty(d.curate) || isempty(d.inputs), return; end
ix = find(strcmp({d.inputs.name}, d.curate), 1);
end

function d = writeBackCurate(d)
% copy the live edits back into the curated input's data, so a set keeps its
% edits when CURATE moves away and a de-selected set draws its edited ticks
ix = curateIx(d);
if isempty(ix) || ~isstruct(d.inputs(ix).data), return; end
if strcmp(d.inputs(ix).type, 'stateStrip')
    d.inputs(ix).data.labels = d.labels(:);
elseif d.nEvents > 0
    d.inputs(ix).data.accepted = logical(d.accepted(:));
end
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
    'defRegion', regAlias(pc.region), 'defOrder', order, ...
    'chInfo', {pick(pc, 'chInfo', [])}, 'yAdjust', pick(pc, 'yAdjust', 1));
end

function r = regAlias(region)
% user-facing top/bottom -> internal wide/narrow (both accepted)
switch lower(region)
    case {'top', 'wide'},      r = 'wide';
    case {'bottom', 'narrow'}, r = 'narrow';
    otherwise, r = region;
end
end

function [ed, accepted, nEv, has] = eventFromData(D)
% the live event state from a set's data struct (peakTime + optional times /
% accepted / state)
ed = []; accepted = logical([]); nEv = 0; has = false;
if ~isstruct(D) || ~isfield(D, 'peakTime') || isempty(D.peakTime), return; end
ed = D; ed.peakTime = D.peakTime(:);
nEv = numel(ed.peakTime); has = true;
if isfield(D, 'accepted') && numel(D.accepted) == nEv
    accepted = logical(D.accepted(:));
else
    % A stored mask whose length has drifted from the event list cannot be
    % mapped back (it is positional), so it is dropped. Say so: silently
    % starting all-accepted looks identical to a fresh detection, and the next
    % save would persist that reset over real curation.
    if isfield(D, 'accepted') && ~isempty(D.accepted)
        warning('guiPath:acceptedLen', ...
            ['stored .accepted has %d entries but the list holds %d events; ' ...
            'ignoring it and starting all-accepted. Saving WILL overwrite ' ...
            'the stored mask.'], numel(D.accepted), nEv);
    end
    accepted = true(nEv, 1);
end
end

function s = evtStatus(data)
% status line for the event stepper: the current event's vigilance state. The
% event list carries .state when it came from ripp / ed (evt_states writes it);
% a bare times vector does not, and then the line stays empty.
s = '';
if ~isstruct(data.ed) || ~isfield(data.ed, 'state') || isempty(data.ed.state)
    return;
end
if data.currIdx < 1 || data.currIdx > numel(data.ed.state), return; end
v = data.ed.state(data.currIdx);
if ismissing(v), s = 'State: undefined'; else, s = ['State: ', char(string(v))]; end
end

function [labels, epochT, epochLen, nstates, names, colors] = stateFromData(D)
% the live state-scoring state from a stateStrip set's data struct
labels = []; epochT = []; epochLen = 1; nstates = 0; names = {}; colors = {};
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

function drawEventMarks(ax, m, clr)
% the current event on a Bottom panel: the peak solid, start / stop dashed, all
% in the event colour. m is [start peak stop] (from eventMarks) or just a peak.
if numel(m) == 3
    xline(ax, m([1, 3]), '--', 'Color', clr);
    xline(ax, m(2), '-', 'Color', clr);
else
    xline(ax, m, '-', 'Color', clr);
end
end

function drawWindowMarks(ax, data, ws, we, xf)
% overlay, on a Bottom signal panel, the in-window marks of every event set in
% its own colour (one NaN-separated line per set, so the object count is set-
% bounded, not event-bounded, and the marks filter to [ws, we] first). The
% curated set's current event is then emphasised on top (peak solid, start /
% stop dashed). State sets never mark (they stay the coloured strip).
yl = ax.YLim;
for k = 1:numel(data.inputs)
    inp = data.inputs(k);
    if ~strcmp(inp.type, 'eventTicks'), continue; end
    if ~isstruct(inp.data) || ~isfield(inp.data, 'peakTime'), continue; end
    pk = inp.data.peakTime(:);
    if strcmp(inp.name, data.curate)
        if numel(data.accepted) == numel(pk), pk = pk(data.accepted); end
    elseif isfield(inp.data, 'accepted') && numel(inp.data.accepted) == numel(pk)
        pk = pk(logical(inp.data.accepted));
    end
    inWin = pk(pk >= ws & pk <= we);
    drawSpanLines(ax, inWin / xf, yl, inp.clr);
end
% emphasise the curated set's current event
if data.hasEvents && data.nEvents > 0 && strcmp(data.mode, 'events')
    ev = data.currIdx;
    if data.ed.peakTime(ev) >= ws && data.ed.peakTime(ev) <= we
        cix = curateIx(data);
        clr = data.clrEvt;
        if ~isempty(cix), clr = data.inputs(cix).clr; end
        drawEventMarks(ax, eventMarks(data.ed, ev) / xf, clr);
    end
end
end

function drawSpanLines(ax, x, yl, clr)
% one full-height vertical line per x (across the panel's current y-limits),
% drawn as a single NaN-separated line object
if isempty(x), return; end
x = x(:)';
X = [x; x; nan(1, numel(x))];
Y = repmat([yl(1); yl(2); NaN], 1, numel(x));
line(ax, X(:), Y(:), 'Color', clr);
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

function inputs = addInput(inputs, inp)
% append INP if its name is not already present (the name is assumed unique-ised)
if isempty(inputs) || ~any(strcmp(inp.name, {inputs.name}))
    inputs(end + 1) = inp;
end
end

function inputs = mergeInputs(inputs, incoming)
% accumulate: append each incoming input whose name is new, and KEEP an existing
% one as-is (it may hold live curation edits, and a preset's data for a shared
% name is the same file). This is what lets loaded signals / sets persist across
% a preset switch, so PRESET only changes the arrangement, not what is available.
for i = 1:numel(incoming)
    inputs = addInput(inputs, incoming(i));
end
end

function nm = uniqueName(base, inputs)
% BASE, or BASE_k for the first k that no existing input already uses
nm = base; k = 1;
while ~isempty(inputs) && any(strcmp(nm, {inputs.name}))
    k = k + 1; nm = sprintf('%s_%d', base, k);
end
end

function inputs = assignEvtColors(inputs)
% give each event set a stable, distinct colour (used for its ticks + its
% window marks). The first is the classic blue, so a single-set session is
% unchanged.
if isempty(inputs), return; end
ix = find(strcmp({inputs.type}, 'eventTicks'));
for j = 1:numel(ix)
    inputs(ix(j)).clr = evtPalette(j);
end
end

function clr = evtPalette(k)
% distinct event-set colours, cycled past the end
P = [0.00 0.20 0.80;      % blue   (matches d.clrEvt / the single-set look)
     0.85 0.33 0.10;      % orange
     0.20 0.60 0.20;      % green
     0.55 0.20 0.60;      % purple
     0.80 0.60 0.10;      % gold
     0.10 0.60 0.65];     % teal
k = mod(k - 1, size(P, 1)) + 1;
clr = P(k, :);
end

%% ========================================================================
%  INPUTS MODEL (pure)
%  ========================================================================

function b = blankInput()
b = struct('name', '', 'type', '', 'data', [], 'fs', NaN, 'ylim', [], ...
    'clr', 'k', 'label', '', 'height', 1, 'defRegion', 'narrow', 'defOrder', 99, ...
    'chInfo', [], ...      % traces only: stack display stats (see traceStack)
    'yAdjust', 1, ...      % live amplitude factor (shift+scroll / shift+-0)
    'saveFcn', []);        % events / states: this set's save handle
end

function out = normalizeInputs(in)
% fill defaults on an inputs struct array (requires .name .type). Per-type
% appearance is set upstream by guiPath_panel, so these are generic fallbacks only.
n = numel(in);
out = repmat(blankInput(), 1, n);
for i = 1:n
    s = in(i);
    if ~isfield(s, 'name') || isempty(s.name), error('guiPath:inputs', 'input %d needs a name', i); end
    if ~isfield(s, 'type') || isempty(s.type), error('guiPath:inputs', 'input %d needs a type', i); end
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
    out(i).chInfo    = pick(s, 'chInfo', []);
    out(i).yAdjust   = pick(s, 'yAdjust', 1);
    out(i).saveFcn   = pick(s, 'saveFcn', []);
end
end

function v = pick(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end

function out = ternary(tf, a, b)
if tf, out = a; else, out = b; end
end

function e = tf2e(tf)
% logical -> an Enable string for a uicomponent
if tf, e = 'on'; else, e = 'off'; end
end

function styleYLabel(ax, txt, source, inputs)
% a panel's y-label: muted and compact. A thin marker strip (eventTicks /
% stateStrip / hypnogram) gets a HORIZONTAL label so its short name sits within
% the strip's height instead of a vertical label overflowing into the neighbour
% panels; a tall signal panel keeps the usual vertical label.
inp = getInput(inputs, source);
if ~isempty(inp) && isStripType(inp.type)
    ylabel(ax, txt, 'Color', [0.15 0.15 0.15], 'FontSize', 9, 'Rotation', 0, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
else
    ylabel(ax, txt, 'Color', [0.15 0.15 0.15], 'FontSize', 9);
end
end

function src = panelSrcFromAx(data, ax)
% the source name of the panel drawn into axis ax (either region), '' if none
src = '';
for fld = {'wideP', 'narrowP'}
    P = data.(fld{1});
    for i = 1:numel(P)
        if isvalid(P(i).ax) && P(i).ax == ax, src = P(i).source; return; end
    end
end
end

function tf = knownSrc(data, src)
% is src the source of any panel, in either region?
tf = (~isempty(data.narrowP) && any(strcmp({data.narrowP.source}, src))) || ...
     (~isempty(data.wideP)   && any(strcmp({data.wideP.source},   src)));
end

function src = firstAmpSrc(data)
% source of the first trace / traces panel (Bottom preferred), '' if none. The
% keyboard default when no panel has been clicked yet; spec is not offered here
% - brightness is a pointed / clicked action, not the no-target fallback.
src = '';
for fld = {'narrowP', 'wideP'}
    P = data.(fld{1});
    for i = 1:numel(P)
        ix = find(strcmp({data.inputs.name}, P(i).source), 1);
        if ~isempty(ix) && any(strcmp(data.inputs(ix).type, {'trace', 'traces'}))
            src = P(i).source; return;
        end
    end
end
end

function src = srcAtPointer(data, cp)
% the source of the panel (either region) whose pixel rect holds the pointer cp
src = '';
for fld = {'narrowP', 'wideP'}
    P = data.(fld{1});
    for i = 1:numel(P)
        ax = P(i).ax;
        if ~isvalid(ax), continue; end
        r = getpixelposition(ax, true);
        if cp(1) >= r(1) && cp(1) <= r(1) + r(3) && ...
                cp(2) >= r(2) && cp(2) <= r(2) + r(4)
            src = P(i).source; return;
        end
    end
end
end

function T = computeTend(inputs)
% session length [s] = the largest time spanned by any input
T = eps;
for i = 1:numel(inputs)
    inp = inputs(i);
    switch inp.type
        case {'trace', 'traces'}
            % size(,1), not numel: a traces panel holds [nSamples x nCh], and
            % numel would scale the session length by the channel count
            if ~isempty(inp.data) && isfinite(inp.fs) && inp.fs > 0
                T = max(T, (size(inp.data, 1) - 1) / inp.fs);
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

function hideOuterXTicks(axs, label)
% the region's Time axis (ticks + label) prints on its BOTTOM panel, whatever it
% is (a marker strip too); every panel above shows its tick marks aligned but no
% numbers or label. buildPlot leaves a blank gap below the bottom panel for the
% labels.
n = numel(axs);
for i = 1:n
    if ~isvalid(axs(i)), continue; end
    if i == n
        axs(i).XAxisLocation = 'bottom';
        xlabel(axs(i), label);
    else
        axs(i).XTickLabel = [];
        xlabel(axs(i), '');
    end
end
end

function tf = isStripType(type)
% a thin marker / label strip: gets a horizontal y-label (styleYLabel) and no
% y-ticks, but still carries the x-axis when it is a region's bottom panel
tf = any(strcmp(type, {'eventTicks', 'stateStrip', 'hypnogram'}));
end

function mute(ax)
% make plotted children non-pickable so the axis ButtonDownFcn fires on click
ch = allchild(ax);
if ~isempty(ch), set(ch, 'HitTest', 'off'); end
end

% EOF
