function hFig = ed_gui(basepath, ed, sSig, varargin)
% ED_GUI Config-driven, general-purpose curation viewer for events (EDs).
%
%   hFig = ED_GUI(basepath, ed, sSig, varargin)
%
%   SUMMARY:
%       Single-session raw-signal viewer on uifigure, built on the +tblgui
%       helper layer. The plot area is described by a PANELS list, each entry
%       naming a data SOURCE and a REGION. Two stacks of panels are shown
%       simultaneously, separated by a black divider:
%           Top    - full-session overview (drawn once; x-range zoomable).
%           Bottom - a window around the cursor t0 (redrawn as t0 moves).
%       (Region keys are 'wide'/'narrow' internally; the UI labels them
%       Top/Bottom.) Any source can go in either region, and each region has
%       its own x-axis unit (ms|s|min|hr).
%
%       Every panel lives in one tiledlayout, so all plot boxes share one left
%       gutter and one width. The Top region prints its x-axis on top and the
%       Bottom on the bottom, so the divider never hides it.
%
%       Navigation: a single cursor t0 drives the Bottom window. Clicking ANY
%       panel (Top or Bottom) moves t0 to the clicked time, so the Bottom
%       re-centres there; prev/next/accept/reject and the event index set t0 to
%       an event. The Top overview is fixed (it does not follow t0); a cursor
%       line and a shaded band on it mark t0 and the Bottom window.
%
%       The Region dropdown picks the ACTIVE region for panel configuration,
%       the X-units dropdown, and the zoom buttons (Top scales the overview
%       x-range, Bottom scales the window).
%
%       Interaction:
%           Left / Right          previous / next event
%           Up / Down             accept / reject (auto-advances)
%           Ctrl+S                save
%           + / -                 zoom the active region in / out
%           0                     reset the active region
%           click any panel       move t0 (the Bottom) to the clicked time
%           scroll over a region  zoom that region (by pointer position)
%       Accept/reject writes ed.accepted; Save writes <basename>.ed.mat.
%       Panels and panel counts are reconfigurable live from the left column.
%
%   INPUTS:
%       basepath    - (Char) Session directory (holds <basename>.ed.mat).
%       ed          - (Struct) ED struct (ed_wrapper output). Requires
%                              .peakTime, .times, .accepted, .info.fs.
%       sSig        - (Struct) sleep_sig struct. If empty, loaded via ed_sigLoad.
%       varargin    - Parameter/Value pairs:
%           'panels'      - (Struct array) Each: .source .region .height .label.
%                           Sources: spec|emgRms|hypnogram|eventTicks|lfp|emg|raster.
%                           Regions: wide|narrow. Default: a sensible ED layout.
%           'specAdapter' - (Struct) Spectrogram for plot_spec (.s/.freq/.tstamps).
%           'winPlot'     - (Num) Bottom-region window full width [s]. {1.0}
%           'basename'    - (Char) Override (defaults to the folder name).
%           'Visible'     - (Char) 'on' (default) | 'off' (headless).
%
%   OUTPUT:
%       hFig        - (uifigure) Handle to the viewer window.
%
%   DEPENDENCIES:
%       tblgui.layout / labeledControl / eventPanel / notify / chooseDialog,
%       plot_spec, plot_hypnogram (style 'strip'), plot_raster, ed_sigLoad,
%       basepaths2vars, as_loadConfig.
%
%   HISTORY:
%       Created:  22 Jun 2026
%       Updated:  23 Jun 2026 - single-tiledlayout plot area; global cursor t0
%                               (Top overview / Bottom window); per-region
%                               x-units; swappable event module; live config.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'basepath', @ischar);
addRequired(p, 'ed', @isstruct);
addRequired(p, 'sSig', @(x) isempty(x) || isstruct(x));
addParameter(p, 'panels', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'specAdapter', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'winPlot', 1.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'basename', '', @ischar);
addParameter(p, 'Visible', 'on', @(x) any(strcmpi(char(x), {'on', 'off'})));

parse(p, basepath, ed, sSig, varargin{:});
basepath    = p.Results.basepath;
ed          = p.Results.ed;
sSig        = p.Results.sSig;
panels      = p.Results.panels;
specAdapter = p.Results.specAdapter;
winPlot     = p.Results.winPlot;
vis         = char(p.Results.Visible);

basename = p.Results.basename;
if isempty(basename)
    [~, basename] = fileparts(basepath);
end

%% ========================================================================
%  LOAD LAYER (signals + events once)
%  ========================================================================

if isempty(sSig)
    [~, ~, ~, ~, specAdapter, sSig] = ed_sigLoad(basepath, 'basename', basename);
end
if isempty(specAdapter) && isfield(sSig, 'spec')
    specAdapter = struct('s', sSig.spec, 'freq', sSig.spec_freq(:), ...
        'tstamps', sSig.spec_tstamps(:));
end

d = struct();
d.basepath   = basepath;
d.basename   = basename;
d.edFile     = fullfile(basepath, [basename, '.ed.mat']);
d.fs         = ed.info.fs;
d.eeg        = sSig.eeg(:);
d.emg        = sSig.emg(:);
if isfield(sSig, 'emg_rms'), d.emgRms = sSig.emg_rms(:); else, d.emgRms = []; end
d.timestamps = (0:numel(d.eeg) - 1)' / d.fs;
d.Tend_s     = max(eps, d.timestamps(end));
d.Tend_hr    = d.Tend_s / 3600;
d.specAdapter = specAdapter;
d.ylimEeg    = prctile(d.eeg, [0.1, 99.9]);
d.ylimEmg    = prctile(d.emg, [0.1, 99.9]);

d.ed       = ed;
d.nEvents  = numel(ed.peakTime);
d.accepted = logical(ed.accepted(:));
d.clrAccept = [0.10 0.55 0.10];
d.clrReject = [0.65 0.15 0.15];

% Optional spikes (raster source appears only if present)
d.hasSpikes = false;
d.spktimes  = {};
try
    vsp = basepaths2vars('basepaths', {basepath}, 'vars', {'spikes'}, 'flgPrnt', false);
    if isfield(vsp, 'spikes') && isfield(vsp.spikes, 'times')
        d.spktimes  = vsp.spikes.times(:);
        d.hasSpikes = ~isempty(d.spktimes);
    end
catch
end

% Optional sleep states (for the hypnogram strip). ss.bouts.times is a cell
% indexed by state (linear index); build an hours cell of length cfg.nstates.
d.hasStates   = false;
d.boutTimesHr = {};
try
    vss = basepaths2vars('basepaths', {basepath}, 'vars', {'sleep_states'}, 'flgPrnt', false);
    if isfield(vss, 'ss') && isfield(vss.ss, 'bouts') && isfield(vss.ss.bouts, 'times')
        bt  = vss.ss.bouts.times;
        cfg = as_loadConfig();
        ns  = cfg.nstates;
        d.boutTimesHr = cell(1, ns);
        for s = 1:ns
            if s <= numel(bt) && ~isempty(bt{s})
                d.boutTimesHr{s} = bt{s} / 3600;
            else
                d.boutTimesHr{s} = zeros(0, 2);
            end
        end
        d.hasStates = true;
    end
catch
end

% View state: a single cursor t0 drives the Bottom window (width win). The Top
% is a fixed overview (its zoom lives in its axes' xlim). Each region keeps its
% own x-axis unit. cfgRegion is the ACTIVE region for config / units / zoom.
d.currIdx    = 1;
d.win        = winPlot;
d.winDefault = winPlot;
if d.nEvents > 0, d.t0 = d.ed.peakTime(1); else, d.t0 = d.Tend_s / 2; end
d.unit       = struct('wide', 'hr', 'narrow', 's');
d.hInd       = gobjects(0);
d.hCenter    = gobjects(0);
d.dirty      = false;
d.cfgRegion  = 'wide';

%% ========================================================================
%  PANELS
%  ========================================================================

if isempty(panels)
    panels = defaultPanels(d);
end
panels = filterAvail(panels, d);
panels = addAxField(panels);
wideP   = panels(strcmp({panels.region}, 'wide'));
narrowP = panels(strcmp({panels.region}, 'narrow'));

% the active region must have panels (a one-region config is legal)
if isempty(wideP) && ~isempty(narrowP), d.cfgRegion = 'narrow'; end

%% ========================================================================
%  FIGURE + LAYOUT
%  ========================================================================

hFig = uifigure('Name', [basename, ' - ED Curation'], ...
    'Position', [50, 50, 1320, 820], 'Visible', vis, ...
    'WindowKeyPressFcn', @onKey, 'WindowScrollWheelFcn', @onScroll);

[~, gPlot, gCtrl, gActions] = tblgui.layout(hFig, 'CtrlWidth', 240, 'CtrlSide', 'left');

% --- Controls: layout configuration (per region) ---
if strcmp(d.cfgRegion, 'wide'), nInitCfg = numel(wideP); else, nInitCfg = numel(narrowP); end
tblgui.labeledControl(gCtrl, 'label', 'LAYOUT', 'FontWeight', 'bold');
d.hRegionDD = tblgui.labeledControl(gCtrl, 'dropdown', 'Region', ...
    'Items', {'Top', 'Bottom'}, 'Value', regDisp(d.cfgRegion), 'ValueChangedFcn', @onRegionSel);
d.hPanelsN = tblgui.labeledControl(gCtrl, 'editnum', '# Panels', ...
    'Limits', [1, 6], 'RoundFractionalValues', 'on', 'Value', max(1, nInitCfg), ...
    'ValueChangedFcn', @onPanelsN);
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

% --- Event module (swappable) in the pinned action strip ---
api = struct('prev', @() navStep(-1), 'next', @() navStep(1), ...
    'accept', @() setAccept(true), 'reject', @() setAccept(false), ...
    'save', @() onSave(), 'setIdx', @(v) setIdx(v));
d.ev = tblgui.eventPanel(gActions, api);

% --- Plot area: one tiledlayout holding every panel ---
d.hPanel  = uipanel(gPlot, 'BorderType', 'none');
d.wideP   = wideP;
d.narrowP = narrowP;
d = buildPlot(d);

hFig.UserData = d;

%% ========================================================================
%  INITIAL RENDER
%  ========================================================================

populateSrc(hFig);
renderWideStatic(hFig);
renderNarrow(hFig);
finalizeInteractions(hFig);

drawnow;
if numel(d.axWide)   > 1, linkaxes(hFig.UserData.axWide, 'x');   end
if numel(d.axNarrow) > 1, linkaxes(hFig.UserData.axNarrow, 'x'); end

updateMarker(hFig);
refreshEvent(hFig);

hFig.CloseRequestFcn = @(~,~) onClose();

%% ========================================================================
%  PLOT CONSTRUCTION (single tiledlayout)
%  ========================================================================

    function dd = buildPlot(dd)
        % (re)build the tiledlayout from dd.wideP / dd.narrowP. Per-panel
        % heights become integer row spans; a 1-row black tile divides the
        % regions. Sets dd.tl, dd.axWide, dd.axNarrow, dd.divider, each .ax.
        if isfield(dd, 'tl') && ~isempty(dd.tl) && isvalid(dd.tl), delete(dd.tl); end
        K = 10;
        nW = numel(dd.wideP); nN = numel(dd.narrowP);
        wSpans = ones(1, nW); for i = 1:nW, wSpans(i) = max(1, round(dd.wideP(i).height * K)); end
        nSpans = ones(1, nN); for i = 1:nN, nSpans(i) = max(1, round(dd.narrowP(i).height * K)); end
        hasDiv = nW > 0 && nN > 0;
        totalRows = sum(wSpans) + hasDiv + sum(nSpans);

        tl = tiledlayout(dd.hPanel, max(1, totalRows), 1, ...
            'TileSpacing', 'none', 'Padding', 'tight');
        dd.tl = tl;

        r = 1;
        dd.axWide = gobjects(1, nW);
        for i = 1:nW
            ax = nexttile(tl, r, [wSpans(i), 1]); r = r + wSpans(i);
            dd.wideP(i).ax = ax; dd.axWide(i) = ax;
        end

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
            drawSource(data, ax, pn.source, a, b, xf);
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
        hideOuterXTicks(data.axWide, ['Time (' data.unit.wide ')'], 'top');
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
        hideOuterXTicks(data.axWide, ['Time (' data.unit.wide ')'], 'top');
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
        showEv = data.nEvents > 0 && data.ed.peakTime(ev) >= ws && data.ed.peakTime(ev) <= we;
        for i = 1:numel(data.narrowP)
            pn = data.narrowP(i); ax = pn.ax;
            cla(ax); hold(ax, 'on');
            drawSource(data, ax, pn.source, ws, we, xf);
            if showEv
                xline(ax, [data.ed.times(ev,1), data.ed.peakTime(ev), data.ed.times(ev,2)] / xf, ...
                    '--b', 'HandleVisibility', 'off');
            end
            ax.XLim = [ws, we] / xf;
            ylabel(ax, pn.label);
            set(get(ax, 'YLabel'), 'Color', [0.15 0.15 0.15]);
            mute(ax);
        end
        axTop = data.narrowP(1).ax;
        if showEv
            if data.accepted(ev)
                title(axTop, sprintf('Event %d/%d   [ACCEPTED]', ev, data.nEvents), 'Color', data.clrAccept);
            else
                title(axTop, sprintf('Event %d/%d   [REJECTED]', ev, data.nEvents), 'Color', data.clrReject);
            end
        else
            title(axTop, sprintf('t = %.2f s', data.t0), 'Color', [0.2 0.2 0.2]);
        end
        hideOuterXTicks(data.axNarrow, ['Time (' data.unit.narrow ')'], 'bottom');
        fig.UserData = data;
    end

    function drawSource(data, ax, source, a, b, xf)
        % draw one source into ax, x in display units (= seconds / xf),
        % covering the seconds range [a, b]
        switch source
            case 'spec'
                if ~isempty(data.specAdapter)
                    plot_spec(data.specAdapter, 'axh', ax, 'saveFig', false, 'xtime', xf);
                end
            case 'emgRms'
                if ~isempty(data.emgRms)
                    te = (0:numel(data.emgRms) - 1)';
                    sel = te >= a & te <= b;
                    plot(ax, te(sel) / xf, data.emgRms(sel), 'k');
                end
            case 'hypnogram'
                if data.hasStates
                    bt = cellfun(@(x) x * 3600 / xf, data.boutTimesHr, 'uni', false);
                    plot_hypnogram('boutTimes', bt, 'style', 'strip', 'hAx', ax);
                end
            case 'eventTicks'
                drawTicks(ax, data, xf);
            case 'lfp'
                drawSig(ax, data.eeg, data.fs, a, b, xf, data.ylimEeg);
            case 'emg'
                drawSig(ax, data.emg, data.fs, a, b, xf, data.ylimEmg);
            case 'raster'
                if data.hasSpikes
                    spk = cellfun(@(s) s(s >= a & s <= b) / xf, data.spktimes, 'uni', false);
                    plot_raster(spk, 'hAx', ax, 'xLim', [a, b] / xf, 'flgLbls', false);
                end
        end
    end

    function drawSig(ax, sig, fs, a, b, xf, yl)
        % windowed raw trace, decimated for display, x in display units
        s1 = max(1, floor(a * fs) + 1);
        s2 = min(numel(sig), ceil(b * fs) + 1);
        if s2 < s1, s2 = s1; end
        rng = s1:s2;
        t = ((rng - 1) / fs) / xf;
        np = numel(rng); maxPts = 20000;
        if np > maxPts
            st = ceil(np / maxPts);
            plot(ax, t(1:st:end), sig(rng(1:st:end)), 'k');
        else
            plot(ax, t, sig(rng), 'k');
        end
        ax.YLim = yl;
    end

    function drawTicks(ax, data, xf)
        pk = data.ed.peakTime(:) / xf;
        drawTickLine(ax, pk(data.accepted),  data.clrAccept);
        drawTickLine(ax, pk(~data.accepted), data.clrReject);
        ylim(ax, [0, 1]);
    end

    function drawTickLine(ax, x, clr)
        if isempty(x), return; end
        x = x(:)';
        X = [x; x; nan(1, numel(x))];
        Y = repmat([0; 1; NaN], 1, numel(x));
        line(ax, X(:), Y(:), 'Color', clr, 'HandleVisibility', 'off');
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
        data = fig.UserData;
        if data.nEvents == 0
            data.ev.refresh(0, 0, false, 'no events');
            return;
        end
        ev = data.currIdx;
        st = ''; if isfield(data.ed, 'state') && numel(data.ed.state) >= ev, st = char(data.ed.state(ev)); end
        amp = NaN; if isfield(data.ed, 'amp') && numel(data.ed.amp) >= ev, amp = data.ed.amp(ev); end
        emgZ = NaN; if isfield(data.ed, 'emgZ') && numel(data.ed.emgZ) >= ev, emgZ = data.ed.emgZ(ev); end
        statusText = sprintf('t=%.2fs  %s  amp=%.1f  emgZ=%.1f', data.ed.peakTime(ev), st, amp, emgZ);
        data.ev.refresh(ev, data.nEvents, data.accepted(ev), statusText);
    end

%% ========================================================================
%  NAVIGATION / STATE (t0 drives the Bottom; the Top is a fixed overview)
%  ========================================================================

    function navStep(step)
        data = hFig.UserData;
        if data.nEvents == 0, return; end
        data.currIdx = min(max(1, data.currIdx + step), data.nEvents);
        data.t0 = data.ed.peakTime(data.currIdx);
        hFig.UserData = data;
        renderNarrow(hFig); updateMarker(hFig); refreshEvent(hFig);
    end

    function setIdx(v)
        data = hFig.UserData;
        if data.nEvents == 0, return; end
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
        if data.nEvents == 0, return; end
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
            opts = availSources(region, data);
            for i = cur + 1:n
                pArr(i) = makePanel(opts{1}, region);
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
        pArr(k).source = src.Value;
        pArr(k).label  = srcLabel(src.Value);
        pArr(k).height = defaultHeight(src.Value);
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
        opts = availSources(region, data);
        for k = 1:numel(pArr)
            val = pArr(k).source;
            if ~any(strcmp(val, opts)), val = opts{1}; end
            dd = uidropdown(g, 'Items', opts, 'Value', val, ...
                'ValueChangedFcn', @(s, ~) onSrcChange(k, s));
            dd.Layout.Row = k; dd.Layout.Column = 1;
        end
    end

    function rebuildPlot()
        d = hFig.UserData;
        d = buildPlot(d);
        hFig.UserData = d;
        renderWideStatic(hFig);
        renderNarrow(hFig);
        finalizeInteractions(hFig);
        drawnow;
        d = hFig.UserData;
        if numel(d.axWide) > 1, linkaxes(d.axWide, 'x'); end
        if numel(d.axNarrow) > 1, linkaxes(d.axNarrow, 'x'); end
        updateMarker(hFig);
        refreshEvent(hFig);
    end

%% ========================================================================
%  SAVE / CLOSE
%  ========================================================================

    function onSave()
        data = hFig.UserData;
        ed = data.ed;                                       %#ok<PROPLC>
        ed.accepted = data.accepted(:);
        save(data.edFile, 'ed', '-v7.3');
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
%  LOCAL HELPERS (pure)
%  ========================================================================

function panels = defaultPanels(d)
P = {};
add = @(s, r, h, l) struct('source', s, 'region', r, 'height', h, 'label', l);
if d.hasStates,             P{end + 1} = add('hypnogram', 'wide', 0.28, 'State'); end
if ~isempty(d.specAdapter), P{end + 1} = add('spec', 'wide', 1.4, 'Freq (Hz)'); end
if ~isempty(d.emgRms),      P{end + 1} = add('emgRms', 'wide', 0.7, 'EMG RMS'); end
P{end + 1} = add('eventTicks', 'wide', 0.4, 'Events');
P{end + 1} = add('lfp', 'narrow', 1.2, 'LFP');
P{end + 1} = add('emg', 'narrow', 0.8, 'EMG');
if d.hasSpikes,             P{end + 1} = add('raster', 'narrow', 1.2, 'Units'); end
panels = [P{:}];
end

function panels = filterAvail(panels, d)
keep = true(1, numel(panels));
for i = 1:numel(panels)
    switch panels(i).source
        case 'raster',    keep(i) = d.hasSpikes;
        case 'hypnogram', keep(i) = d.hasStates;
        case 'spec',      keep(i) = ~isempty(d.specAdapter);
        case 'emgRms',    keep(i) = ~isempty(d.emgRms);
    end
end
panels = panels(keep);
end

function panels = addAxField(panels)
% guarantee an .ax field on every entry (and on an empty array), so panel
% struct arrays stay concatenable when the user adds panels later
if isempty(panels)
    panels = struct('source', {}, 'region', {}, 'height', {}, 'label', {}, 'ax', {});
    return;
end
for i = 1:numel(panels)
    panels(i).ax = gobjects(1);
end
end

function pn = makePanel(source, region)
pn = struct('source', source, 'region', region, ...
    'height', defaultHeight(source), 'label', srcLabel(source), 'ax', gobjects(1));
end

function s = availSources(~, d)
% every data-available source is offered in BOTH regions
s = {};
if ~isempty(d.specAdapter), s{end + 1} = 'spec'; end
if ~isempty(d.emgRms),      s{end + 1} = 'emgRms'; end
if d.hasStates,             s{end + 1} = 'hypnogram'; end
s{end + 1} = 'eventTicks';
s{end + 1} = 'lfp';
s{end + 1} = 'emg';
if d.hasSpikes,             s{end + 1} = 'raster'; end
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
% show x-tick labels + the axis label only on the region's outer axis (top of
% the Top region, bottom of the Bottom region), so the divider never hides it
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

function h = defaultHeight(src)
switch src
    case 'spec',       h = 1.4;
    case 'emgRms',     h = 0.7;
    case 'hypnogram',  h = 0.28;
    case 'eventTicks', h = 0.4;
    case 'lfp',        h = 1.2;
    case 'emg',        h = 0.8;
    case 'raster',     h = 1.2;
    otherwise,         h = 1.0;
end
end

function l = srcLabel(src)
switch src
    case 'spec',       l = 'Freq (Hz)';
    case 'emgRms',     l = 'EMG RMS';
    case 'hypnogram',  l = 'State';
    case 'eventTicks', l = 'Events';
    case 'lfp',        l = 'LFP';
    case 'emg',        l = 'EMG';
    case 'raster',     l = 'Units';
    otherwise,         l = src;
end
end

function mute(ax)
% make plotted children non-pickable so the axis ButtonDownFcn fires on click
ch = allchild(ax);
if ~isempty(ch), set(ch, 'HitTest', 'off'); end
end

% EOF
