function sel = gui_loadDialog(parent, basepath)
% GUI_LOADDIALOG  Progressive "Load data" dialog for guiPath_curate.
%
%   sel = gui_loadDialog(parent, basepath) opens a small always-on-top
%   dialog that reveals inputs as the user chooses, and returns a selection
%   struct (or [] on cancel). The flow is: pick a TYPE (what the data is), then a
%   SOURCE (where it comes from); only then does the source-specific input
%   appear - a base-workspace variable list (editable, so a name or var.field can
%   be typed), a file browser plus its variable list, or a channel box. An fs box
%   shows only for a raw trace; a Top/Bottom placement only for a signal.
%
%   OUTPUT sel (fields):
%       .type   'trace' | 'spec' | 'hypnogram' | 'raster' | 'eventTicks' | 'stateStrip'
%       .from   'ws' | 'file' | 'bin'
%       .value  for ws  -> the variable name (optionally var.field) [char]
%               for bin -> the channel spec (e.g. '5' or '[5 6 7]')  [char]
%               for file-> the loaded value (inline)
%       .var    variable name (ws / file), '' for bin
%       .file   full path (file), '' otherwise
%       .fs     sampling rate [Hz] for a trace, else []
%       .region 'top' | 'bottom'
%
%   The host (guiPath_curate.onLoadUnified) turns sel into a loadCore call. Curation
%   types (events / states) cannot come from a binary channel.
%
%   See also guiPath_curate, guiPath_load, gui_chooseDialog.
%
%   HISTORY:
%       05 Jul 2026 - progressive Load dialog (replaces the always-on Load fields).
%       05 Jul 2026 - dropped WindowStyle 'modal' (its input grab could outlive
%                     teardown and freeze the host); explicit delete + drawnow.

narginchk(1, 2);
if nargin < 2, basepath = ''; end
sel = [];

st = struct('file', '', 'fileVars', {{}}, 'session', [], 'binFile', '');   % captured by callbacks
dyn = struct();

d = uifigure('Name', 'Load data', 'Position', [100, 100, 400, 300]);
% kept 'normal', NOT 'modal': a modal uifigure's input grab can outlive the
% window during teardown and collide with the alert the caller shows next,
% leaving the main figure frozen. uiwait already blocks the caller's code.
try, d.WindowStyle = 'alwaysontop'; catch, end
host = ancestor(parent, 'figure');
if ~isempty(host) && isvalid(host)                          % centre over the host
    pr = host.Position;
    d.Position = [pr(1) + (pr(3) - 400) / 2, pr(2) + (pr(4) - 300) / 2, 400, 300];
else
    movegui(d, 'center');
end
cleaner = onCleanup(@() delete(d(isvalid(d))));

gl = uigridlayout(d, [6, 2], 'RowHeight', {'fit', 'fit', 'fit', 'fit', 'fit', 'fit'}, ...
    'ColumnWidth', {110, '1x'}, 'Padding', 12, 'RowSpacing', 8, 'ColumnSpacing', 6);

lt = uilabel(gl, 'Text', 'Type', 'FontWeight', 'bold'); lt.Layout.Row = 1; lt.Layout.Column = 1;
ddType = uidropdown(gl, 'Items', {'Trace', 'Spectrogram', 'Hypnogram', 'Raster', 'Events', 'States'}, ...
    'ValueChangedFcn', @(~, ~) onType());
ddType.Layout.Row = 1; ddType.Layout.Column = 2;

lf = uilabel(gl, 'Text', 'Source', 'FontWeight', 'bold'); lf.Layout.Row = 2; lf.Layout.Column = 1;
ddFrom = uidropdown(gl, 'Items', {'Workspace', 'File', 'Binary channel'}, ...
    'ValueChangedFcn', @(~, ~) onFrom());
ddFrom.Layout.Row = 2; ddFrom.Layout.Column = 2;

dynP = uipanel(gl, 'BorderType', 'none'); dynP.Layout.Row = 3; dynP.Layout.Column = [1, 2];

lfs = uilabel(gl, 'Text', 'fs (Hz)'); lfs.Layout.Row = 4; lfs.Layout.Column = 1;
edFs = uieditfield(gl, 'numeric', 'Value', 1250, 'Limits', [eps, Inf]);
edFs.Layout.Row = 4; edFs.Layout.Column = 2;

lreg = uilabel(gl, 'Text', 'Place in'); lreg.Layout.Row = 5; lreg.Layout.Column = 1;
ddReg = uidropdown(gl, 'Items', {'Top', 'Bottom'}, 'Value', 'Bottom');
ddReg.Layout.Row = 5; ddReg.Layout.Column = 2;

bg = uigridlayout(gl, [1, 2], 'Padding', 0, 'ColumnWidth', {'1x', '1x'}, 'ColumnSpacing', 6);
bg.Layout.Row = 6; bg.Layout.Column = [1, 2];
uibutton(bg, 'Text', 'Load', 'ButtonPushedFcn', @(~, ~) onOk());
uibutton(bg, 'Text', 'Cancel', 'ButtonPushedFcn', @(~, ~) uiresume(d));

onType();                              % initialise dependent controls
d.CloseRequestFcn = @(~, ~) uiresume(d);
uiwait(d);                             % sel is set by onOk (else stays [])
if isvalid(d), delete(d); end          % tear the dialog down NOW, then flush, so
drawnow;                               % the caller's next alert never races it

% =====================================================================
    function onType()
        isCur = any(strcmp(ddType.Value, {'Events', 'States'}));
        keepItems(ddFrom, ifelse(isCur, {'Workspace', 'File'}, {'Workspace', 'File', 'Binary channel'}));
        setRow(lfs, edFs, strcmp(ddType.Value, 'Trace'));   % fs only for a raw trace
        setRow(lreg, ddReg, ~isCur);                        % placement only for a signal
        onFrom();
    end

    function onFrom()
        rebuildDyn();
        if strcmp(ddFrom.Value, 'Binary channel')
            ensureSession();
            if ~isempty(st.session), edFs.Value = st.session.extracellular.srLfp; end
        end
    end

    function rebuildDyn()
        delete(dynP.Children);
        g = uigridlayout(dynP, [3, 2], 'RowHeight', {'fit', 'fit', 'fit'}, 'ColumnWidth', {110, '1x'}, ...
            'Padding', [0, 4, 0, 4], 'RowSpacing', 6, 'ColumnSpacing', 6);
        switch ddFrom.Value
            case 'Workspace'
                a = uilabel(g, 'Text', 'Variable'); a.Layout.Row = 1; a.Layout.Column = 1;
                vars = wsVars();
                dyn.wsVar = uidropdown(g, 'Editable', 'on', 'Items', vars, 'Value', firstOr(vars));
                dyn.wsVar.Layout.Row = 1; dyn.wsVar.Layout.Column = 2;
                hint = uilabel(g, 'Text', 'pick, or type  var  or  var.field', 'FontColor', [0.5 0.5 0.5]);
                hint.Layout.Row = 2; hint.Layout.Column = 2;
            case 'File'
                a = uilabel(g, 'Text', 'File'); a.Layout.Row = 1; a.Layout.Column = 1;
                dyn.fileLbl = uilabel(g, 'Text', fileShort(st.file));
                dyn.fileLbl.Layout.Row = 1; dyn.fileLbl.Layout.Column = 2;
                b = uibutton(g, 'Text', 'Browse...', 'ButtonPushedFcn', @(~, ~) onBrowse());
                b.Layout.Row = 2; b.Layout.Column = 1;
                dyn.fileVar = uidropdown(g, 'Items', orEmpty(st.fileVars));
                dyn.fileVar.Layout.Row = 2; dyn.fileVar.Layout.Column = 2;
            otherwise    % Binary channel
                a = uilabel(g, 'Text', 'Channel(s)'); a.Layout.Row = 1; a.Layout.Column = 1;
                dyn.binCh = uieditfield(g, 'text', 'Value', '', 'Placeholder', 'e.g. 5 or [5 6 7]');
                dyn.binCh.Layout.Row = 1; dyn.binCh.Layout.Column = 2;
                bb = uibutton(g, 'Text', 'Browse...', 'ButtonPushedFcn', @(~, ~) onBrowseBin());
                bb.Layout.Row = 2; bb.Layout.Column = 1;
                dyn.binLbl = uilabel(g, 'Text', binShort());
                dyn.binLbl.Layout.Row = 2; dyn.binLbl.Layout.Column = 2;
                ensureSession();
                if isempty(st.session)
                    txt = 'session.mat not found (needs nCh / fs)';
                else
                    txt = sprintf('%d channels @ %g Hz', st.session.extracellular.nChannels, ...
                        st.session.extracellular.srLfp);
                end
                hint = uilabel(g, 'Text', txt, 'FontColor', [0.5 0.5 0.5]);
                hint.Layout.Row = 3; hint.Layout.Column = 2;
        end
    end

    function onBrowse()
        [fn, fp] = uigetfile('*.mat', 'Load .mat');
        if isequal(fn, 0), return; end
        st.file = fullfile(fp, fn);
        try, info = whos('-file', st.file); st.fileVars = {info.name}; catch, st.fileVars = {}; end
        rebuildDyn();
    end

    function onBrowseBin()
        [fn, fp] = uigetfile({'*.lfp;*.dat;*.bin', 'Binary (*.lfp, *.dat, *.bin)'; '*.*', 'All files'}, ...
            'Select a binary file');
        if isequal(fn, 0), return; end
        st.binFile = fullfile(fp, fn);
        if isfield(dyn, 'binLbl') && isvalid(dyn.binLbl), dyn.binLbl.Text = binShort(); end
    end

    function s = binShort()
        if isempty(st.binFile), s = '<basename>.lfp (default)'; else, [~, n, e] = fileparts(st.binFile); s = [n, e]; end
    end

    function onOk()
        s = struct('type', typeCode(ddType.Value), 'from', '', 'value', [], ...
            'var', '', 'file', '', 'fs', [], 'region', lower(ddReg.Value));
        if strcmp(ddType.Value, 'Trace'), s.fs = edFs.Value; end
        switch ddFrom.Value
            case 'Workspace'
                nm = strtrim(dyn.wsVar.Value);
                if isempty(nm), uialert(d, 'Pick or type a variable name.', 'Load'); return; end
                s.from = 'ws'; s.value = nm; s.var = nm;
            case 'File'
                if isempty(st.file), uialert(d, 'Browse to a .mat file.', 'Load'); return; end
                vn = dyn.fileVar.Value;
                if isempty(vn), uialert(d, 'Pick a variable in the file.', 'Load'); return; end
                try, S = load(st.file, vn); catch ME, uialert(d, ME.message, 'Load'); return; end
                s.from = 'file'; s.value = S.(vn); s.var = vn; s.file = st.file;
            otherwise
                ch = strtrim(dyn.binCh.Value);
                if isempty(ch) || isempty(str2num(ch)) %#ok<ST2NM>
                    uialert(d, 'Type a channel number (e.g. 5 or [5 6 7]).', 'Load'); return;
                end
                s.from = 'bin';
                if isempty(st.binFile), s.value = ch; else, s.value = [ch, '>', st.binFile]; end
        end
        sel = s;
        uiresume(d);
    end

    function ensureSession()
        if ~isempty(st.session) || isempty(basepath), return; end
        try
            v = basepaths2vars('basepaths', {basepath}, 'vars', {'session'}, 'flgPrnt', false);
            if isfield(v, 'session'), st.session = v.session; end
        catch
        end
    end
end     % MAIN

% =====================================================================
%  HELPERS (pure)
% =====================================================================
function keepItems(dd, items)
% set a dropdown's Items, keeping the current Value if still valid, else first
dd.Items = items;
if ~any(strcmp(dd.Value, items)), dd.Value = items{1}; end
end

function setRow(lbl, ctrl, tf)
% show/hide a label + its control (the 'fit' row collapses when both are hidden)
v = 'off'; if tf, v = 'on'; end
lbl.Visible = v; ctrl.Visible = v;
end

function vars = wsVars()
try, vars = evalin('base', 'who'); catch, vars = {}; end
vars = vars(:)';
if isempty(vars), vars = {''}; end
end

function s = firstOr(c)
if isempty(c), s = ''; else, s = c{1}; end
end

function c = orEmpty(c)
if isempty(c), c = {''}; end
end

function s = fileShort(f)
if isempty(f), s = '(none)'; else, [~, n, e] = fileparts(f); s = [n, e]; end
end

function out = ifelse(tf, a, b)
if tf, out = a; else, out = b; end
end

function t = typeCode(label)
switch label
    case 'Trace',       t = 'trace';
    case 'Spectrogram', t = 'spec';
    case 'Hypnogram',   t = 'hypnogram';
    case 'Raster',      t = 'raster';
    case 'Events',      t = 'eventTicks';
    case 'States',      t = 'stateStrip';
    otherwise,          t = 'trace';
end
end
