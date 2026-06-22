function hFig = utypes_gui(varargin)

% UTYPES_GUI Interactive visualization of unit types.
%
% Three coordinated uifigure windows (scatter + traces + waveforms) that
% share selection and grouping. The scatter embeds tblGUI_scatHist; the
% traces and waveforms embed tblGUI_xy. A "Push Units" button (added to the
% scatter window's action area) saves the curated unit types via utypes_push.
%
% INPUT (Optional Key-Value Pairs):
%   basepaths    (cell array) Full paths to recording folders.
%   tblUnit      (table) Pre-computed unit table. If empty, loads it.
%   tAxis        (numeric) Time axis for the Traces window.
%
% OUTPUT:
%   hFig         (figure handle) Handle to the Scatter Plot figure.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addOptional(p, 'basepaths', {}, @(x) iscell(x));
addOptional(p, 'tblUnit', table(), @istable);
addParameter(p, 'tAxis', [], @isnumeric);

parse(p, varargin{:});
basepaths = p.Results.basepaths;
tblUnit = p.Results.tblUnit;
tAxis = p.Results.tAxis;

%% ========================================================================
%  LOAD DATA
%  ========================================================================

% If table empty, create preliminary table
if isempty(tblUnit)
    tblUnit = mcu_tblVivo('basepaths', basepaths);
end

% Add Waveform column
tWv = [];
if ~ismember('Waveform', tblUnit.Properties.VariableNames)
    try
        [tblWv, tWv] = swv_tbl('basepaths', basepaths, 'flgPlot', false);
        tblUnit.Waveform = tblWv.Waveform;
    catch
        warning('Could not load waveforms (swv_tbl failed). Proceeding without them.');
    end
else
    tWv = linspace(-0.75, 0.8, size(tblUnit.Waveform, 2));
end

%% ========================================================================
%  PARAMS
%  ========================================================================

xVar = 'TP';
yVar = 'BLidor';
szVar = 'FR';
grpVar = 'UnitType';
dotAlpha = 0.5;

cfg = mcu_cfg();
clr = cfg.clr.unit;

posScat  = [50, 400, 800, 600];
posTrace = [900, 550, 700, 450];
posWv    = [900, 50, 700, 450];

% Handles (assigned below; coordinator callbacks close over them)
hFigScat = [];
hFigTrace = [];
hFigWv = [];
guiReady = false;

%% ========================================================================
%  WINDOWS
%  ========================================================================

cbkScat = @(indices) onSelect(indices, 'scatter');
cbkGrp  = @(varName, activeCats, src) onGroupChange(varName, activeCats, src);

% --- WINDOW 1: SCATTER ---
hTabScat = uifigure('Name', 'Scatter Plot', 'Position', posScat);
hFigScat = tblGUI_scatHist(tblUnit, ...
    'xVar', xVar, 'yVar', yVar, 'szVar', szVar, 'grpVar', grpVar, ...
    'clr', clr, 'alpha', dotAlpha, ...
    'Parent', hTabScat, 'SelectionCallback', cbkScat, 'GroupByCallback', cbkGrp);

% --- WINDOW 2: TRACES ---
if ~isempty(tAxis)
    hTabTraces = uifigure('Name', 'Traces', 'Position', posTrace);
    cbkTrace = @(indices) onSelect(indices, 'traces');
    hFigTrace = tblGUI_xy(tAxis, tblUnit, 'Parent', hTabTraces, 'yVar', [], ...
        'SelectionCallback', cbkTrace, 'GroupByCallback', cbkGrp);
end

% --- WINDOW 3: WAVEFORMS ---
if ismember('Waveform', tblUnit.Properties.VariableNames)
    hTabWv = uifigure('Name', 'Waveforms', 'Position', posWv);
    cbkWv = @(indices) onSelect(indices, 'waveforms');
    hFigWv = tblGUI_xy(tWv, tblUnit, 'Parent', hTabWv, 'yVar', 'Waveform', ...
        'SelectionCallback', cbkWv, 'GroupByCallback', cbkGrp);
end

%% ========================================================================
%  ACTION BUTTON
%  ========================================================================

% Add "Push Units" to the scatter window's reserved action area (so it never
% collides with tblGUI_scatHist's own Select / Save buttons).
dScat = hFigScat.UserData;
tblgui.labeledControl(dScat.gActions, 'button', '', 'Text', 'Push Units', ...
    'ButtonPushedFcn', @(~, ~) onPushUnits(basepaths, hFigScat));

hFig = hFigScat;
guiReady = true;

% Initial cross-window sync to the scatter's grouping
try
    [~, allCats] = tblgui.selectedCats(dScat.chkGrp);
    onGroupChange(dScat.ddGrp.Value, allCats, hFigScat);
catch
end

%% ========================================================================
%  COORDINATOR CALLBACKS
%  ========================================================================

    function onSelect(indices, sourceName)
        if ~guiReady, return; end
        targetFigs = {hFigScat, hFigTrace, hFigWv};
        targetNames = {'scatter', 'traces', 'waveforms'};
        for i = 1:length(targetFigs)
            h = targetFigs{i};
            if isempty(h) || ~isvalid(h) || strcmp(sourceName, targetNames{i}), continue; end
            try
                data = h.UserData;
                if isfield(data, 'highlightFcn')
                    data.highlightFcn(indices);
                end
            catch
            end
        end
    end

    function onGroupChange(varName, activeCats, srcHandle)
        if ~guiReady, return; end
        targetFigs = {hFigScat, hFigTrace, hFigWv};
        for i = 1:length(targetFigs)
            h = targetFigs{i};
            if isempty(h) || ~isvalid(h) || h == srcHandle, continue; end
            try
                data = h.UserData;
                if isfield(data, 'setGroupVarFcn')
                    data.setGroupVarFcn(varName, activeCats);
                end
            catch
            end
        end
    end

end

function onPushUnits(basepaths, hContainer)
data = hContainer.UserData;
if isfield(data, 'tbl')
    utypes_push(basepaths, data.tbl);
    tblgui.notify(hContainer, 'Units saved successfully!', 'success');
else
    tblgui.notify(hContainer, 'Could not retrieve table data from GUI.', 'error');
end
end
