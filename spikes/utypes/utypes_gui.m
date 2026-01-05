function hFig = utypes_gui(varargin)

% UTYPES_GUI Interactive visualization of unit types.
%
% INPUT (Optional Key-Value Pairs):
%   basepaths    (cell array) Full paths to recording folders.
%   tblUnit      (table) Pre-computed unit table. If empty, loads it.
%   tAxis        (numeric) Time axis for the Traces tab.
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
    tblUnit = mcu_tblVivo(basepaths);
end

% Add Waveform column
if ~ismember('Waveform', tblUnit.Properties.VariableNames)
    try
        [tblWv, tWv] = swv_tbl('basepaths', basepaths, 'flgPlot', false);
        tblUnit.Waveform = tblWv.Waveform;
    catch
        warning('Could not load waveforms (swv_tbl failed). Proceeding without them.');
    end
else
    % If table passed with Waveform but we need time axis
    tWv = linspace(-0.75, 0.8, size(tblUnit.Waveform, 2));
end

%% ========================================================================
%  PARAMS
%  ========================================================================

% Hardcoded Configuration for scatter plot
% Hardcoded Configuration for scatter plot
xVar = 'TP';
yVar = 'BLidor';
szVar = 'FR';
grpVar = 'UnitType';
dotAlpha = 0.5;

% Colors from mcu_cfg
cfg = mcu_cfg();
clr = cfg.clr.unit;

%% ========================================================================
%  PLOT
%  ========================================================================

posScat = [50, 400, 800, 600];
posTrace = [900, 550, 700, 450];
posWv = [900, 50, 700, 450];

% Initialize handles
hFigScat = [];
hFigTrace = [];
hFigWv = [];
guiReady = false;

% --- WINDOW 1: SCATTER ---
hTabScat = figure('Name', 'Scatter Plot', 'NumberTitle', 'off', ...
    'Position', posScat, 'MenuBar', 'none', 'ToolBar', 'figure');

cbkScat = @(indices) onSelect(indices, 'scatter');
cbkGrp = @(varName, activeCats, src) onGroupChange(varName, activeCats, src);

hFigScat = tblGUI_scatHist(tblUnit, ...
    'xVar', xVar, 'yVar', yVar, 'szVar', szVar, 'grpVar', grpVar, ...
    'clr', clr, 'alpha', dotAlpha, ...
    'Parent', hTabScat, 'SelectionCallback', cbkScat, ...
    'GroupByCallback', cbkGrp);

% --- WINDOW 2: TRACES ---
if ~isempty(tAxis)
    hTabTraces = figure('Name', 'Traces', 'NumberTitle', 'off', ...
        'Position', posTrace, 'MenuBar', 'none', 'ToolBar', 'figure');

    cbkTrace = @(indices) onSelect(indices, 'traces');

    hFigTrace = tblGUI_xy(tAxis, tblUnit, 'Parent', hTabTraces, 'yVar', [], ...
        'SelectionCallback', cbkTrace, 'GroupByCallback', cbkGrp);
end

% --- WINDOW 3: WAVEFORMS ---
if ismember('Waveform', tblUnit.Properties.VariableNames)
    hTabWv = figure('Name', 'Waveforms', 'NumberTitle', 'off', ...
        'Position', posWv, 'MenuBar', 'none', 'ToolBar', 'figure');

    cbkWv = @(indices) onSelect(indices, 'waveforms');

    % Creates waveform gui directly by calling tblGUI_xy
    hFigWv = tblGUI_xy(tWv, tblUnit, 'Parent', hTabWv, 'yVar', 'Waveform', ...
        'SelectionCallback', cbkWv, 'GroupByCallback', cbkGrp);
end


% -------------------------------------------------------------------------

% Button for Saving (on Scatter)
uicontrol('Parent', hTabScat, 'Style', 'pushbutton', ...
    'String', 'Push Units', ...
    'Units', 'normalized', 'Position', [0.01, 0.11, 0.18, 0.05], ...
    'Callback', @(src, evt) onPushUnits(src, basepaths, hFigScat));

hFig = hFigScat; % Return main handle
guiReady = true; % Enable callbacks

% Force initial sync to Scatter's Config
try
    d = hFigScat.UserData;
    idx = get(d.ddGrp, 'Value');
    items = get(d.ddGrp, 'String');
    initGrp = items{idx};

    if isfield(d, 'chkGrp') && ~isempty(d.chkGrp)
        selectedIdx = arrayfun(@(x) get(x, 'Value'), d.chkGrp);
        allCats = arrayfun(@(x) string(get(x, 'String')), d.chkGrp);
        activeCats = cellstr(allCats(logical(selectedIdx)));
        onGroupChange(initGrp, activeCats, hFigScat);
    end
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
            name = targetNames{i};

            if isempty(h) || ~isvalid(h) || strcmp(sourceName, name), continue; end

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
            if isempty(h) || ~isvalid(h), continue; end

            if h == srcHandle, continue; end

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

function onPushUnits(src, basepaths, hContainer)
data = hContainer.UserData;
if isfield(data, 'tbl')
    fetTbl = data.tbl;
    utypes_push(basepaths, fetTbl);
    msgbox('Units saved successfully!', 'Success');
else
    errordlg('Could not retrieve table data from GUI.', 'Error');
end
end
