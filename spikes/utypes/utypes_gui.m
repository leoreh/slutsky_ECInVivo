function hFig = utypes_gui(varargin)

% UTYPES_GUI Interactive visualization of unit types.
%
% This GUI provides a comprehensive view of neuronal unit classification,
% featuring three synchronized WINDOWS:
%   1. Scatter Plot: Feature distribution visualization (e.g. TP vs Lidor).
%   2. Traces: Functional responses or firing rate traces over time.
%   3. Waveforms: Average spike waveforms.
%
% Selection AND Grouping are synchronized across all windows.
%
% INPUT (Optional Key-Value Pairs):
%   basepaths    (cell array) Full paths to recording folders.
%   uTbl         (table) Pre-computed unit table. If provided, loading is skipped.
%   flgSave      (logical) If true, adds a "Push Units" button to save classification.
%   cfg          (struct) Configuration details for tblGUI_scatHist.
%   tAxis        (numeric) Time axis for the Traces tab.
%
% OUTPUT:
%   hFig         (figure handle) Handle to the Scatter Plot figure (main).
%
% DEPENDENCIES:
%   basepaths2vars, v2tbl, tblGUI_scatHist, tblGUI_xy, plot_wv
%

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addOptional(p, 'basepaths', {}, @(x) iscell(x));
addOptional(p, 'uTbl', table(), @istable);
addOptional(p, 'flgSave', false, @islogical);
addOptional(p, 'cfg', struct(), @isstruct);
addParameter(p, 'tAxis', [], @isnumeric);

parse(p, varargin{:});
basepaths = p.Results.basepaths;
uTbl = p.Results.uTbl;
flgSave = p.Results.flgSave;
cfg = p.Results.cfg;
tAxis = p.Results.tAxis;

%% ========================================================================
%  LOAD DATA
%  ========================================================================

if isempty(uTbl)

    % Load required variables
    vars = {'swv_metrics', 'st_metrics', 'fr', 'units'};
    v = basepaths2vars('basepaths', basepaths, 'vars', vars);

    % Prepare table metadata
    tagFiles = struct();
    tagFiles.Name = get_mname(basepaths);
    [~, fNames] = fileparts(basepaths);
    tagFiles.File = fNames;

    % VarMap (Features)
    varMap = struct();
    varMap.mfr = 'fr.mfr';
    varMap.royer = 'st.royer';
    varMap.lidor = 'st.lidor';
    varMap.tp = 'swv.tp';
    varMap.asym = 'swv.asym';
    varMap.hpk = 'swv.hpk';

    % Create table using v2tbl
    uTbl = v2tbl('v', v, 'varMap', varMap, 'tagFiles', tagFiles);

    % Assign UnitType from loaded units struct
    units = catfields([v.units], 2, true);
    uTbl.UnitType = vertcat(units.type{:});

end

%% ========================================================================
%  PLOT PARAMS
%  ========================================================================

% Default plotting configuration
if ~isfield(cfg, 'xVar'), cfg.xVar = 'TP'; end
if ~isfield(cfg, 'yVar'), cfg.yVar = 'BLidor'; end
if ~isfield(cfg, 'szVar'), cfg.szVar = 'FR'; end
if ~isfield(cfg, 'grpVar'), cfg.grpVar = 'UnitType'; end
if ~isfield(cfg, 'alpha'), cfg.alpha = 0.5; end

% Default colors from mcu_cfg
if ~isfield(cfg, 'clr')
    cfgMcu = mcu_cfg();
    cfg.clr = cfgMcu.clr.unit([3 : -1 : 1], :);
end

%% ========================================================================
%  PLOT - 3 Separate Figure Windows
%  ========================================================================

screenSz = get(0, 'ScreenSize');
% width = screenSz(3); height = screenSz(4);
% Define positions
posScat = [50, 400, 800, 600];
posTrace = [900, 550, 700, 450];
posWv = [900, 50, 700, 450];

% Initialize handles to prevent callback errors during creation
hFigScat = [];
hFigTrace = [];
hFigWv = [];
guiReady = false;

% --- WINDOW 1: SCATTER ---
hTabScat = figure('Name', 'Scatter Plot', 'NumberTitle', 'off', ...
    'Position', posScat, 'MenuBar', 'none', 'ToolBar', 'figure');

cbkScat = @(indices) onSelect(indices, 'scatter');
cbkGrp = @(varName, activeCats, src) onGroupChange(varName, activeCats, src);

hFigScat = tblGUI_scatHist(uTbl, 'cfg', cfg, ...
    'Parent', hTabScat, 'SelectionCallback', cbkScat, ...
    'GroupByCallback', cbkGrp); % Passed to tblGUI_scatHist if modified

% --- WINDOW 2: TRACES ---
hFigTrace = [];
if ~isempty(tAxis)
    hTabTraces = figure('Name', 'Traces', 'NumberTitle', 'off', ...
        'Position', posTrace, 'MenuBar', 'none', 'ToolBar', 'figure');

    cbkTrace = @(indices) onSelect(indices, 'traces');

    hFigTrace = tblGUI_xy(tAxis, uTbl, 'Parent', hTabTraces, 'yVar', [], ...
        'SelectionCallback', cbkTrace, 'GroupByCallback', cbkGrp);
end

% --- WINDOW 3: WAVEFORMS ---
hTabWv = figure('Name', 'Waveforms', 'NumberTitle', 'off', ...
    'Position', posWv, 'MenuBar', 'none', 'ToolBar', 'figure');

cbkWv = @(indices) onSelect(indices, 'waveforms');

[~, hFigWv] = plot_wv('basepaths', basepaths, 'UnitType', uTbl.UnitType, ...
    'Parent', hTabWv, 'SelectionCallback', cbkWv, 'GroupByCallback', cbkGrp);

% -------------------------------------------------------------------------

% Button for Saving (on Scatter)
if flgSave
    uicontrol('Parent', hTabScat, 'Style', 'pushbutton', ... % hFigScat is container
        'String', 'Push Units', ...
        'Units', 'normalized', 'Position', [0.01, 0.11, 0.18, 0.05], ...
        'Callback', @(src, evt) onPushUnits(src, basepaths, hFigScat));
end

hFig = hFigScat; % Return main handle
guiReady = true; % Enable callbacks

% Force initial sync to Scatter's Config
% Get initial Grp from Scatter user data if available
try
    d = hFigScat.UserData;
    idx = get(d.ddGrp, 'Value');
    items = get(d.ddGrp, 'String');
    initGrp = items{idx};

    % Gather initial categories
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
        % Sync Selection
        fprintf('Selection from: %s\n', sourceName);

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
        % Sync Grouping Variable
        fprintf('Group Change to: %s\n', varName);

        targetFigs = {hFigScat, hFigTrace, hFigWv};

        for i = 1:length(targetFigs)
            h = targetFigs{i};
            if isempty(h) || ~isvalid(h), continue; end

            % Skip if this handle is the source (or contains the source)
            % srcHandle might be the container or a child of it.
            % But 'h' is the figure/container stored.
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
    push_units(basepaths, fetTbl);
    msgbox('Units saved successfully!', 'Success');
else
    errordlg('Could not retrieve table data from GUI.', 'Error');
end
end

% EOF
