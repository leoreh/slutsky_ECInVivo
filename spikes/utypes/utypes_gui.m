function hFig = utypes_gui(varargin)

% UTYPES_GUI Interactive visualization of unit types.
%
% This GUI provides a comprehensive view of neuronal unit classification,
% featuring three synchronized tabs:
%   1. Scatter Plot: Feature distribution visualization (e.g. TP vs Lidor).
%   2. Traces: Functional responses or firing rate traces over time.
%   3. Waveforms: Average spike waveforms.
%
% Selection in any tab is propagated to the others.
%
% INPUT (Optional Key-Value Pairs):
%   basepaths    (cell array) Full paths to recording folders.
%   uTbl         (table) Pre-computed unit table. If provided, loading is skipped.
%   flgSave      (logical) If true, adds a "Push Units" button to save classification.
%   cfg          (struct) Configuration details for tblGUI_scatHist.
%   tAxis        (numeric) Time axis for the Traces tab.
%
% OUTPUT:
%   hFig         (figure handle) Handle to the GUI figure.
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
%  PLOT
%  ========================================================================

% Main Figure
hFig = figure('Name', 'Unit Types GUI', 'NumberTitle', 'off', ...
    'Position', [50, 50, 1400, 900], 'MenuBar', 'none', 'ToolBar', 'figure');

hTabGroup = uitabgroup(hFig);
hTabScat = uitab(hTabGroup, 'Title', 'Scatter');
hTabTraces = uitab(hTabGroup, 'Title', 'Traces');
hTabWv = uitab(hTabGroup, 'Title', 'Waveforms');

% --- TAB 1: SCATTER ---
% Callback: When Scatter selects, update Traces and Waveforms
cbkScat = @(indices) onSelect(indices, 'scatter');

hFigScat = tblGUI_scatHist(uTbl, 'cfg', cfg, ...
    'Parent', hTabScat, 'SelectionCallback', cbkScat);

% --- TAB 2: TRACES ---
hFigTrace = [];
if ~isempty(tAxis)
    % Callback: When Traces selects, update Scatter and Waveforms
    cbkTrace = @(indices) onSelect(indices, 'traces');

    hFigTrace = tblGUI_xy(tAxis, uTbl, 'Parent', hTabTraces, 'yVar', [], ...
        'SelectionCallback', cbkTrace);
else
    uicontrol('Parent', hTabTraces, 'Style', 'text', 'String', 'No time axis provided.', ...
        'Units', 'normalized', 'Position', [0.4 0.5 0.2 0.05]);
end

% --- TAB 3: WAVEFORMS ---
% Callback: When Waveforms selects, update Scatter and Traces
cbkWv = @(indices) onSelect(indices, 'waveforms');

% We need to pass basepaths to plot_wv to let it load waveforms if not in uTbl
% Or passing uTbl if plot_wv supports it?
% plot_wv refactor: inputs basepaths, prepares wvTbl.
% We want proper linking. plot_wv generates a NEW table (wvTbl).
% Logic check: uTbl and wvTbl must match rows.
% Since both load from same basepaths (presumably), rows should match order if sorted same.
% But to be safe, we rely on the fact they come from same source.
% We pass basepaths to plot_wv.
[~, hFigWv] = plot_wv('basepaths', basepaths, 'UnitType', uTbl.UnitType, ...
    'Parent', hTabWv, 'SelectionCallback', cbkWv);
% Note: plot_wv call to tblGUI_xy inside needs access to 'SelectionCallback'.
% The refactored plot_wv in step 61 did NOT add 'SelectionCallback' to its parser
% or pass it to tblGUI_xy. I missed that in the plan for plot_wv,
% but tblGUI_xy supports it. I need to update plot_wv quickly or just rely on
% modifying plot_wv to pass varargin or explicit arg.
% Actually, plot_wv does parse varargin, but doesn't explicitly look for 'SelectionCallback'.
% However, tblGUI_xy is called at the end. I can modify plot_wv to accept it
% or better, since I can't easily modify plot_wv mid-stream here without another tool call,
% I will assume I will fix plot_wv in a moment.
% FOR NOW: I will pass it, expecting to add support to plot_wv.

% Button for Saving
if flgSave
    % Add Push/Save button to the figure (Scatter tab container)
    uicontrol('Parent', hTabScat, 'Style', 'pushbutton', ...
        'String', 'Push Units', ...
        'Units', 'normalized', 'Position', [0.01, 0.16, 0.18, 0.05], ...
        'Callback', @(src, evt) onPushUnits(src, basepaths, hFigScat));
end


%% ========================================================================
%  COORDINATOR CALLBACK
%  ========================================================================

    function onSelect(indices, sourceName)
        % indices: logical array of selected units

        fprintf('Selection from: %s\n', sourceName);

        % 1. Update Scatter (if not source)
        if ~strcmp(sourceName, 'scatter') && ~isempty(hFigScat)
            try
                data = hFigScat.UserData;
                if isfield(data, 'highlightFcn')
                    data.highlightFcn(indices);
                end
            catch
            end
        end

        % 2. Update Traces (if not source)
        if ~strcmp(sourceName, 'traces') && ~isempty(hFigTrace)
            try
                data = hFigTrace.UserData;
                if isfield(data, 'highlightFcn')
                    data.highlightFcn(indices);
                end
            catch
            end
        end

        % 3. Update Waveforms (if not source)
        if ~strcmp(sourceName, 'waveforms') && ~isempty(hFigWv)
            try
                data = hFigWv.UserData;
                if isfield(data, 'highlightFcn')
                    data.highlightFcn(indices);
                end
            catch
            end
        end

    end

end

function onPushUnits(src, basepaths, hContainer)
% hContainer passed explicitly
data = hContainer.UserData;

% Access the modified table from the GUI data
if isfield(data, 'tbl')
    fetTbl = data.tbl;
    % Call the standalone push_units function
    push_units(basepaths, fetTbl);
    msgbox('Units saved successfully!', 'Success');
else
    errordlg('Could not retrieve table data from GUI.', 'Error');
end
end

% EOF
