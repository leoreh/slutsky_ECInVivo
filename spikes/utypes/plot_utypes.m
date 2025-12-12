function hFig = plot_utypes(varargin)

% PLOT_UTYPES Visualization of unit types using tblGUI_scatHist
%
% This function replaces the legacy waveform plotting. It loads unit data
% (metrics and classification) and launches the interactive tblGUI_scatHist.
%
% INPUT (Optional Key-Value Pairs):
%   basepaths    (cell array) Full paths to recording folders.
%   uTbl         (table) Pre-computed unit table. If provided, loading is skipped.
%   flgSave      (logical) If true, adds a "Push Units" button to save classification.
%   cfg          (struct) Configuration details for tblGUI_scatHist.
%
% OUTPUT:
%   hFig         (figure handle) Handle to the GUI figure.
%
% DEPENDENCIES:
%   basepaths2vars, v2tbl, tblGUI_scatHist
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

if isempty(tAxis)
    % Standard single plot mode
    hFig = tblGUI_scatHist(uTbl, 'cfg', cfg);
    hScatContainer = hFig;
else
    % Dual tab mode with tAxis
    hFig = figure('Name', 'Unit Types & Traces', 'NumberTitle', 'off', ...
        'Position', [100, 100, 1200, 800], 'MenuBar', 'none', 'ToolBar', 'figure');

    hTabGroup = uitabgroup(hFig);
    hTab1 = uitab(hTabGroup, 'Title', 'Scatter');
    hTab2 = uitab(hTabGroup, 'Title', 'Traces');

    % --- TAB 1: Scatter ---
    % Define callback for linking
    cbk = @(indices) onSelectUnit(indices, hTab2);

    hScatContainer = tblGUI_scatHist(uTbl, 'cfg', cfg, ...
        'Parent', hTab1, 'SelectionCallback', cbk);

    % --- TAB 2: XY Traces ---
    % We assume uTbl has 'frMat' or similar suitable for plotting.
    % tblGUI_xy auto-selects if yVar is empty, or uses cfg if we passed it?
    % Let's try to be smart or default to empty yVar (auto-select).
    % If cfg has yVar, it might be 'lidor' (scalar) which is wrong for XY.
    % So we force yVar to empty or let user/auto handle it.
    % Actually, we should probably check if uTbl has 'frMat' or 'fr_strd'.

    tblGUI_xy(tAxis, uTbl, 'Parent', hTab2, 'yVar', []);

end

if flgSave
    % Add Push/Save button to the figure
    % Position it above the existing "Save Table" button in tblGUI_scatHist
    % We parent it to the scatter container (Figure or Tab)
    uicontrol('Parent', hScatContainer, 'Style', 'pushbutton', ...
        'String', 'Push Units', ...
        'Units', 'normalized', 'Position', [0.01, 0.16, 0.18, 0.05], ...
        'Callback', @(src, evt) onPushUnits(src, basepaths, hScatContainer));
end

end

function onSelectUnit(indices, hTabXY)
% Callback from Scatter to highlight XY
try
    data = hTabXY.UserData;
    if isfield(data, 'highlightFcn')
        data.highlightFcn(indices);
    end
catch ME
    warning(ME.identifier, 'Failed to highlight traces: %s', ME.message);
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
