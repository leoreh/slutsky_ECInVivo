function hFig = plot_utypes(varargin)

% PLOT_UTYPES Visualization of unit types using plot_tblScatter
%
% This function replaces the legacy waveform plotting. It loads unit data
% (metrics and classification) and launches the interactive plot_tblScatter.
%
% INPUT (Optional Key-Value Pairs):
%   basepaths    (cell array) Full paths to recording folders.
%   uTbl         (table) Pre-computed unit table. If provided, loading is skipped.
%   flgSave      (logical) If true, adds a "Push Units" button to save classification.
%   cfg          (struct) Configuration details for plot_tblScatter.
%
% OUTPUT:
%   hFig         (figure handle) Handle to the GUI figure.
%
% DEPENDENCIES:
%   basepaths2vars, v2tbl, plot_tblScatter
%

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addOptional(p, 'basepaths', {}, @(x) iscell(x));
addOptional(p, 'uTbl', table(), @istable);
addOptional(p, 'flgSave', false, @islogical);
addOptional(p, 'cfg', struct(), @isstruct);

parse(p, varargin{:});
basepaths = p.Results.basepaths;
uTbl = p.Results.uTbl;
flgSave = p.Results.flgSave;
cfg = p.Results.cfg;

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
if ~isfield(cfg, 'xVar'), cfg.xVar = 'tp'; end
if ~isfield(cfg, 'yVar'), cfg.yVar = 'lidor'; end
if ~isfield(cfg, 'szVar'), cfg.szVar = 'mfr'; end
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

hFig = plot_tblScatter(uTbl, 'cfg', cfg);

if flgSave
    % Add Push/Save button to the figure
    % Position it above the existing "Save Table" button in plot_tblScatter
    uicontrol('Parent', hFig, 'Style', 'pushbutton', ...
        'String', 'Push Units', ...
        'Units', 'normalized', 'Position', [0.01, 0.16, 0.18, 0.05], ...
        'Callback', @(src, evt) onPushUnits(src, basepaths));
end

end

function onPushUnits(src, basepaths)
hFig = ancestor(src, 'figure');
data = guidata(hFig);

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
