function hFig = plot_utypes(varargin)

% PLOT_UTYPES Visualization of unit types using plot_tblGUI
%
% This function replaces the legacy waveform plotting. It loads unit data
% (metrics and classification) and launches the interactive plot_tblGUI.
%
% INPUT (Optional Key-Value Pairs):
%   basepaths    (cell array) Full paths to recording folders.
%   uTbl         (table) Pre-computed unit table. If provided, loading is skipped.
%   cfg          (struct) Configuration details for plot_tblGUI.
%
% OUTPUT:
%   hFig         (figure handle) Handle to the GUI figure.
%
% DEPENDENCIES:
%   basepaths2vars, v2tbl, plot_tblGUI
%

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addOptional(p, 'basepaths', {}, @(x) iscell(x));
addOptional(p, 'uTbl', table(), @istable);
addOptional(p, 'cfg', struct(), @isstruct);

parse(p, varargin{:});
basepaths = p.Results.basepaths;
uTbl = p.Results.uTbl;
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

hFig = plot_tblGUI(uTbl, 'cfg', cfg);

end

% EOF
