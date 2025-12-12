function [wvTbl, hFig] = plot_wv(varargin)

% PLOT_WV displays average spike waveforms using tblGUI_xy.
% It loads waveform data, prepares a table, and launches the interactive GUI.
%
% INPUT (Optional Key-Value Pairs):
%   basepaths    (cell array) Full paths to recording folders. If empty,
%                uses current directory.
%   flgRaw       (logical) Flag to use raw waveforms instead of pre-processed
%                metrics. If true, waveforms are interpolated and scaled.
%                {false}
%   b2uv         (numeric) Conversion factor from bits to microvolts for
%                raw waveforms. {0.195}
%   UnitType     (categorical) Vector of unit classifications. If empty,
%                loads 'units.type'. {[]}
%
% OUTPUT:
%   wvTbl        (table) The created data table with waveforms.
%   hFig         (handle) Handle to the GUI figure.
%
% DEPENDENCIES:
%   basepaths2vars, v2tbl, tblGUI_xy
%

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addOptional(p, 'basepaths', {}, @(x) iscell(x));
addOptional(p, 'flgRaw', false, @islogical);
addOptional(p, 'b2uv', 0.195, @(x) isnumeric(x) && (isempty(x) || isscalar(x)));
addOptional(p, 'UnitType', [], @(x) iscategorical(x) || isnumeric(x));
addParameter(p, 'Parent', [], @(x) isempty(x) || isgraphics(x));
addParameter(p, 'SelectionCallback', [], @(x) isempty(x) || isa(x, 'function_handle'));

parse(p, varargin{:});
basepaths = p.Results.basepaths;
flgRaw = p.Results.flgRaw;
b2uv = p.Results.b2uv;
UnitType = p.Results.UnitType;
hParent = p.Results.Parent;
selCbk = p.Results.SelectionCallback;

%% ========================================================================
%  LOAD DATA
%  ========================================================================

% Determine variables to load
if flgRaw
    vars = {'swv_raw', 'units'};
else
    vars = {'swv_metrics', 'units'};
end

v = basepaths2vars('basepaths', basepaths, 'vars', vars);

%% ========================================================================
%  PREP METADATA
%  ========================================================================

% Prepare table metadata (Mouse, File, UnitID)
tagFiles = struct();
tagFiles.Mouse = get_mname(basepaths);
[~, fNames] = fileparts(basepaths);
tagFiles.File = fNames;

% Create base table using v2tbl
% We initially just want the structure, metadata and UnitID
% We will manually added the waveform because v2tbl handles scalar/string columns best
% or we can try to map it if we can. But manual is safer for matrix columns.
% Actually, v2tbl might not handle 'swv_raw' which is a cell of matrices well if we map it directly.
% So we'll use v2tbl to get the backbone.

varMap = struct();
varMap.UnitType = 'units.type';
wvTbl = v2tbl('v', v, 'varMap', varMap, 'tagFiles', tagFiles);

% Handle UnitType
if isempty(UnitType)
    units = catfields([v.units], 2, true);
    if isfield(units, 'type') && ~isempty(units.type)
        wvTbl.UnitType = vertcat(units.type{:});
    else
        wvTbl.UnitType = repmat(categorical({'Unknown'}), height(wvTbl), 1);
    end
else
    wvTbl.UnitType = UnitType(:);
end

%% ========================================================================
%  PROCESS WAVEFORMS
%  ========================================================================

if flgRaw
    % Grab raw waveforms
    swv = [v(:).swv_raw];
    % Mean across spikes for each unit if it is a matrix [samples x nSpikes]
    % Usually swv_raw is { [nSamp x nSpikes], ... } per unit?
    % Or cell array of units where each is a mean?
    % Checking previous code: swv = [v(:).swv_raw]; swv = cellfun(@(x) mean(x, 2, 'omitnan'), swv, 'uni', false);
    % So it assumes cell array of matrices.

    swv = cellfun(@(x) mean(x, 2, 'omitnan'), swv, 'uni', false);
else
    % Grab L2 normalized waveforms from metrics
    swv = catfields([v(:).swv], 'cell');
    swv = cellfun(@(x) num2cell(x, 2), swv.wv, 'UniformOutput', false);
    swv = vertcat(swv{:});

end

% Ensure column vector
swv = swv(:);

% Interpolate to 32 points
swv = cellfun(@(x) interp1(linspace(0, 1, length(x)), x, linspace(0, 1, 32)', 'spline'),...
    swv, 'UniformOutput', false);

% Convert to Matrix [nUnits x nSamples]
wv = cell2padmat(swv, 2)';

% Apply Scaling if Raw
if flgRaw && ~isempty(b2uv)
    wv = wv * b2uv;
end

% Assign to Table
wvTbl.Waveform = wv;

%% ========================================================================
%  PLOT
%  ========================================================================

% Time Axis (ms). Assumes 20kHz sampling
xVal = linspace(-0.75, 0.8, 32);

% Launch GUI
hFig = tblGUI_xy(xVal, wvTbl, 'yVar', 'Waveform', 'Parent', hParent, ...
    'SelectionCallback', selCbk);

end


% EOF
