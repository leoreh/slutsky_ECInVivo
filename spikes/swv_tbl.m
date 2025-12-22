function [tblWv, xVal, hFig] = swv_tbl(varargin)

% SWV_TBL Loads waveform data and optionally plots it.
%
% INPUT (Optional Key-Value Pairs):
%   basepaths    (cell array) Full paths to recording folders.
%   flgPlot      (logical) Flag to display GUI. {true}
%   flgRaw       (logical) Use raw waveforms. {false}
%   b2uv         (numeric) Conversion factor. {0.195}
%
% OUTPUT:
%   tblWv        (table) Waveform data.
%   xVal         (numeric) Time axis.
%   hFig         (handle) Figure handle.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addOptional(p, 'basepaths', {}, @(x) iscell(x));
addOptional(p, 'flgPlot', true, @islogical);
addOptional(p, 'flgRaw', false, @islogical);
addOptional(p, 'b2uv', 0.195, @isnumeric);

parse(p, varargin{:});
basepaths = p.Results.basepaths;
flgPlot = p.Results.flgPlot;
flgRaw = p.Results.flgRaw;
b2uv = p.Results.b2uv;

hFig = [];

%% ========================================================================
%  LOAD DATA
%  ========================================================================

if flgRaw
    vars = {'swv_raw', 'units'};
else
    vars = {'swv_metrics', 'units'};
end

v = basepaths2vars('basepaths', basepaths, 'vars', vars);

%% ========================================================================
%  PREP TABLE
%  ========================================================================

% Metadata
tagFiles = struct();
tagFiles.Mouse = get_mname(basepaths);
[~, fNames] = fileparts(basepaths);
tagFiles.File = fNames;

% Base table
varMap = struct();
varMap.UnitType = 'units.type';
tblWv = v2tbl('v', v, 'varMap', varMap, 'tagFiles', tagFiles);

% Handle missing UnitType
if ~ismember('UnitType', tblWv.Properties.VariableNames)
    tblWv.UnitType = repmat(categorical({'Unknown'}), height(tblWv), 1);
end


%% ========================================================================
%  PROCESS WAVEFORMS
%  ========================================================================

if flgRaw
    swv = [v(:).swv_raw];
    % Mean across spikes [nSamp x nSpikes]
    swv = cellfun(@(x) mean(x, 2, 'omitnan'), swv, 'uni', false);
else
    swv = catfields([v(:).swv], 'cell');
    swv = cellfun(@(x) num2cell(x, 2), swv.wv, 'UniformOutput', false);
    swv = vertcat(swv{:});
end

% Flatten
swv = swv(:);

% Interpolate to 32 points
swv = cellfun(@(x) interp1(linspace(0, 1, length(x)), x, linspace(0, 1, 32)', 'spline'),...
    swv, 'UniformOutput', false);

% Convert to Matrix [nUnits x nSamples]
wv = cell2padmat(swv, 2)';

% Scaling
if flgRaw && ~isempty(b2uv)
    wv = wv * b2uv;
end

% Add to table
tblWv.Waveform = wv;

% Time Axis
xVal = linspace(-0.75, 0.8, 32);

%% ========================================================================
%  PLOT
%  ========================================================================

if flgPlot
    hFig = tblGUI_xy(xVal, tblWv, 'yVar', 'Waveform');
end

end     % EOF
