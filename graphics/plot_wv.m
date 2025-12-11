function hAx = plot_wv(varargin)

% PLOT_WV displays average spike waveforms for differentiated cell types.
% It visualizes average spike waveforms for Regular Spiking (RS) vs
% Fast Spiking (FS) units with optional scaling and interpolation.
%
% INPUT (Optional Key-Value Pairs):
%   basepaths    (cell array) Full paths to recording folders. If empty,
%                uses current directory.
%   flgRaw       (logical) Flag to use raw waveforms instead of pre-processed
%                metrics. If true, waveforms are averaged and interpolated.
%                {true}
%   flgOther     (logical) Flag to include "Other" (unclassified) units in plots.
%                {false}
%   clr          (numeric matrix) Nx3 matrix of RGB colors for each cell type.
%                If empty, uses default colors. {[]}
%   b2uv         (numeric) Conversion factor from bits to microvolts for
%                waveform scaling. If empty, no scaling is applied. {0.195}
%   UnitType     (categorical or numeric) Vector of unit classifications (1 x nUnits).
%                If categorical: Uses categories for legend.
%                If numeric: 1=RS, 2=FS.
%                If empty: Attempts to load 'units.type' from data. {[]}
%   hAx          (axes handle) Optional axes handle to plot into. If empty,
%                creates new figure and axes. {[]}
%
% OUTPUT:
%   hAx          (axes handle) Handle to the axes containing the plot.
%
% DEPENDENCIES:
%   basepaths2vars, catfields, cell2padmat, plot_axSize, plot_stdShade, mcu_cfg

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addOptional(p, 'basepaths', {}, @(x) iscell(x));
addOptional(p, 'flgRaw', true, @islogical);
addOptional(p, 'flgOther', false, @islogical);
addOptional(p, 'clr', [], @(x) isnumeric(x) && (isempty(x) || size(x, 2) == 3));
addOptional(p, 'b2uv', 0.195, @(x) isnumeric(x) && (isempty(x) || isscalar(x)));
addOptional(p, 'UnitType', [], @(x) iscategorical(x) || isnumeric(x));
addOptional(p, 'hAx', [], @(x) isempty(x) || ishandle(x) && strcmp(get(x, 'Type'), 'axes'));

parse(p, varargin{:});
basepaths = p.Results.basepaths;
flgRaw = p.Results.flgRaw;
flgOther = p.Results.flgOther;
clr = p.Results.clr;
b2uv = p.Results.b2uv;
UnitType = p.Results.UnitType;
hAx = p.Results.hAx;

%% ========================================================================
%  PREPARATIONS
%  ========================================================================

% Load data
if flgRaw
    vars = {'swv_raw'};
else
    vars = {'swv_metrics'};
end

if isempty(UnitType)
    vars = [vars, {'units'}];
end

v = basepaths2vars('basepaths', basepaths, 'vars', vars);

% Get unit classifications if not provided
if isempty(UnitType)
    units = catfields([v.units], 2);

    if isfield(units, 'type')
        UnitType = units.type;
    elseif isfield(units, 'clean')          % Legacy
        UnitType = zeros(1, size(units.clean, 2));
        UnitType(units.clean(1, :)) = 1;    % RS
        UnitType(units.clean(2, :)) = 2;    % FS
    end
end

% Handle categorical vs numeric configuration
if iscategorical(UnitType)
    % Extract categories for legend
    cats = categories(UnitType);

    % If flgOther is false, remove 'Other' category if it exists
    if ~flgOther
        idxRemove = strcmpi(cats, 'Other');
    end

    % Define unit groups based on categories
    txtUnit = cats;
    nTypes = length(cats);

    % Convert UnitType to indices for the loop
    [unitIdx, unitCats] = grp2idx(UnitType);
else
    % Numeric handling (legacy 1=RS, 2=FS)
    txtUnit = {'RS', 'FS'};
    if flgOther
        UnitType(UnitType == 0) = 3;
        txtUnit{3} = 'Other';
    else
        % Ensure we don't plot 0 if flgOther is false
        UnitType(UnitType == 0) = NaN;
    end
    unitIdx = UnitType;
    nTypes = length(txtUnit);
end

% Default colors if not provided
if isempty(clr)
    cfg = mcu_cfg();
    clr = cfg.clr.unit;
end

% Update legend text with unit count and prepare iteration list
plotParams = struct('idx', {}, 'lbl', {}, 'clr', {});
count = 0;

for iType = 1 : nTypes

    % Check if we should skip 'Other'
    if ~flgOther && (strcmpi(txtUnit{iType}, 'Other') || iType == 3)
        continue;
    end

    % Logic mask for this group
    idxMsk = (unitIdx == iType);

    % Map to color (linear mapping)
    if iType <= size(clr, 1)
        currClr = clr(iType, :);
    else
        currClr = [0 0 0];
    end

    nUnitsType = sum(idxMsk);
    if nUnitsType == 0
        continue;
    end

    count = count + 1;
    plotParams(count).idx = idxMsk; % Logic mask for this group
    plotParams(count).lbl = sprintf('%s (n=%d)', txtUnit{iType}, nUnitsType);
    plotParams(count).clr = currClr;
end


%% ========================================================================
%  PLOT
%  ========================================================================

% Create or use figure/axes
if isempty(hAx)
    [~, hAx] = plot_axSize('szOnly', false);
end
hold on

% Get waveforms
nfiles = length(v);
swv = cell(nfiles, 1);

if flgRaw
    % grab raw waveforms
    swv = [v(:).swv_raw];
    swv = cellfun(@(x) mean(x, 2, 'omitnan'), swv, 'uni', false);

    % Fig params
    txtY = 'Amplitude (µV)';
else
    % Process each file separately
    for ifile = 1 : nfiles
        % convert each unit's waveform to a cell
        nUnits = size(v(ifile).swv.wv, 1);
        swv{ifile} = mat2cell(v(ifile).swv.wv, ones(nUnits, 1), size(v(ifile).swv.wv, 2));
    end
    % concatenate all cells
    swv = cat(1, swv{:});

    % Fig params
    txtY = 'Norm. Amplitude (a.u.)';
end

% store waveform length and interpolate all to 32 samples
wvLen = cellfun(@length, swv, 'uni', true);
wv32 = wvLen == 32;

% Interpolate all waveforms to 32 samples
swv = cellfun(@(x) interp1(linspace(0, 1, length(x)), x, linspace(0, 1, 32)', 'spline'),...
    swv, 'UniformOutput', false);

% concate
wv = cell2padmat(swv, 2)';

% apply b2uv scaling
if ~isempty(b2uv) && flgRaw
    wv(wv32, :) = wv(wv32, :) * b2uv;
end

% waveform params
xVal = [-15 : 16] / 20000 * 1000;

% Plot waveforms
plottedLbls = {};
for iGrp = 1 : length(plotParams)

    % grab average per unit type
    wvUnit = wv(plotParams(iGrp).idx, :);
    if isempty(wvUnit)
        continue
    end

    % plot with std shade
    plot_stdShade('dataMat', wvUnit, 'xVal', xVal, ...
        'hAx', hAx, 'clr', plotParams(iGrp).clr, 'alpha', 0.5);

    plottedLbls{end+1} = plotParams(iGrp).lbl;
end

% Format axes
xlim([xVal(1), xVal(end)])
xlabel('Time (ms)')
ylabel(txtY)
legend(hAx, plottedLbls, 'Location', 'southeast')

end

% EOF
