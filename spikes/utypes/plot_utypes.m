function hAx = plot_utypes(varargin)

% PLOT_UTYPES ensures clear visualization of neuronal unit classifications
% (Regular Spiking-RS vs Fast Spiking-FS).
% It supports two types of plots:
%   1. Waveform Plots: Displays average spike waveforms for differentiated cell
%      types with optional scaling and interpolation.
%   2. Scatter Plots: Illustrates the separation between cell types using
%      various metrics (e.g., tp, lidor).
%
% METHODOLOGY:
% Waveforms in 'wv' mode are averaged and interpolated to 32 samples.
% Optional scaling (b2uv) converts bits to microvolts.
% Scatter plots visualize the clustering of RS vs FS units.
%
% INPUT (Optional Key-Value Pairs):
%   basepaths    (cell array) Full paths to recording folders. If empty,
%                uses current directory.
%   flgRaw       (logical) Flag to use raw waveforms instead of pre-processed
%                metrics. If true, waveforms are averaged and interpolated.
%                {true}
%   flgOther     (logical) Flag to include "Other" (unclassified) units in plots.
%                {false}
%   plotType     (char) Type of plot to generate:
%                'wv' - Waveform plot
%                'scatter' - Scatter plot of cell separation
%                'scatter3' - 3D scatter plot using tp, lidor, and mfr
%                {'wv'}
%   clr          (numeric matrix) Nx3 matrix of RGB colors for each cell type.
%                If empty, uses default colors. {[]}
%   b2uv         (numeric) Conversion factor from bits to microvolts for
%                waveform scaling. If empty, no scaling is applied. {[]}
%   fetTbl       (table) Table containing features for plotting. Required for
%                scatter/scatter3 plots.
%   xFet         (char) Field from fetTbl for x-axis. {'tp'}
%   yFet         (char) Field from fetTbl for y-axis. {'lidor'}
%   zFet         (char) Field from fetTbl for z-axis. {'asym'}
%   unitType     (numeric) Vector of unit classifications (1 x nUnits):
%                0 - unclassified
%                1 - Regular Spiking (RS)
%                2 - Fast Spiking (FS)
%                If empty, will use units.clean from basepaths2vars. {[]}
%   hAx          (axes handle) Optional axes handle to plot into. If empty,
%                creates new figure and axes. {[]}
%
% OUTPUT:
%   hAx          (axes handle) Handle to the axes containing the plot.
%
% DEPENDENCIES:
%   basepaths2vars, catfields, cell2padmat, plot_axSize, plot_stdShade
%
% HISTORY:
%   Aug 2024 (AI Assisted) - Added unitType input option, modified to use
%                            unitType instead of units.clean
%   Aug 2024 (AI Assisted) - Added hAx input/output option

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments using inputParser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'fetTbl', table(), @istable);
addOptional(p, 'basepaths', {}, @(x) iscell(x));
addOptional(p, 'flgRaw', true, @islogical);
addOptional(p, 'flgOther', false, @islogical);
addOptional(p, 'plotType', 'wv', @(x) ismember(lower(x), {'wv', 'scatter', 'scatter3'}));
addOptional(p, 'clr', [], @(x) isnumeric(x) && (isempty(x) || size(x, 2) == 3));
addOptional(p, 'b2uv', 0.195, @(x) isnumeric(x) && (isempty(x) || isscalar(x)));
addOptional(p, 'xFet', 'tp', @ischar);
addOptional(p, 'yFet', 'lidor', @ischar);
addOptional(p, 'zFet', 'asym', @ischar);
addOptional(p, 'unitType', [], @(x) isnumeric(x) && (isempty(x) || isvector(x)));
addOptional(p, 'hAx', [], @(x) isempty(x) || ishandle(x) && strcmp(get(x, 'Type'), 'axes'));

parse(p, varargin{:});
fetTbl = p.Results.fetTbl;
basepaths = p.Results.basepaths;
flgRaw = p.Results.flgRaw;
flgOther = p.Results.flgOther;
plotType = lower(p.Results.plotType);
clr = p.Results.clr;
b2uv = p.Results.b2uv;
xFet = p.Results.xFet;
yFet = p.Results.yFet;
zFet = p.Results.zFet;
unitType = p.Results.unitType;
hAx = p.Results.hAx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data if necessary
v = [];
if strcmp(plotType, 'wv')

    if flgRaw
        vars = {'swv_raw'};
    else
        vars = {'swv_metrics'};
    end

    if isempty(unitType)
        vars = [vars, {'units'}];
    end

    v = basepaths2vars('basepaths', basepaths, 'vars', vars);
else

    % Prepare data for scatter
    szVal  = normalize(10 .^ fetTbl.mfr, "range", [10 180])';
    xVal = fetTbl.(xFet);
    yVal = fetTbl.(yFet);
    zVal = fetTbl.(zFet);
end

% get unit classifications
if isempty(unitType)
    units = catfields([v.units], 2);
    unitType = zeros(1, size(units.clean, 2));
    unitType(units.clean(1, :)) = 1;  % RS
    unitType(units.clean(2, :)) = 2;  % FS
end
nUnits = length(unitType);

% Map 'Other' (0) to 3 if requested
if flgOther
    unitType(unitType == 0) = 3;
end

% Default colors if not provided
if isempty(clr)
    clr = mcu_clr();
    clr = clr.unitType;
end

% Update legend text with unit count
txtUnit = {'RS', 'FS'};
if flgOther
    txtUnit{3} = 'Other';
end

for iUnit = 1 : length(txtUnit)
    nUnitsType = sum(unitType == iUnit);
    txtUnit{iUnit} = sprintf('%s (n=%d)', txtUnit{iUnit}, nUnitsType);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot based on type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create or use figure/axes
if isempty(hAx)
    [~, hAx] = plot_axSize('szOnly', false);
end
hold on

switch plotType
    case 'wv'
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
            wv(wv32, :) = wv(wv32, :) * 0.195;
        end

        % waveform params
        xVal = [-15 : 16] / 20000 * 1000;

        % Plot waveforms
        for iUnit = 1 : length(txtUnit)
            % grab average per unit type
            wvUnit = wv(unitType == iUnit, :);
            if isempty(wvUnit)
                continue
            end

            % plot with std shade
            hPlt(iUnit) = plot_stdShade('dataMat', wvUnit, 'xVal', xVal, ...
                'hAx', hAx, 'clr', clr(iUnit, :), 'alpha', 0.5);
        end

        % Format axes
        xlim([xVal(1), xVal(end)])
        xlabel('Time (ms)')
        ylabel('Amplitude (a.u.)')
        legend(hAx, txtUnit(~cellfun(@isempty, txtUnit)), 'Location', 'southeast')

    case 'scatter'
        for iUnit = 1 : length(txtUnit)
            unitIdx = unitType == iUnit;
            scatter(xVal(unitIdx), yVal(unitIdx),...
                szVal(unitIdx), clr(iUnit, :),...
                'filled', 'MarkerFaceAlpha', 0.5)
        end

        % Format axes
        xlabel(xFet, 'Interpreter', 'none')
        ylabel(yFet, 'Interpreter', 'none')
        legend(hAx, txtUnit, 'Location', 'best')

    case 'scatter3'
        for iUnit = 1 : length(txtUnit)
            unitIdx = unitType == iUnit;
            scatter3(xVal(unitIdx),...
                yVal(unitIdx),...
                zVal(unitIdx),...
                szVal(unitIdx),...
                clr(iUnit, :), 'filled', 'MarkerFaceAlpha', 0.3)
        end
        view(3)  % Set 3D view
        camproj(hAx, 'perspective');

        % Format axes
        xlabel(xFet, 'Interpreter', 'none')
        ylabel(yFet, 'Interpreter', 'none')
        zlabel(zFet, 'Interpreter', 'none')
        legend(hAx, txtUnit, 'Location', 'best')
        grid on

end

% % Assert Size
% plot_axSize('hFig', gcf, 'szOnly', false, 'axShape', 'square');

end

% EOF
