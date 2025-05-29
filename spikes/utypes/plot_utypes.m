function hAx = plot_utypes(varargin)
% PLOT_UTYPES Generates plots for cell classification based on waveform and spike metrics.
%
% SUMMARY:
% This function creates visualization plots for neuronal cell classification
% based on waveform characteristics and spike train metrics. It supports two
% types of plots:
%   1. Waveform Plots: Displays average spike waveforms for different cell
%      types (e.g., putative pyramidal cells vs. interneurons), with optional
%      scaling and interpolation to standardize waveform length.
%   2. Scatter Plots: Illustrates the separation between cell types using
%      various metrics (e.g., trough-to-peak duration vs. burstiness index).
%
% METHODOLOGY:
% Waveforms can be provided as raw data or pre-processed metrics. Raw
% waveforms are averaged across spikes and interpolated to a standard length
% (32 samples). Optional scaling (b2uv) can be applied to convert to
% microvolts. For scatter plots, any combination of waveform and spike train
% metrics can be used to visualize cell type separation.
%
% INPUT (Optional Key-Value Pairs):
%   basepaths    (cell array) Full paths to recording folders. If empty,
%                uses current directory.
%   flgRaw       (logical) Flag to use raw waveforms instead of pre-processed
%                metrics. If true, waveforms are averaged and interpolated.
%                {true}
%   plotType     (char) Type of plot to generate:
%                'wv' - Waveform plot
%                'scatter' - Scatter plot of cell separation
%                'scatter3' - 3D scatter plot using tp, lidor, and mfr
%                {'wv'}
%   clr          (numeric matrix) Nx3 matrix of RGB colors for each cell type.
%                If empty, uses default colors. {[]}
%   b2uv         (numeric) Conversion factor from bits to microvolts for
%                waveform scaling. If empty, no scaling is applied. {[]}
%   swvFld       (char) Field from swv_metrics to use for x-axis in scatter
%                plot (e.g., 'tp' for trough-to-peak). {'tp'}
%   stFld        (char) Field from st_metrics to use for y-axis in scatter
%                plot (e.g., 'lidor' for burstiness index). {'lidor'}
%   unitIdx      (numeric) Vector of unit classifications (1 x nUnits):
%                0 - unclassified
%                1 - putative pyramidal cell (pPYR)
%                2 - putative interneuron (pINT)
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
%   Aug 2024 (AI Assisted) - Added unitIdx input option, modified to use
%                            unitIdx instead of units.clean
%   Aug 2024 (AI Assisted) - Added hAx input/output option

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments using inputParser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'basepaths', @(x) iscell(x));
addOptional(p, 'flgRaw', true, @islogical);
addOptional(p, 'plotType', 'wv', @(x) ismember(lower(x), {'wv', 'scatter', 'scatter3'}));
addOptional(p, 'clr', [], @(x) isnumeric(x) && (isempty(x) || size(x, 2) == 3));
addOptional(p, 'b2uv', 0.195, @(x) isnumeric(x) && (isempty(x) || isscalar(x)));
addOptional(p, 'swvFld', 'tp', @ischar);  % field from swv_metrics for x-axis
addOptional(p, 'stFld', 'lidor', @ischar);  % field from st_metrics for y-axis
addOptional(p, 'unitIdx', [], @(x) isnumeric(x) && (isempty(x) || isvector(x)));
addOptional(p, 'hAx', [], @(x) isempty(x) || ishandle(x) && strcmp(get(x, 'Type'), 'axes'));

parse(p, varargin{:});
basepaths = p.Results.basepaths;
flgRaw = p.Results.flgRaw;
plotType = lower(p.Results.plotType);
clr = p.Results.clr;
b2uv = p.Results.b2uv;
swvFld = p.Results.swvFld;
stFld = p.Results.stFld;
unitIdx = p.Results.unitIdx;
hAx = p.Results.hAx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load state vars
vars = {'swv_metrics', 'st_metrics', 'units', 'fr'};
if flgRaw
    vars = [vars, {'swv_raw'}];
end
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

% Get metrics
swv = catfields([v.swv], 2);
st = catfields([v.st], 2);
fr = catfields([v.fr], 1);

% get unit classifications
if isempty(unitIdx) || isscalar(unitIdx)
    if isscalar(unitIdx)
        unitFld = ['clean', num2str(unitIdx)];
    else
        unitFld = clean;
    end
    % concatenate units
    units = catfields([v.units], 2);
    % convert units.clean to unitIdx
    unitIdx = zeros(1, size(units.(unitFld), 2));
    unitIdx(units.(unitFld)(1, :)) = 1;  % pPYR
    unitIdx(units.(unitFld)(2, :)) = 2;  % pINT
end
nUnits = length(unitIdx);

% Default colors if not provided
if isempty(clr)
    clr(1, :) = [10 / 255 10 / 255 80 / 255];       % pPYR
    clr(2, :) = [180 / 255 80 / 255 80 / 255];      % pINT
end
txtUnit = {'pPYR', 'pINT'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot based on type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create or use figure/axes
if isempty(hAx)
    [~, hAx] = plot_axSize('szOnly', false);
end

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
        for iUnit = 1 : 2
            % grab average per unit type
            wvUnit = wv(unitIdx == iUnit, :);
            if isempty(wvUnit)
                continue
            end

            % plot with std shade
            hPlt(iUnit) = plot_stdShade('dataMat', wvUnit, 'xVal', xVal, ...
                'axh', hAx, 'clr', clr(iUnit, :), 'alpha', 0.5);
            
            % Update legend text with unit count
            nUnitsType = sum(unitIdx == iUnit);
            txtUnit{iUnit} = sprintf('%s (n=%d)', txtUnit{iUnit}, nUnitsType);
        end
        
        % Format axes
        xlim([xVal(1), xVal(end)])
        xlabel('Time (ms)')
        ylabel('Amplitude (a.u.)')
        legend(hAx, txtUnit(~cellfun(@isempty, txtUnit)), 'Location', 'southeast')
        
        % Add non-classified units
        wvUnit = wv(unitIdx == 0, :);
        if ~isempty(wvUnit) && 0
            plot_stdShade('dataMat', wvUnit, 'xVal', xVal, ...
                'axh', hAx, 'clr', [1 0 0], 'alpha', 0.5);
        end

    case 'scatter'
        
        % Select data (manually)
        szVal  = normalize(fr.mfr, "range", [10 180])';
        yVal  = st.(stFld);
        xVal  = swv.(swvFld);
        % yVal  = swv.tailSlope;
        
        for iUnit = 1 : 2
            scatter(xVal(unitIdx == iUnit), yVal(unitIdx == iUnit),...
                szVal(unitIdx == iUnit), clr(iUnit, :),...
                'filled', 'MarkerFaceAlpha', 0.2)
        end

        % Format axes
        % xlim([0, max(swv.(swvFld))])
        xlabel([upper(swvFld(1)) swvFld(2:end) ' (ms)'])
        ylabel([upper(stFld(1)) stFld(2:end) ' (a.u.)'])
        legend(hAx, txtUnit, 'Location', 'northwest')

    case 'scatter3'
        
        % Select data (manually)
        szVal  = normalize(fr.mfr, "range", [20 180])';
        zVal  = st.(stFld);
        yVal  = swv.(swvFld);
        xVal  = asinh(swv.tailSlope);
        xVal  = swv.tp;

        % Plot scatter3
        for iUnit = 1 : 2
            scatter3(xVal(unitIdx == iUnit),...
                yVal(unitIdx == iUnit),...
                zVal(unitIdx == iUnit),...
                szVal(unitIdx == iUnit),...
                clr(iUnit, :), 'filled', 'MarkerFaceAlpha', 0.3)
        end
        view(3)  % Set 3D view
        camproj(hAx, 'perspective');

        % Format axes
        % zlim([0, max(swv.(swvFld))])
        % ylim([min(st.(stFld)), max(st.(stFld))])
        % xlim([0, max(fr.mfr)])
        % zlabel(['Waveform'])
        ylabel('Y')
        xlabel('X')
        zlabel('Z')
        % legend(hAx, txtUnit, 'Location', 'northwest')
        grid on
        % set(hAx, 'XScale', 'log')

end

% % Assert Size
% plot_axSize('hFig', gcf, 'szOnly', false, 'axShape', 'square');

end

% EOF



