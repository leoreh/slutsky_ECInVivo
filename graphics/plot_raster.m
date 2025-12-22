function [hAx, hPlt] = plot_raster(spktimes, varargin)
% PLOT_RASTER Creates a raster plot from spike times.
%
%   [hAx, hPlt] = PLOT_RASTER(SPKTIMES, ...) efficiently generates a raster
%   plot for neural data. It supports multiple plot styles and is optimized
%   for speed using NaN-separated line plotting.
%
%   Based on 'plotSpikeRaster' by Jeffrey Chiou (2014).
%
%   INPUTS:
%       spktimes    - (cell) Mx1 cell array of spike times, where M is the
%                     number of units/trials. Each cell contains a vector
%                     of spike times in seconds.
%       varargin    - (param/value) Optional parameters:
%                     'hAx'          : (axes) Target axes {gca}
%                     'plotType'     : (char) 'vertline', 'horzline', 'scatter'
%                     'clr'          : (opt) Plot clr {[0.2 0.2 0.2]}
%                     'lineWidth'    : (num) Line width for 'line' plots {0.5}
%                     'lineHeight'   : (num) Vertical height (0-1) {0.9}
%                     'markerSz'     : (num) Marker size for 'scatter' {5}
%                     'spkDur'       : (num) Duration of spike line (s) {0.001}
%                     'xLim'         : (vec) [min max] x-axis limits
%                     'flgLbls'      : (log) Add labels {true}
%
%   OUTPUTS:
%       hAx         - (axes) Handle to the axes.
%       hPlt        - (handle) Handle to the plot object.
%

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'spktimes', @iscell);
addParameter(p, 'hAx', [], @(x) isempty(x) || isa(x, 'matlab.graphics.axis.Axes'));
addParameter(p, 'plotType', 'vertline', @(x) any(validatestring(x, {'vertline', 'horzline', 'scatter'})));
addParameter(p, 'clr', [0.2 0.2 0.2], @(x) isnumeric(x) || ischar(x));
addParameter(p, 'lineWidth', 0.5, @isnumeric);
addParameter(p, 'markerSz', 5, @isnumeric);
addParameter(p, 'spkDur', 0.001, @isnumeric);
addParameter(p, 'lineHeight', 0.9, @isnumeric);
addParameter(p, 'xLim', [], @isnumeric);
addParameter(p, 'flgLbls', true, @islogical);

parse(p, spktimes, varargin{:});
hAx         = p.Results.hAx;
plotType    = p.Results.plotType;
clr         = p.Results.clr;
lineWidth   = p.Results.lineWidth;
lineHeight  = p.Results.lineHeight;
markerSz    = p.Results.markerSz;
spkDur      = p.Results.spkDur;
xLim        = p.Results.xLim;
flgLbls     = p.Results.flgLbls;


%% ========================================================================
%  PREPARATIONS
%  ========================================================================

% Ensure column vector of cells
spktimes = spktimes(:);
nUnits   = length(spktimes);

% Validate content and force row vectors within cells
% (plotSpikeRaster had a bug here, we fix it by ensuring row vectors)
for iUnit = 1:nUnits
    if ~isempty(spktimes{iUnit})
        spktimes{iUnit} = spktimes{iUnit}(:)';
    end
end

if isempty(hAx)
    hAx = gca;
end

% Set hold on
hold(hAx, 'on');


%% ========================================================================
%  PLOTTING LOGIC
%  ========================================================================

xPoints = [];
yPoints = [];

% Pre-calculate total spikes to pre-allocate if needed,
% but MATLAB handles dynamic growth reasonably well for vectors.
% For maximum speed, we follow the NaN separator approach.

nTotalSpikes = sum(cellfun(@length, spktimes));

if nTotalSpikes == 0
    warning('No spikes to plot.');
    hPlt = [];
    return;
end

switch plotType
    % ---------------------------------------------------------------------
    % Vertical Lines
    % ---------------------------------------------------------------------
    case 'vertline'

        % Pre-allocate (3 points per spike: top, bottom, NaN)
        xPoints = nan(nTotalSpikes * 3, 1);
        yPoints = nan(nTotalSpikes * 3, 1);

        halfH = lineHeight / 2;
        idx   = 1;

        for iUnit = 1:nUnits
            st = spktimes{iUnit};
            nSpks = length(st);

            if nSpks == 0, continue; end

            % X coords: [t, t, NaN]
            x_blk = [st; st; nan(1, nSpks)];

            % Y coords: [unit-h, unit+h, NaN]
            y_blk = [(iUnit - halfH) * ones(1, nSpks);
                (iUnit + halfH) * ones(1, nSpks);
                nan(1, nSpks)];

            % Fill main vector
            nPts = nSpks * 3;
            xPoints(idx : idx + nPts - 1) = x_blk(:);
            yPoints(idx : idx + nPts - 1) = y_blk(:);

            idx = idx + nPts;
        end

        hPlt = plot(hAx, xPoints, yPoints, 'color', clr, 'LineWidth', lineWidth);

        % -----------------------------------------------------------------
        % Horizontal Lines
        % -----------------------------------------------------------------
    case 'horzline'

        xPoints = nan(nTotalSpikes * 3, 1);
        yPoints = nan(nTotalSpikes * 3, 1);
        idx = 1;

        for iUnit = 1:nUnits
            st = spktimes{iUnit};
            nSpks = length(st);

            if nSpks == 0, continue; end

            % X coords: [t, t+dur, NaN]
            x_blk = [st; st + spkDur; nan(1, nSpks)];

            % Y coords: [unit, unit, NaN]
            y_blk = [iUnit * ones(1, nSpks);
                iUnit * ones(1, nSpks);
                nan(1, nSpks)];

            nPts = nSpks * 3;
            xPoints(idx : idx + nPts - 1) = x_blk(:);
            yPoints(idx : idx + nPts - 1) = y_blk(:);

            idx = idx + nPts;
        end

        hPlt = plot(hAx, xPoints, yPoints, 'color', clr, 'LineWidth', lineWidth);

        % -----------------------------------------------------------------
        % Scatter
        % -----------------------------------------------------------------
    case 'scatter'

        % Concatenate all times and unit indices
        % This is vectorized and fast

        % Creates a vector of unit indices matching the spike times
        unitIds = cellfun(@(x, id) id * ones(size(x)), spktimes, num2cell(1:nUnits)', 'UniformOutput', false);

        xPoints = [spktimes{:}];
        yPoints = [unitIds{:}];

        hPlt = plot(hAx, xPoints, yPoints, '.', 'color', clr, 'markerSz', markerSz);

end


%% ========================================================================
%  FORMATTING
%  ========================================================================

set(hAx, 'YDir', 'reverse');
set(hAx, 'TickDir', 'out');

if ~isempty(xLim)
    xlim(hAx, xLim);
else
    % Auto-scale X to fit data tight
    minT = min(xPoints, [], 'omitnan');
    maxT = max(xPoints, [], 'omitnan');
    if ~isempty(minT) && ~isempty(maxT)
        xlim(hAx, [minT, maxT]);
    end
end

ylim(hAx, [0, nUnits + 1]);

if flgLbls
    xlabel(hAx, 'Time (s)');
    ylabel(hAx, 'Unit');
end

end     % EOF
