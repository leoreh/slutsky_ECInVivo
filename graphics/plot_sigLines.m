function plot_sigLines(axh, barIdx, barLbl, varargin)
% PLOT_SIGLINES Adds significance notation (lines and labelsbarIdx) between specified bars.
%
% INPUTS:
%   axh         Axes handle containing the bar plot.
%   barIdx      Cell array. Each cell contains a 2-element vector
%               specifying the 1-based indices of bars to connect
%               (from left to right). E.g., {[1,2], [3,5]}
%   barLbl      Cell array of strings, corresponding to labels for barIdx.
%               E.g., {'*', 'p < .01'}
%   VARARGIN:
%     'LineColor'       (default 'k')
%     'LineWidth'       (default 1)
%     'lineOffset'      (default 0.06 of Y-axis range, from bar top to line)
%     'txtOffset'       (default 0.005 of Y-axis range, from line to text)
%     'fontSize'        (default get(axh, 'FontSize'))
%     'lineGap'         (default 0.03 of Y-axis range, between stacked notations)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING & INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section defines and parses input arguments using MATLAB's inputParser.
% Default values are assigned, and initial parameters are set.
p = inputParser;
addRequired(p, 'axh', @(x) isa(x, 'matlab.graphics.axis.Axes'));
addRequired(p, 'barIdx', @iscell);
addRequired(p, 'barLbl', @iscell);

addParameter(p, 'LineColor', 'k');
addParameter(p, 'LineWidth', 1);
addParameter(p, 'lineOffset', 0.1);
addParameter(p, 'txtOffset', -0.01);
addParameter(p, 'fontSize', get(axh, 'FontSize'));
addParameter(p, 'lineGap', 0.06);

parse(p, axh, barIdx, barLbl, varargin{:});
opts = p.Results;

if length(barIdx) ~= length(barLbl)
    error('barIdx and barLbl must have the same number of elements.');
end

hold(axh, 'on');
yLimit = get(axh, 'YLim');
yRange = diff(yLimit);

% Convert fractional offsets to data units
opts.lineOffset = opts.lineOffset * yRange;
opts.txtOffset = opts.txtOffset * yRange;
opts.lineGap = opts.lineGap * yRange;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BAR INFORMATION EXTRACTION & PREPARATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section retrieves information about the bars in the plot,
% including their x-positions, y-heights, and widths.
% Bar details are then sorted by their x-coordinates to match visual order.
hBars = findobj(axh, 'Type', 'Bar');
if isempty(hBars)
    warning('No bar objects found in the axes.');
    return;
end

% Consolidate bar data and sort by x-position. If multiple series from one
% 'bar' call, findobj returns them in reverse order of creation. We want
% them in series order (1, 2, ... N) if 'SeriesIndex' is available.
barDetails = [];
if length(hBars) > 1 && isprop(hBars(1), 'SeriesIndex')
    try
        seriesIndices = [hBars.SeriesIndex];
        [~, sortOrder] = sort(seriesIndices);
        hBars = hBars(sortOrder);
    catch
        hBars = flipud(hBars); % Fallback if SeriesIndex fails or not present
    end
elseif length(hBars) > 1 % If plotted individually, flipud is a common need for findobj order
    hBars = flipud(hBars);
end

for iBar = 1:length(hBars)
    hb = hBars(iBar);
    xDataSeries = hb.XData; % Group centers for this series
    yDataSeries = hb.YData;
    xOffset = hb.XOffset;   % Offset of this series within group

    for iX = 1:length(xDataSeries)
        barXCenter = xDataSeries(iX) + xOffset;
        barYHeight = yDataSeries(iX);
        barWidthVal = NaN; % Initialize barWidthVal

        % Get bar width (more robust with XEndPoints if available)
        % Fallback: Estimate width (less accurate)
        if isprop(hb, 'BarWidth') && ~isempty(hb.BarWidth)
            currentBarWidth = hb.BarWidth;
            nominalWidthFraction = currentBarWidth;

            uniqueFiniteX = unique(xDataSeries(isfinite(xDataSeries)));
            minGroupSpacing = min(diff(sort(uniqueFiniteX)));
            barWidthVal = nominalWidthFraction * minGroupSpacing / length(hBars);
        end

        barDetails = [barDetails; struct('x', barXCenter, 'y', barYHeight, 'w', barWidthVal)];
    end
end

% Sort bars by x-coordinate to match visual order
[~, sortIdx] = sort([barDetails.x]);
barDetails = barDetails(sortIdx);
nBars = length(barDetails);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOOP THROUGH SIGNIFICANCE PAIRS & DRAWING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section iterates through each pair of bars/midpoints specified for
% significance notation. For each pair, it calculates coordinates,
% determines the vertical position considering stacking, and draws the
% horizontal line and text label.
yTop = 0;
for iBar = 1:length(barIdx)

    % --- Calculate x and y for the first point ---
    if mod(barIdx{iBar}(1), 1) ~= 0 % Non-integer, implies midpoint
        barLow = barDetails(floor(barIdx{iBar}(1)));
        barHigh = barDetails(ceil(barIdx{iBar}(1)));
        xCoords(1) = (barLow.x + barHigh.x) / 2;
        yCoords(1) = max(barLow.y, barHigh.y);
    else % Integer index
        xCoords(1) = barDetails(barIdx{iBar}(1)).x;
        yCoords(1) = barDetails(barIdx{iBar}(1)).y;
    end

    % --- Calculate x and y for the second point ---
    if mod(barIdx{iBar}(2), 1) ~= 0 % Non-integer, implies midpoint
        barLow = barDetails(floor(barIdx{iBar}(2)));
        barHigh = barDetails(ceil(barIdx{iBar}(2)));
        xCoords(2) = (barLow.x + barHigh.x) / 2;
        yCoords(2) = max(barLow.y, barHigh.y);
    else % Integer index
        xCoords(2) = barDetails(barIdx{iBar}(2)).x;
        yCoords(2) = barDetails(barIdx{iBar}(2)).y;
    end

    % determine height of line
    if iBar == 1
        yHeight = max([yCoords, yTop]) + opts.lineOffset;
    else
        yHeight = max([yCoords, yTop]) + opts.lineGap;
    end

    % Draw horizontal line
    plot(axh, xCoords, [yHeight, yHeight], ...
        'Color', opts.LineColor,...
        'LineWidth', opts.LineWidth,...
        'HandleVisibility', 'off');

    % Add text label
    textX = mean(xCoords);
    textY = yHeight + opts.txtOffset;
    th = text(axh, textX, textY, barLbl{iBar}, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', opts.fontSize, ...
        'Color', opts.LineColor, ...
        'HandleVisibility', 'off');

    % Update maxYOverall and activeLevels for stacking
    drawnow; % Ensure text extent is calculated
    txtExtent = get(th, 'Extent');
    yTop = txtExtent(2) + txtExtent(4); % yPos + yHeight
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADJUST Y-LIMITS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section adjusts the Y-axis limits of the plot if the newly added
% significance notations extend beyond the current limits, ensuring all
% elements are visible.
yPad = 0;
if yTop > yLimit(2)
    set(axh, 'YLim', [yLimit(1), yTop + yPad * yRange]); % Add 5% padding
end

hold(axh, 'off');
end