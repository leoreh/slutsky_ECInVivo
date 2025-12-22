function [hFig, hAx] = plot_axSize(varargin)
% plot_defaults - Opens/configures a figure and axes with specified plot area dimensions.
%
% Figure dimensions are determined by the requested axis plot area dimensions,
% taking into account space needed for titles and labels.
%
% INPUT (Optional name-value pairs):
%   'hFig'        - Existing figure handle. If not provided, a new figure is created.
%   'axShape'     - Shape of the axes plot area: {'square', 'tall', 'wide'}.
%                   Default is 'square'. Used if only one of 'axWidth' or 'axHeight'
%                   is provided, or if neither is provided (to set default plot area shape).
%   'axWidth'     - Desired width of the axes plot area in pixels.
%   'axHeight'    - Desired height of the axes plot area in pixels.
%   'szOnly'      - Only adjust figure size without default aesthetics
%   'flgPos'      - If true, moves the figure to monitor 2. Default is false.
%
% OUTPUT:
%   hFig          - Handle to the figure.
%   hAx           - Handle to the axes.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse input arguments
p = inputParser;
addParameter(p, 'hFig', [], @(x) isa(x, 'matlab.ui.Figure') || isempty(x));
addParameter(p, 'axShape', 'square', @(x) ischar(x) && ismember(x, {'square', 'tall', 'wide'}));
addParameter(p, 'axWidth', [], @(x) isnumeric(x));
addParameter(p, 'axHeight', [], @(x) isnumeric(x));
addParameter(p, 'szOnly', true, @islogical);
addParameter(p, 'flgPos', false, @islogical);
addParameter(p, 'flgFullscreen', false, @islogical);

parse(p, varargin{:});
hFig            = p.Results.hFig;
axShape         = p.Results.axShape;
axWidth         = p.Results.axWidth;
axHeight        = p.Results.axHeight;
szOnly          = p.Results.szOnly;
flgPos          = p.Results.flgPos;
flgFullscreen   = p.Results.flgFullscreen;

% Define aspect ratios and default dimension
tallRatio = 0.62; % width/height for tall shape
wideRatio = 1.62; % width/height for wide shape
defDim = 300;     % Default dimension for plot area calculation

% Initialize figure if not provided
if isempty(hFig)
    hFig = figure;
else
    hFig = hFig;
    figure(hFig); % Bring to front if existing
end

% Create or get axes
hAx = findobj(hFig, 'Type', 'Axes');
if isempty(hAx)
    hAx = gca; % Get current axes or create one if none exists
else
    hAx = hAx(1); % Use the first axes if multiple exist
    axes(hAx);    % Make it the current axes
end
axUnits = get(hAx, 'Units');
set(hAx, 'Units', 'pixels');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT AESTHETICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default properties early, as they affect TightInset
if ~szOnly
    set(hFig, 'Color', 'w');

    fntSize = 16;
    FntName = 'Arial';
    set(hAx, 'FontName', FntName, 'FontSize', fntSize);
    set(hAx, 'XTickLabelRotation', 0);

    box(hAx, 'OFF');
    hold(hAx, 'on');

    % Set font and font size for the title, axis labels, colorbar and legend
    hTtl = get(hAx, 'Title');
    set(hTtl, 'FontName', FntName, 'FontSize', fntSize + 4);

    hLbl = get(hAx, 'XLabel');
    set(hLbl, 'FontName', 'Arial', 'FontSize', fntSize + 4);
    hLbl = get(hAx, 'YLabel');
    set(hLbl, 'FontName', 'Arial', 'FontSize', fntSize + 4);

    hLgnd = get(hAx, 'Legend');
    set(hLgnd, 'FontName', 'Arial', 'FontSize', fntSize);

    hCb = findobj(hFig, 'Type', 'Colorbar');
    set(hCb, 'FontName', 'Arial', 'FontSize', fntSize);
    hCbLbl = get(hCb, 'Label');
    set(hCbLbl, 'FontName', 'Arial', 'FontSize', fntSize + 4);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETERMINE AXIS SIZE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For tall shape, width = height * tallRatio => height = width / tallRatio
% For wide shape, width = height * wideRatio => height = width / wideRatio
if ~isempty(axWidth) && ~isempty(axHeight)
    axWidth = axWidth;
    axHeight = axHeight;
elseif ~isempty(axWidth) && isempty(axHeight)
    switch axShape
        case 'square'
            axHeight = axWidth;
        case 'tall'
            axHeight = axWidth / tallRatio;
        case 'wide'
            axHeight = axWidth / wideRatio;
    end
elseif ~isempty(axHeight) && isempty(axWidth)
    switch axShape
        case 'square'
            axWidth = axHeight;
        case 'tall'
            axWidth = axHeight * tallRatio;
        case 'wide'
            axWidth = axHeight * wideRatio;
    end
else
    % Neither axWidth nor axHeight provided explicitly, use axShape for
    % default plot area
    switch axShape
        case 'square'
            axWidth = defDim;
            axHeight = defDim;
        case 'tall' % Default height is default_plot_dim, width is proportional
            axHeight = defDim;
            axWidth = defDim * tallRatio;
        case 'wide' % Default width is default_plot_dim, height is proportional
            axWidth = defDim;
            axHeight = defDim / wideRatio;
    end
end

% Check for colorbar and get its width if it exists
hCb = findobj(hFig, 'Type', 'Colorbar');
cbWidth = 0;
if ~isempty(hCb)
    cbUnits = get(hCb, 'Units');
    set(hCb, 'Units', 'pixels');
    cbPos = get(hCb, 'Position');

    % Get colorbar label to account for its width
    hCbLbl = get(hCb, 'Label');
    if ~isempty(hCbLbl)
        cbLblUnits = get(hCbLbl, 'Units');
        set(hCbLbl, 'Units', 'pixels');
        cbLblExtent = get(hCbLbl, 'Extent');
        set(hCbLbl, 'Units', cbLblUnits);
        % Add label width plus padding for the label
        cbWidth = cbPos(3) + cbLblExtent(3) + 15; % Increased padding for label
    else
        cbWidth = cbPos(3) + 5; % Just add small padding if no label
    end
    set(hCb, 'Units', cbUnits);
end
axWidth = axWidth - cbWidth;

% Set axes Position to the target plot area dimensions.
% The x,y coordinates [1, 1] are temporary and assume the current figure is large enough
% not to clip content that affects TightInset.
set(hAx, 'Position', [1, 1, axWidth, axHeight]);

% Force MATLAB to render/update the layout based on current properties.
drawnow;

% Get TightInset: [left, bottom, right, top] margins in pixels
ti = get(hAx, 'TightInset');

% Calculate figure dimensions needed to encompass the plot area and its decorations
figWidth = axWidth + ti(1) + ti(3) + cbWidth + 15; % Add colorbar width to total figure width
% figWidth = axWidth + ti(1) + ti(3) + cbWidth; % Add colorbar width to total figure width
figHeight = axHeight + ti(2) + ti(4);

% Set figure size. Preserve original screen position (left, bottom) of the figure.
fhUnits = get(hFig, 'Units');
set(hFig, 'Units', 'pixels');
figPos = get(hFig, 'Position'); % To preserve screen x,y for the figure
set(hFig, 'Position', [figPos(1), figPos(2), figWidth, figHeight]);

% Position figure on 2nd screen
if flgPos
    monitors = get(0, 'MonitorPositions');  % Get positions of all monitors
    % Check if 2nd monitor exists
    if size(monitors, 1) > 1
        mnt = monitors(2, :);      % Select second screen [left, bottom, width, height]
    else
        mnt = monitors(1, :);      % Fallback to first screen
    end

    % flgPos case: Offset position on target monitor
    hFig.Position = [mnt(1) + 400, mnt(2) + 200,...
        figWidth, figHeight];
end

% Force MATLAB to render/update the figure with its new size.
drawnow;

% Position the axes within the now-resized figure
set(hAx, 'Position', [ti(1), ti(2), axWidth, axHeight]);

% If colorbar exists, position it relative to the axis
if ~isempty(hCb)
    set(hCb, 'Units', 'pixels');
    % Position colorbar to the right of the axis with minimal gap
    % Use the axis position as reference
    axPos = get(hAx, 'Position');
    set(hCb, 'Position', [axPos(1) + axPos(3) + 5, axPos(2), cbPos(3), axPos(4)]);
    set(hCb, 'Units', cbUnits);
end

% Restore original units for figure and axes
set(hFig, 'Units', fhUnits);
set(hAx, 'Units', axUnits);

% Set to cover the entire target monitor
if flgFullscreen
    hFig.WindowState = 'maximized';
end

end     % EOF