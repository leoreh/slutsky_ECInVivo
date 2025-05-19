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

parse(p, varargin{:});
hFig            = p.Results.hFig;
axShape         = p.Results.axShape;
axWidth         = p.Results.axWidth;
axHeight        = p.Results.axHeight;
szOnly          = p.Results.szOnly;

% Define aspect ratios and default dimension
tallRatio = 0.62; % width/height for tall shape
wideRatio = 1.62; % width/height for wide shape
defDim = 400;     % Default dimension for plot area calculation

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

    % Set font and font size for the title, axis labels, and legend
    hTtl = get(hAx, 'Title');
    set(hTtl, 'FontName', FntName, 'FontSize', fntSize + 4);

    hLbl = get(hAx, 'XLabel');
    set(hLbl, 'FontName', 'Arial', 'FontSize', fntSize + 4);
    hLbl = get(hAx, 'YLabel');
    set(hLbl, 'FontName', 'Arial', 'FontSize', fntSize + 4);

    hLgnd = get(hAx, 'Legend');
    set(hLgnd, 'FontName', 'Arial', 'FontSize', fntSize);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADJUST FIGURE SIZE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set axes Position to the target plot area dimensions.
% The x,y coordinates [1, 1] are temporary and assume the current figure is large enough
% not to clip content that affects TightInset. TightInset calculation depends primarily
% on the Width and Height of the Position property, not its absolute X,Y location,
% provided content isn't clipped.
set(hAx, 'Position', [1, 1, axWidth, axHeight]);

% Force MATLAB to render/update the layout based on current properties.
% This is crucial for an accurate TightInset calculation.
drawnow; 

% Get TightInset: [left, bottom, right, top] margins in pixels
% These margins represent the space needed around the plotArea for labels, titles, etc.
ti = get(hAx, 'TightInset'); 

% Calculate figure dimensions needed to encompass the plot area and its decorations
figWidth = axWidth + ti(1) + ti(3);
figHeight = axHeight + ti(2) + ti(4);

% Set figure size. Preserve original screen position (left, bottom) of the figure.
fhUnits = get(hFig, 'Units');
set(hFig, 'Units', 'pixels');
figPos = get(hFig, 'Position'); % To preserve screen x,y for the figure
set(hFig, 'Position', [figPos(1), figPos(2), figWidth, figHeight]);

% Force MATLAB to render/update the figure with its new size.
drawnow; 

% Position the axes within the now-resized figure.
% The axes plot area (Position) should start at ti(1) from the figure's client area left edge,
% and ti(2) from the figure's client area bottom edge.
set(hAx, 'Position', [ti(1), ti(2), axWidth, axHeight]);

% Restore original units for figure and axes
set(hFig, 'Units', fhUnits);
set(hAx, 'Units', axUnits);

end % EOF