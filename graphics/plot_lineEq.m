function [hRef, hZero] = plot_lineEq(varargin)
% PLOT_LINEEQ Plots an equality line (y=x) and optional zero lines.
%
%   PLOT_LINEEQ() plots a y=x line on the current axis, adjusting limits
%   to be symmetric and cover all data.
%
%   PLOT_LINEEQ(..., 'Name', Value) specifies additional options:
%       'hAx'       - (handle) Target axis handle. Default: gca.
%       'lims'      - (scalar) Explicit limit value [-lims, lims].
%                     If empty, calculated from data in hAx.
%       'flgZero'   - (logical) Plot x=0 and y=0 lines. Default: determines
%                     from data (true if data spans zero).
%       'flgSqr'    - (logical) Make axis square. Default: false.
%
%   INPUTS:
%       hAx         - Handle to target axes.
%
%   OUTPUTS:
%       hRef        - (handle) Equality line object.
%       hZero       - (handle) Vector of zero line objects [xLine, yLine].
%
%   See also: PLOT_LINEREG, PLOT_SCAT

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addParameter(p, 'hAx', [], @(x) isempty(x) || isgraphics(x));
addParameter(p, 'lims', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'flgZero', [], @(x) isempty(x) || islogical(x));
addParameter(p, 'flgSqr', false, @islogical);
parse(p, varargin{:});

hAx     = p.Results.hAx;
lims    = p.Results.lims;
flgZero = p.Results.flgZero;
flgSqr  = p.Results.flgSqr;

if isempty(hAx), hAx = gca; end

%% ========================================================================
%  DATA EXTRACTION & LIMITS
%  ========================================================================

vals = [];
if isempty(lims) || isempty(flgZero)
    % Extract data from all children to find max range
    chld = get(hAx, 'Children');
    
    for iChld = 1:length(chld)
        if isprop(chld(iChld), 'XData') && ~isempty(chld(iChld).XData)
            vals = [vals; chld(iChld).XData(:)]; %#ok<AGROW>
        end
        if isprop(chld(iChld), 'YData') && ~isempty(chld(iChld).YData)
            vals = [vals; chld(iChld).YData(:)]; %#ok<AGROW>
        end
    end
    
    % If no data, default to current limits or [0 1]
    if isempty(vals)
        currX = xlim(hAx);
        currY = ylim(hAx);
        vals = [currX(:); currY(:)];
    end
end

% Default flgZero: True if data spans zero (has negative and positive values)
if isempty(flgZero)
    flgZero = min(vals) < 0 && max(vals) > 0;
end

% Calculate symmetric limits if not provided
if isempty(lims)
    maxVal = max(abs(vals), [], 'omitnan');
    if maxVal == 0, maxVal = 1; end
    lims = maxVal * 1.1; 
end

%% ========================================================================
%  PLOTTING
%  ========================================================================

hold(hAx, 'on');

% Set Limits
xlim(hAx, [-lims, lims]);
ylim(hAx, [-lims, lims]);

% Plot Equality Line (y=x)
hRef = plot(hAx, [-lims, lims], [-lims, lims], '--k', ...
    'LineWidth', 1.5, ... 
    'HandleVisibility', 'off', ...
    'Tag', 'EqualityLine');

% Plot Zero Lines
hZero = gobjects(0);
if flgZero
    hZero(1) = xline(hAx, 0, '-.', 'Color', repmat(0.5, 1, 3), ...
        'HandleVisibility', 'off', 'Tag', 'ZeroLineX');
    hZero(2) = yline(hAx, 0, '-.', 'Color', repmat(0.5, 1, 3), ...
        'HandleVisibility', 'off', 'Tag', 'ZeroLineY');
end

% Set Square Axis
if flgSqr
    axis(hAx, 'square');
end

% Ensure consistent style
grid(hAx, 'off');
box(hAx, 'on');

end
