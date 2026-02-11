function [hFit, fitStats] = plot_lineReg(x, y, varargin)
% PLOT_LINEREG Adds a linear or orthogonal regression line to the current axis.
%
%   [hFit, fitStats] = PLOT_LINEREG(x, y) computes and plots a regression line.
%
%   fitStats = PLOT_LINEREG(..., 'Name', Value) specifies additional options:
%       'hAx'    - (handle)  Target axis handle. Default: gca.
%       'clr'    - (char/vec) Color of the line. Default: 'k'.
%       'type'   - (char)    Regression type: 'linear' (default) or 'ortho'.
%       'flgTxt' - (logical) Appends stats to DisplayName (for legend). Default: true.
%
%   INPUTS:
%       x, y    - (numeric) Data vectors.
%
%   OUTPUTS:
%       hFit    - (handle)  Line object.
%       fitStats   - (struct)  Regression statistics (slope, intercept, R2/angle).
%
%   See also: PCA, POLYFIT, PLOT

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'x', @isnumeric);
addRequired(p, 'y', @isnumeric);
addParameter(p, 'hAx', [], @(x) isempty(x) || isgraphics(x));
addParameter(p, 'clr', 'k');
addParameter(p, 'type', 'linear', @(x) any(validatestring(x, {'linear', 'ortho'})));
addParameter(p, 'flgTxt', true, @islogical);
parse(p, x, y, varargin{:});

hAx = p.Results.hAx;
clr = p.Results.clr;
regType = p.Results.type;
flgTxt = p.Results.flgTxt;

if isempty(hAx), hAx = gca; end

%% ========================================================================
%  PREPARATION
%  ========================================================================

% Prepare Data (Remove NaNs)
x = x(:);
y = y(:);
mk = ~isnan(x) & ~isnan(y);
x = x(mk);
y = y(mk);

fitStats = struct('slope', NaN, 'intercept', NaN, 'R2', NaN, 'pVal', NaN);
hFit  = gobjects(0); 

if isempty(x) || length(x) < 2
    return;
end

axes(hAx);
hold(hAx, 'on');

% Determine range for plotting line
xLims = xlim(hAx);
if all(xLims == [0 1]) && all(ylim(hAx) == [0 1]) && ~ishold(hAx)
    % Likely empty axis, use data range
    xRange = [min(x), max(x)];
    margin = diff(xRange) * 0.1;
    xEval = [min(x)-margin, max(x)+margin];
else
    % Use current axis limits
    xEval = xLims;
end

%% ========================================================================
%  PLOT & STATS
%  ========================================================================

statsStr = '';

switch lower(regType)
    case 'linear'
        % Linear Regression (Least Squares)
        [pPoly, ~] = polyfit(x, y, 1);
        slope = pPoly(1);
        intercept = pPoly(2);
        
        yEval = polyval(pPoly, xEval);
        
        % Calculate fitStats
        yfit = polyval(pPoly, x);
        yresid = y - yfit;
        SSresid = sum(yresid.^2);
        SStotal = (length(y)-1) * var(y);
        R2 = 1 - SSresid/SStotal;

        fitStats.slope = slope;
        fitStats.intercept = intercept;
        fitStats.R2 = R2;
        
        % Format String
        % Slope: .2f, R2: .2f
        statsStr = sprintf('Slope = %.2f, R^2 = %.2f', slope, R2);
        
    case 'ortho'
        % Orthogonal Regression (PCA)
        dataXY = [x, y];
        [coeff, ~, ~] = pca(dataXY);
        v1 = coeff(:, 1); % First Principal Component
        mu = mean(dataXY);
        
        slope = v1(2) / v1(1);
        yInt  = mu(2) - slope * mu(1);
        
        yEval = slope .* xEval + yInt;
        
        fitStats.slope = slope;
        fitStats.intercept = yInt;
        fitStats.angle = atan2d(v1(2), v1(1));
        
        % Format String
        statsStr = sprintf('Slope = %.2f (%.1f%c)', ...
            fitStats.slope, fitStats.angle, char(176));
end

% Draw Line
hFit = plot(hAx, xEval, yEval, '-', 'Color', clr, 'LineWidth', 2, ...
    'HandleVisibility', 'on', 'DisplayName', 'Fit'); 

% Update DisplayName with fitStats if requested
if flgTxt
    set(hFit, 'DisplayName', statsStr);
else
    hFit.HandleVisibility = 'off';
end

end
