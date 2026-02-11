function [hFit, fitStats] = plot_lineReg(x, y, varargin)
% PLOT_LINEREG Adds a linear or orthogonal regression line to the current axis.
%
%   Adapts the regression model based on the axis scale:
%   - Linear-Linear: Fits y = mx + c (Linear).
%   - Log-Log:       Fits log(y) = m*log(x) + c (Power Law: y = 10^c * x^m).
%   - Semilog-Y:     Fits log(y) = mx + c (Exponential).
%   - Semilog-X:     Fits y = m*log(x) + c (Logarithmic).
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

% Prepare Data (Remove NaNs and Infs)
x = x(:);
y = y(:);
mk = ~isnan(x) & ~isnan(y) & ~isinf(x) & ~isinf(y);
x = x(mk);
y = y(mk);

fitStats = struct('slope', NaN, 'intercept', NaN, 'R2', NaN, 'pVal', NaN);
hFit  = gobjects(0); 

if isempty(x) || length(x) < 2
    return;
end

axes(hAx);
hold(hAx, 'on');

% --- Determine Scales & Transform ---
isLogX = strcmp(hAx.XScale, 'log');
isLogY = strcmp(hAx.YScale, 'log');

% Filter for Positive Values if Log Scale
if isLogX, mk = x > 0; x = x(mk); y = y(mk); end
if isLogY, mk = y > 0; x = x(mk); y = y(mk); end

if isempty(x) || length(x) < 2
    return;
end

% Transform Data for Fitting
% We fit the model in the transformed space so the line appears straight.
if isLogX, x = log10(x); end
if isLogY, y = log10(y); end

% --- Determine Range for Plotting ---
% Get current limits (in raw space)
xLims = xlim(hAx);
xEvalRaw = xLims;

% Generate Evaluation Points (Dense for curves)
if isLogX
    xEvalLog = linspace(log10(xEvalRaw(1)), log10(xEvalRaw(2)), 100);
    xEvalFit = xEvalLog; % This is what we pass to polyval (transformed space)
    xPlot    = 10.^xEvalLog; % This is x-coordinate for plot
else
    xEvalFit = linspace(xEvalRaw(1), xEvalRaw(2), 100);
    xPlot    = xEvalFit;
end

%% ========================================================================
%  PLOT & STATS
%  ========================================================================

statsStr = '';

switch lower(regType)
    case 'linear'
        % Linear Regression (Least Squares)
        % Fits y = mx + c in the transformed space
        [pPoly, ~] = polyfit(x, y, 1);
        slope = pPoly(1);
        intercept = pPoly(2);
        
        yEvalFit = polyval(pPoly, xEvalFit);
        
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
        % Slope in transformed space (e.g. exponent in Power Law)
        statsStr = sprintf('Slope = %.2f, R^2 = %.2f', slope, R2);
        
    case 'ortho'
        % Orthogonal Regression (PCA)
        dataXY = [x, y];
        [coeff, ~, ~] = pca(dataXY);
        v1 = coeff(:, 1); % First Principal Component
        mu = mean(dataXY);
        
        slope = v1(2) / v1(1);
        yInt  = mu(2) - slope * mu(1);
        
        yEvalFit = slope .* xEvalFit + yInt;
        
        fitStats.slope = slope;
        fitStats.intercept = yInt;
        fitStats.angle = atan2d(v1(2), v1(1));
        
        % Format String
        statsStr = sprintf('Slope = %.2f (%.1f%c)', ...
            fitStats.slope, fitStats.angle, char(176));
end

% Transform Y back for plotting
if isLogY
    yPlot = 10.^yEvalFit;
else
    yPlot = yEvalFit;
end

% Draw Line
hFit = plot(hAx, xPlot, yPlot, '-', 'Color', clr, 'LineWidth', 2, ...
    'HandleVisibility', 'on', 'DisplayName', 'Fit'); 

% Update DisplayName with fitStats if requested
if flgTxt
    set(hFit, 'DisplayName', statsStr);
else
    hFit.HandleVisibility = 'off';
end

end
