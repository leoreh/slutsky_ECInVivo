function Stats = plot_linReg(x, y, varargin)
% PLOT_LINREG Adds a linear or orthogonal regression line to the current axis.
%
%   Stats = PLOT_LINREG(x, y) computes and plots a linear regression line
%   for the data vectors x and y.
%
%   Stats = PLOT_LINREG(..., 'Name', Value) specifies additional options:
%       'hAx'    - (handle)  Target axis handle. Default: gca.
%       'clr'    - (char/vec) Color of the line and text. Default: 'k'.
%       'type'   - (char)    Regression type: 'linear' (default) or 'ortho'.
%       'flgTxt' - (logical) Whether to annotate the slope. Default: true.
%
%   INPUTS:
%       x, y    - (numeric) Data vectors.
%
%   OUTPUTS:
%       Stats   - (struct)  Regression statistics (slope, intercept, etc.).
%
%   See also: PCA, POLYFIT, PLOT

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'x', @isnumeric);
addRequired(p, 'y', @isnumeric);
addParameter(p, 'hAx', gca, @isgraphics);
addParameter(p, 'clr', 'k');
addParameter(p, 'type', 'linear', @(x) any(validatestring(x, {'linear', 'ortho'})));
addParameter(p, 'flgTxt', true, @islogical);
parse(p, x, y, varargin{:});

hAx = p.Results.hAx;
clr = p.Results.clr;
regType = p.Results.type;
flgTxt = p.Results.flgTxt;

%% ========================================================================
%  PREPERATIONS
%  ========================================================================

% Prepare Data (Remove NaNs)
x = x(:);
y = y(:);
mk = ~isnan(x) & ~isnan(y);
x = x(mk);
y = y(mk);

Stats = struct('slope', NaN, 'intercept', NaN, 'R2', NaN, 'pVal', NaN);

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
%  PLOT
%  ========================================================================

switch lower(regType)
    case 'linear'
        % Linear Regression (Least Squares)
        % y = mx + c
        [pPoly, S] = polyfit(x, y, 1);
        slope = pPoly(1);
        intercept = pPoly(2);
        
        yEval = polyval(pPoly, xEval);
        
        % Calculate Stats
        yfit = polyval(pPoly, x);
        yresid = y - yfit;
        SSresid = sum(yresid.^2);
        SStotal = (length(y)-1) * var(y);
        R2 = 1 - SSresid/SStotal;

        % Stats Structure Init
        Stats.slope = slope;
        Stats.intercept = intercept;
        Stats.R2 = R2;
        Stats.hLine = [];
        Stats.hTxt  = [];
        
        % Plot
        hLine = plot(hAx, xEval, yEval, '-', 'Color', clr, 'LineWidth', 2, ...
            'HandleVisibility', 'off');
        Stats.hLine = hLine;
        
        % Annotate
        if flgTxt
           
           txtStr = sprintf('Slope: %.2f\nR^2: %.2f', slope, R2);
           
           % Simple placement logic: Top Left relative to limits
           yLims = ylim(hAx);
           xPos = xLims(1) + diff(xLims) * 0.05;
           yPos = yLims(2) - diff(yLims) * 0.05;
           
           hTxt = text(xPos, yPos, txtStr, 'Color', clr, 'FontSize', 9, ...
               'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
               'FontWeight', 'bold');
           Stats.hTxt = hTxt;
        end
        
    case 'ortho'
        % Orthogonal Regression (PCA)
        dataXY = [x, y];
        [coeff, ~, ~] = pca(dataXY);
        v1 = coeff(:, 1); % First Principal Component
        mu = mean(dataXY);
        
        slope = v1(2) / v1(1);
        yInt  = mu(2) - slope * mu(1);
        
        % Calculate Plot Points
        % y = slope * x + yInt
        yEval = slope .* xEval + yInt;
        
        % Stats Structure Init
        Stats.slope = slope;
        Stats.intercept = yInt;
        Stats.angle = atan2d(v1(2), v1(1));
        Stats.hLine = [];
        Stats.hTxt  = [];
        
        % Plot
        hLine = plot(hAx, xEval, yEval, '-', 'Color', clr, 'LineWidth', 2, ...
            'HandleVisibility', 'off');
        Stats.hLine = hLine;
        
        % Annotate
        if flgTxt
            regAngle = atan2d(v1(2), v1(1));
            
            txtStr = sprintf('Slope: %.2f (%.1f\\circ)', slope, regAngle);
            
            % Placement: Top Left generally works, or adapt to data
            yLims = ylim(hAx);
            xPos = xLims(1) + diff(xLims) * 0.05;
            yPos = yLims(2) - diff(yLims) * 0.05;
                       
            hTxt = text(xPos, yPos, txtStr, 'Color', clr, 'FontSize', 9, ...
                'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
                'FontWeight', 'bold');
            Stats.hTxt = hTxt;
        end
end

end
