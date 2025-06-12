function hndSct = plot_scatterCorr(x, y, varargin)
% PLOT_SCATTERCORR Creates a scatter plot with correlation statistics and histograms.
%
% SUMMARY:
% This function creates a scatter plot with marginal histograms and correlation
% statistics. It supports outlier removal, data validation, and automatic log
% scaling based on data skewness.
%
% INPUT (Required):
%   x           - Vector of x-axis data
%   y           - Vector of y-axis data
%
% INPUT (Optional Key-Value Pairs):
%   xLbl        - X-axis label {''}
%   yLbl        - Y-axis label {''}
%   flgOtl      - Logical flag to remove outliers {false}
%   plotType    - Type of plot: 'scatter' or 'scatterHist' {'scatter'}
%   hndAx       - Handle to axes for plotting {[]}
%
% OUTPUT:
%   hndSct      - Handle to the plot object
%
% DEPENDENCIES:
%   None
%
% HISTORY:
%   Aug 2024 (AI Assisted) - Initial version

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING & INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define input parameters and their defaults/validation functions
p = inputParser;
addRequired(p, 'x', @isnumeric);
addRequired(p, 'y', @isnumeric);
addParameter(p, 'xLbl', '', @ischar);
addParameter(p, 'yLbl', '', @ischar);
addParameter(p, 'flgOtl', false, @islogical);
addParameter(p, 'plotType', 'scatter', @(x) ismember(x, {'scatter', 'scatterHist'}));
addParameter(p, 'hndAx', [], @(x) isempty(x) || isgraphics(x));

% Parse input arguments
parse(p, x, y, varargin{:});
x = p.Results.x;
y = p.Results.y;
xLbl = p.Results.xLbl;
yLbl = p.Results.yLbl;
flgOtl = p.Results.flgOtl;
plotType = p.Results.plotType;
hndAx = p.Results.hndAx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA VALIDATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get valid indices and log scaling flags
[x, y, flgLog] = validate_data(x, y, flgOtl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMMON CALCULATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate statistics
[r, p] = corr(x, y, 'Type', 'Spearman', 'Rows', 'complete');

% Create statistics string
stats_str = sprintf('n = %d; r = %.2f; p = %.3f', ...
    sum(~isnan(x) & ~isnan(y)), r, p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT CREATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch plotType
    case 'scatter'
        % Select axes if provided
        if isempty(hndAx)
            hndAx = gca;
        end
        
        % Create simple scatter plot
        hndSct = scatter(hndAx, x, y, 20, [0.3 0.3 0.3], '.');
        hold(hndAx, 'on');
        
        % Apply log scaling if needed
        if flgLog(1)
            set(hndAx, 'XScale', 'log');
        end
        if flgLog(2)
            set(hndAx, 'YScale', 'log');
        end
        
        % Add title and labels
        title(hndAx, stats_str);
        xlabel(hndAx, xLbl);
        ylabel(hndAx, yLbl);
        
    case 'scatterHist'
        % Create scatter histogram
        binData = histcounts(x, 'BinMethod', 'fd');
        nBins = length(binData);
        
        sh_args = {'Color', [0.3 0.3 0.3], ...
            'MarkerStyle', '.', ...
            'MarkerSize', 20,...
            'NumBins', nBins};

        if ~isempty(hndAx) && isa(hndAx.Parent, 'matlab.graphics.layout.TiledChartLayout')
            parent = hndAx.Parent;
            delete(hndAx); % Free up the tile for scatterhistogram
            hndSct = scatterhistogram(parent, x, y, sh_args{:});
        else
            hndSct = scatterhistogram(x, y, sh_args{:});
        end
        
        % Add title and labels
        title(hndSct, {stats_str});
        xlabel(hndSct, xLbl);
        ylabel(hndSct, yLbl);
        
        % Adjust appearance
        hndSct.HistogramDisplayStyle = 'bar';
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION: VALIDATE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, y, flgLog] = validate_data(x, y, flgOtl)
% VALIDATE_DATA Validates and optionally removes outliers from input data.
% Also tests skewness to determine if log scaling should be applied.
%
% INPUT:
%   x           - Vector of x-axis data
%   y           - Vector of y-axis data
%   flgOtl      - Logical flag to remove outliers
%
% OUTPUT:
%   x           - Cleaned x data
%   y           - Cleaned y data
%   flgLog      - 2-element logical vector [x_log, y_log] indicating if log scale should be used

% Ensure column vectors
x = x(:);
y = y(:);

% Check for NaN and Inf
validIdx = ~isnan(x) & ~isnan(y) & ~isinf(x) & ~isinf(y);

% Remove outliers if requested
if flgOtl
    [~, xOtl] = rmoutliers(x(validIdx), 'percentiles', [1 99]);
    [~, yOtl] = rmoutliers(y(validIdx), 'percentiles', [1 99]);
    
    % Update valid indices
    validIdx(validIdx) = ~xOtl & ~yOtl;
end

% Apply valid indices
x = x(validIdx);
y = y(validIdx);

% Test skewness and determine log scaling. Apply log scale if skewness > 2
% and values are non-negative. Zeros are handled by replacing them with
% half the minimum non-zero value, if they are not too frequent.
flgLog = [false, false];
thrSkew = 2;
thr0 = 0.1;

% Process x-data
if min(x) >= 0 % Data must be non-negative
    % Handle zeros if they are a small proportion of the data
    if any(x == 0) && (sum(x == 0) < thr0 * length(x))
        min_val = min(x(x > 0));
        if ~isempty(min_val)
            x(x == 0) = min_val / 2;
        end
    end
    
    % Check for log scale if all values are positive
    if all(x > 0)
        skewX = skewness(x);
        flgLog(1) = abs(skewX) > thrSkew;
    end
end

% Process y-data
if min(y) >= 0 % Data must be non-negative
    % Handle zeros if they are a small proportion of the data
    if any(y == 0) && (sum(y == 0) < thr0 * length(y))
        min_val = min(y(y > 0));
        if ~isempty(min_val)
            y(y == 0) = min_val / 2;
        end
    end
    
    % Check for log scale if all values are positive
    if all(y > 0)
        skewY = skewness(y);
        flgLog(2) = abs(skewY) > thrSkew;
    end
end

end 