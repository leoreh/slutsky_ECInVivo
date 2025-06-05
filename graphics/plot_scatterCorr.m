function hAx = plot_scatterCorr(x, y, varargin)
% PLOT_SCATTERCORR Creates a scatter plot with correlation statistics and histograms.
%
% SUMMARY:
% This function creates a scatter plot with marginal histograms and correlation
% statistics. It supports outlier removal and data validation.
%
% INPUT (Required):
%   x           - Vector of x-axis data
%   y           - Vector of y-axis data
%
% INPUT (Optional Key-Value Pairs):
%   xLbl        - X-axis label {''}
%   yLbl        - Y-axis label {''}
%   flgOtl       - Logical flag to remove outliers {false}
%   hAx         - Axes handle to plot into {[]}
%
% OUTPUT:
%   hAx         - Handle to the axes containing the plot
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
addParameter(p, 'hAx', [], @(x) isempty(x) || ishandle(x) && strcmp(get(x, 'Type'), 'axes'));

% Parse input arguments
parse(p, x, y, varargin{:});
x = p.Results.x;
y = p.Results.y;
xLbl = p.Results.xLbl;
yLbl = p.Results.yLbl;
flgOtl = p.Results.flgOtl;
hAx = p.Results.hAx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA VALIDATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get valid indices
[x, y] = validate_data(x, y, flgOtl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT CREATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

binData = histcounts(x, 'BinMethod', 'fd');
nBins = length(binData);

% Create scatter histogram
hndSct = scatterhistogram(x, y, ...
    'Color', [0.3 0.3 0.3], ...
    'MarkerStyle', '.', ...
    'MarkerSize', 20,...
    'NumBins', nBins);

% Calculate statistics
[r, p] = corr(x, y, 'Type', 'Spearman', 'Rows', 'complete');

% Add statistics to title
stats_str = sprintf('n = %d; r = %.2f; p = %.3f', ...
    sum(~isnan(x) & ~isnan(y)), r, p);
title({stats_str});

% Add labels
xlabel(xLbl);
ylabel(yLbl);

% Adjust appearance
hndSct.HistogramDisplayStyle = 'bar';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION: VALIDATE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, y] = validate_data(x, y, flgOtl)
% VALIDATE_DATA Validates and optionally removes outliers from input data.
%
% INPUT:
%   x           - Vector of x-axis data
%   y           - Vector of y-axis data
%   flgOtl      - Logical flag to remove outliers
%
% OUTPUT:
%   validIdx    - Logical indices of valid data points
%   x           - Cleaned x data
%   y           - Cleaned y data

% Ensure column vectors
x = x(:);
y = y(:);

% Check for NaN and Inf
validIdx = ~isnan(x) & ~isnan(y) & ~isinf(x) & ~isinf(y);

% Remove outliers if requested
if flgOtl
    % Use MATLAB's built-in outlier detection
    [~, xOtl] = rmoutliers(x(validIdx), 'percentiles', [1 99]);
    [~, yOtl] = rmoutliers(y(validIdx), 'percentiles', [1 99]);
    
    % Update valid indices
    validIdx(validIdx) = ~xOtl & ~yOtl;
end

% Apply valid indices
x = x(validIdx);
y = y(validIdx);

end 