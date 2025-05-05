function ph = plot_stdShade(varargin)

% plots the mean and standard deviation as a shaded area of the data matrix.
% allows specifying the color and transparency for the plot.
%
% 26 may 24 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'dataMat', [], @isnumeric);
addOptional(p, 'xVal', @isnumeric);
addOptional(p, 'axh', @(x) isa(x, 'matlab.graphics.axis.Axes'));
addOptional(p, 'clr', [0, 0, 1], @(x) isnumeric(x) && length(x)==3);
addOptional(p, 'alpha', 0.5, @isnumeric);

parse(p, varargin{:});
dataMat = p.Results.dataMat;
xVal = p.Results.xVal;
axh = p.Results.axh;
clr = p.Results.clr;
alpha = p.Results.alpha;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check orientation
if size(dataMat, 2) == length(xVal)
    % Preferred: columns = x-axis
    dataMatPlot = dataMat;
elseif size(dataMat, 1) == length(xVal)
    % Transpose: rows = x-axis
    dataMatPlot = dataMat';
    warning('plot_stdShade:Transposed', ...
        'Transposing dataMat so that columns correspond to x-axis values.');
else
    error('plot_stdShade:DimensionMismatch', ...
        'Neither dimension of dataMat matches length(xVal).');
end

% ensure the axis handle is held
axes(axh);
hold on;

% calculate the mean and standard deviation
mData = mean(dataMatPlot, 1, 'omitnan');
n = sum(~isnan(dataMatPlot), 1);  % number of non-NaN points per column
sData = std(dataMatPlot, 0, 1, 'omitnan') ./ sqrt(n);

% omit nan values
validIdx = ~isnan(mData) & ~isnan(sData);
mData = mData(validIdx);
sData = sData(validIdx);
xVal = xVal(validIdx);
xVal = xVal(:);

% correct special case where std is zero
sData(sData == 0) = eps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot the mean
ph = plot(axh, xVal, mData, 'Color', clr, 'LineWidth', 2);

% plot the shaded area for standard deviation
% Check if y-axis is log scale
if strcmp(get(axh, 'YScale'), 'log')
    lowerBound = max(mData - sData, eps);
    upperBound = max(mData + sData, eps);
else
    lowerBound = mData - sData;
    upperBound = mData + sData;
end

% Ensure all are column vectors
xVal = xVal(:);
upperBound = upperBound(:);
lowerBound = lowerBound(:);

fillh = fill([xVal; flipud(xVal)], [upperBound; flipud(lowerBound)], ...
    clr, 'FaceAlpha', alpha, 'EdgeColor', 'none', 'Tag', 'sePatch');

end

% EOF