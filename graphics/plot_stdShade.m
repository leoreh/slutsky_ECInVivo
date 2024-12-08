function plot_stdShade(varargin)

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

% ensure the axis handle is held
axes(axh);
hold on;

% calculate the mean and standard deviation
mData = mean(dataMat, 2, 'omitnan');
sData = std(dataMat, 0, 2, 'omitnan');

% omit nan values
validIdx = ~isnan(mData) & ~isnan(sData);
mData = mData(validIdx);
sData = sData(validIdx);
xVal = xVal(validIdx)';

% correct special case where std is zero
sData(sData == 0) = eps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot the mean
plot(axh, xVal, mData, 'Color', clr, 'LineWidth', 2);

% plot the shaded area for standard deviation
fillh = fill([xVal; flipud(xVal)], [mData + sData; flipud(mData - sData)], ...
    clr, 'FaceAlpha', alpha, 'EdgeColor', 'none');
set(fillh, 'HandleVisibility', 'off');

end

% EOF