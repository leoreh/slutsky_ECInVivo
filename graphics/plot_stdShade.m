function [hLine, hFill] = plot_stdShade(varargin)

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
addOptional(p, 'hAx', @(x) isa(x, 'matlab.graphics.axis.Axes'));
addOptional(p, 'clr', [0, 0, 1], @(x) isnumeric(x) && length(x)==3);
addOptional(p, 'alpha', 0.5, @isnumeric);
addOptional(p, 'varbose', false, @islogical);

parse(p, varargin{:});
dataMat = p.Results.dataMat;
xVal = p.Results.xVal;
hAx = p.Results.hAx;
clr = p.Results.clr;
alpha = p.Results.alpha;
varbose = p.Results.varbose;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assure columns = x-axis
if size(dataMat, 1) == length(xVal)
    dataMat = dataMat';
    if varbose; warning('[STDSHADE]: Transposing dataMat'); end

elseif size(dataMat, 2) ~= length(xVal)
    error('plot_stdShade:DimensionMismatch', ...
        'Neither dimension of dataMat matches length(xVal).');
end

% calculate the mean and standard eror of mean
mData = mean(dataMat, 1, 'omitnan');
n = sum(~isnan(dataMat), 1);  % number of non-NaN points per column
sData = std(dataMat, 0, 1, 'omitnan') ./ sqrt(n);

% correct special case where std is zero
sData(sData == 0) = eps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axes(hAx);
hold on;

% plot the mean - plot function handles NaNs correctly by default
hLine = plot(hAx, xVal, mData, 'Color', clr, 'LineWidth', 2);

% plot the shaded area for standard deviation
% Identify contiguous non-NaN blocks to avoid filling across NaN gaps
idxValid = ~isnan(mData);
blocks = diff([0, idxValid, 0]);
idxStart = find(blocks == 1);
idxEnd = find(blocks == -1) - 1;

hFill = gobjects(length(idxStart), 1); % pre-allocate graphics array

for iBlock = 1 : length(idxStart)
    
    % Extract data for the current block
    block_xVal = xVal(idxStart(iBlock) : idxEnd(iBlock));
    block_mData = mData(idxStart(iBlock) : idxEnd(iBlock));
    block_sData = sData(idxStart(iBlock) : idxEnd(iBlock));

    % Check if y-axis is log scale
    if strcmp(get(hAx, 'YScale'), 'log')
        lowerBound = max(block_mData - block_sData, eps);
        upperBound = max(block_mData + block_sData, eps);
    else
        lowerBound = block_mData - block_sData;
        upperBound = block_mData + block_sData;
    end
    
    % Ensure all are column vectors
    block_xVal = block_xVal(:);
    upperBound = upperBound(:);
    lowerBound = lowerBound(:);
    
    % Plot the fill for this block
    hFill(iBlock) = fill([block_xVal; flipud(block_xVal)], [upperBound; flipud(lowerBound)], ...
        clr, 'FaceAlpha', alpha, 'EdgeColor', 'none', 'Tag', 'sePatch',...
        'HandleVisibility', 'off');
end

end

% EOF




% % plot the shaded area for standard deviation
% % Check if y-axis is log scale
% if strcmp(get(hAx, 'YScale'), 'log')
%     lowerBound = max(mData - sData, eps);
%     upperBound = max(mData + sData, eps);
% else
%     lowerBound = mData - sData;
%     upperBound = mData + sData;
% end
% 
% % Ensure all are column vectors
% xVal = xVal(:);
% upperBound = upperBound(:);
% lowerBound = lowerBound(:);
% 
% hFill = fill([xVal; flipud(xVal)], [upperBound; flipud(lowerBound)], ...
%     clr, 'FaceAlpha', alpha, 'EdgeColor', 'none', 'Tag', 'sePatch',...
%     'HandleVisibility', 'off');
% 
% end
% 
% % EOF