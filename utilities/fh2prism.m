function [fhData, varNames] = fh2prism(fh)
% extracts mean and standard error data from a figure handle containing
% either shaded error plots created by plot_stdShade or bar plots with error bars. 
% Returns a structure with data ready for plotting in Prism.
% Note the order of the column pairs does not necassarily fit that of the
% original plotting
%
% 11 jan 24 LH
% modified to handle bar plots and extract variable names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varNames = {};

% First determine if we're dealing with patches or bars
patchObjs = findobj(fh, 'Type', 'patch');
barObjs = findobj(fh, 'Tag', 'barError');

if ~isempty(patchObjs)
    % Handle patch objects (line plots with shaded errors)
    nObjects = length(patchObjs);
    % extract data from first patch to get dimensions
    [x, m, se] = get_patchData(patchObjs(1));
    nPoints = length(x);
    
    % initialize output matrix
    % columns: x, mean1, se1, nan, mean2, se2, nan, ...
    fhData = nan(nPoints, 1 + nObjects * 2);
    fhData(:, 1) = x;
    
    % fill in mean, se for each patch
    for iObj = 1 : nObjects
        [~, m, se] = get_patchData(patchObjs(iObj));
        colIdx = 2 + (iObj - 1) * 2; % starting column for this object
        fhData(:,colIdx:colIdx+1) = [m se];
    end
    
elseif ~isempty(barObjs)
    % Handle bar objects
    nObjects = length(barObjs);
    % extract data from first bar to get dimensions
    [x, m, se] = get_barData(barObjs(1));
    nPoints = length(x);
    
    % initialize output matrix
    % columns: x, mean1, se1, nan, mean2, se2, nan, ...
    fhData = nan(nPoints, 1 + nObjects * 2);
    fhData(:, 1) = x;
    
    % fill in mean, se for each bar group
    for iObj = 1 : nObjects
        [~, m, se] = get_barData(barObjs(iObj));
        colIdx = 2 + (iObj - 1) * 2; % starting column for this object
        fhData(:,colIdx:colIdx+1) = [m se];
    end
end

% varNames = get_varNames(fh);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, meanVal, SE] = get_patchData(patchObj)
% extracts mean and SE from a patch object created by plot_stdShade
vertices = patchObj.Vertices;
nPoints = size(vertices,1)/2;
% extract x values (first half since they go forward then backward)
x = vertices(1:nPoints,1);
% get bounds (upper goes forward, lower goes backward)
upperBound = vertices(1:nPoints,2);
lowerBound = flipud(vertices(nPoints+1:end,2));
% calculate mean and SE
SE = (upperBound - lowerBound)/2;

ax = get(patchObj, 'Parent');
lineObjs = findobj(ax, 'Type', 'Line');
% First, try to match by color (assuming each patch and its mean line share the same color)
patchColor = get(patchObj, 'FaceColor');
meanVal = [];
for k = 1:length(lineObjs)
    if norm(lineObjs(k).Color - patchColor) < 1e-3
        meanVal = lineObjs(k).YData(:);
        break;
    end
end
% Fallback: use axes Children order (if the line was plotted immediately before the patch)
if isempty(meanVal)
    children = get(ax, 'Children');
    idx = find(children == patchObj, 1);
    if idx > 1 && strcmp(get(children(idx-1), 'Type'), 'line')
        meanVal = children(idx-1).YData(:);
    else
        meanVal = (upperBound + lowerBound) / 2; % fallback if no matching line is found
    end
end


end

function [x, meanVal, SE] = get_barData(barObj)
% Extracts data from error bar objects that were tagged with 'barError'

% Get x positions
x = barObj.XData;
if size(x,1) < size(x,2)  % ensure column vector
    x = x';
end

% Get mean values (the actual data points)
meanVal = barObj.YData;
if size(meanVal,1) < size(meanVal,2)
    meanVal = meanVal';
end

% Get standard error values
SE = barObj.YPositiveDelta;  % assuming symmetric error bars
if size(SE,1) < size(SE,2)
    SE = SE';
end

% Handle cases where data might be in a matrix (multiple groups)
if size(meanVal,2) > 1
    meanVal = meanVal(:);
    SE = SE(:);
    x = repmat(x, size(meanVal,2), 1);
end
end


function varNames = get_varNames(fh)
% Returns a cell array of variable names from the figure
varNames = {};

% Get the legend object
legendObj = findobj(fh, 'Type', 'Legend');
if ~isempty(legendObj)
    % Get strings from legend
    varNames = legendObj.String;
end

% get axes labels
axObj = findobj(fh, 'Type', 'Axes');
if ~isempty(axObj)
    xLabel = get_varNames(get_varNames(axObj(1), 'XLabel'), 'String');
    yLabel = get_varNames(get_varNames(axObj(1), 'YLabel'), 'String');
    varNames = [varNames, {xLabel, yLabel}];
end

% Convert to vertical cell array
varNames = varNames(:);

end