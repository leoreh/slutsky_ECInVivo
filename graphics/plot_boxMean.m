function bh = plot_boxMean(varargin)
% plots a box plot or bar (mean + se) of the data mat and the average.
% handles both single group and multiple group data
%
% INPUT:
%   dataMat     matrix or cell array of matrices. For grouped data, each cell
%               contains a matrix where columns are subgroups
%   xVal        x-axis values (optional, default: 1:n)
%   clr         color matrix, one row per group
%   alphaIdx    transparency values
%   allPnts     logical, whether to plot individual points
%   plotType    'allPnts', 'box' or 'bar'
%   hAx         axis handle
%   grpNames    cell array of group names (optional)
%
% 10 jan 24 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'dataMat', [], @(x) isnumeric(x) || iscell(x));
addOptional(p, 'xVal', [], @isnumeric);
addOptional(p, 'clr', [0, 0, 0]);
addOptional(p, 'alphaIdx', []);
addOptional(p, 'plotType', 'box', @ischar);
addOptional(p, 'hAx', []);
addOptional(p, 'grpNames', []);

parse(p, varargin{:})
dataMat = p.Results.dataMat;
xVal = p.Results.xVal;
clr = p.Results.clr;
alphaIdx = p.Results.alphaIdx;
plotType = p.Results.plotType;
hAx = p.Results.hAx;
grpNames = p.Results.grpNames;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert to cell if matrix input
if isnumeric(dataMat)
    dataMat = {dataMat};
end

% get dimensions
n_grps = length(dataMat);
n_subgrps = cellfun(@(x) size(x,2), dataMat);

% prepare x values
if isempty(xVal)
    xVal = 1 : max(n_subgrps);
end

% calculate offsets for grouped data
if n_grps > 1
    grp_width = 0.8 / n_grps;  % total width of 0.8 for all groups
    offsets = linspace(-0.4 + grp_width/2, 0.4 - grp_width/2, n_grps);
else
    grp_width = 0.8;
    offsets = 0;
end

% prepare colors
if isempty(clr)
    clr = repmat([0, 0, 0], n_grps, 1);
elseif size(clr, 1) == 1
    clr = repmat(clr, n_grps, 1);
end

% prepare transparency
if isempty(alphaIdx)
    alphaIdx = linspace(0.5, 1, n_grps);
end
if isscalar(alphaIdx)
    alphaIdx = repmat(alphaIdx, 1, n_grps);
end

% prepare figure axis
if isempty(hAx)
    fh = figure;
    hAx = gca;
end
hold(hAx, 'on')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot each group
for igrp = 1 : n_grps
    curr_data = dataMat{igrp};
    curr_x = xVal + offsets(igrp);

    % add individual points
    if strcmp(plotType, 'allPnts')
        for isub = 1 : size(curr_data, 1)
            valid_data = curr_data(isub, :);
            valid_data = valid_data(~isnan(valid_data));

            jitter = (grp_width/2) .* rand(length(valid_data), 1) - grp_width/4;
            xval = ones(size(valid_data)) * curr_x(isub) + jitter;

            plot(hAx, xval, valid_data, '.',...
                'MarkerSize', 5, 'Color', clr(igrp, :),...
                'MarkerEdgeColor', clr(igrp, :), 'HandleVisibility', 'off');
        end

        % add mean markers
        bh(igrp) = plot(hAx, curr_x, mean(curr_data, 2, 'omitnan'), 'kd',...
            'markerfacecolor', clr(igrp, :), 'MarkerSize', 10,...
            'HandleVisibility', 'off');

    elseif strcmp(plotType, 'box')
        % box plot
        bh = boxplot(hAx, curr_data', 'Positions', curr_x,...
            'PlotStyle', 'traditional', 'Whisker', 2,...
            'Labels', {}, 'Color', clr(igrp,:));

        % color boxes
        bh = findobj(hAx, 'Tag', 'Box');
        bh = flipud(bh((end - size(curr_data, 1) + 1) : end));

        for ibox = 1:length(bh)
            patch(get(bh(ibox), 'XData'), get(bh(ibox), 'YData'),...
                clr(igrp,:), 'FaceAlpha', alphaIdx(igrp))
        end

    elseif strcmp(plotType, 'bar')
        % calculate stats
        means = mean(curr_data, 2, 'omitnan');
        n = sum(~isnan(curr_data), 2);
        sem = std(curr_data, [], 2, 'omitnan') ./ sqrt(n);

        % plot bars
        bh(igrp) = bar(hAx, curr_x, means, grp_width);
        bh(igrp).CData = clr(igrp,:);
        bh(igrp).FaceColor = 'flat';
        bh(igrp).FaceAlpha = alphaIdx(igrp);

        % add error bars
        errorbar(curr_x, means, [], sem, 'k', 'linestyle', 'none',...
            'tag', 'barError');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finalize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set axis properties
if ~isempty(grpNames)
    legend(grpNames, 'Location', 'best')
end

xlim(hAx, [min(xVal)-0.5, max(xVal)+0.5])
xticks(hAx, xVal)

% remove temporary boxplot legends
hLgd = findobj(hAx, 'Tag', 'Box');
for iLgd = 1 : length(hLgd)
    set(hLgd(iLgd), 'HandleVisibility', 'off');
end

end