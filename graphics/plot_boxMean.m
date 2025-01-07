function plot_boxMean(varargin)

% plots a box plot or bar (mean + se) of the data mat and the average.
% assumes columns are different groups. can color code the boxes / bars.
%
% 09 jan 22 LH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'dataMat', [], @isnumeric);
addOptional(p, 'clr', [0, 0, 0]);
addOptional(p, 'alphaIdx', []);
addOptional(p, 'allPnts', false, @islogical);
addOptional(p, 'plotType', 'box', @ischar);
addOptional(p, 'axh', []);

parse(p, varargin{:})
dataMat = p.Results.dataMat;
clr = p.Results.clr;
alphaIdx = p.Results.alphaIdx;
allPnts = p.Results.allPnts;
plotType = p.Results.plotType;      % can be 'box' or 'bar'
axh = p.Results.axh;      

% prepare colors
ngrp = size(dataMat, 2);
if isempty(clr)  
    clr = repmat([0, 0, 0], ngrp, 1);
end
if size(clr, 1) == 1
    clr = repmat(clr, ngrp, 1);
end

% prepare transparency
if isempty(alphaIdx)
    alphaIdx = linspace(0.5, 1, ngrp);
end

% prepare figure axis
if isempty(axh)
    fh = figure;
    axh = gca;
end
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(plotType, 'box')

    boxplot(axh, dataMat, 'PlotStyle', 'traditional', 'Whisker', 2);
    bh = findobj(gca, 'Tag', 'Box');
    bh = flipud(bh);
    for ibox = 1 : length(bh)
        patch(get(bh(ibox), 'XData'), get(bh(ibox), 'YData'),...
            clr(ibox, :), 'FaceAlpha', alphaIdx(ibox))
    end

elseif strcmp(plotType, 'bar')
    means = mean(dataMat, 1, 'omitnan');
    sem = std(dataMat, [], 1, 'omitnan') / sqrt(size(dataMat, 1)); 
    bh = bar(axh, 1 : ngrp, means);
    bh.CData = clr;
    bh.FaceColor = 'flat';
    errorbar(1 : ngrp, means, sem, 'k', 'linestyle', 'none'); 
    bh.FaceColorMode = 'manual'

end

if allPnts
    for igrp = 1 : ngrp
        jitter = (0.2) .* rand(size(dataMat, 1), 1) - 0.1;
        xval = ones(size(dataMat(:, igrp), 1), 1) * igrp + jitter;
        ph = plot(axh, xval, dataMat(:, igrp),...
            '.', 'MarkerSize', 5);
        ph.Color = clr(igrp, :);
        ph.MarkerEdgeColor = 'k';
    end
elseif strcmp(plotType, 'box')
    plot(axh, 1 : ngrp, mean(dataMat, 1, 'omitnan'), 'kd',...
        'markerfacecolor', 'k', 'MarkerSize', 10)
end

% adjust ylimit
% dataRange = [min(dataMat(:)), max(dataMat(:))];
% pad = diff(dataRange) * 0.5;
% ylim([dataRange(1) - pad, dataRange(2) + pad]);

% adjust xlimit
xlim([1 - 0.5, ngrp + 0.5])
xticks([1 : ngrp]);

% add statistical comparison
pVal = stat_compare1D(dataMat, 'axh', axh);

end