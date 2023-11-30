function plot_boxMean(varargin)

% plots a box plot (with patch color) of the data mat and the average.
% assumes columns are different groups
%
% 09 jan 22 LH  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'dataMat', [], @isnumeric);
addOptional(p, 'clr', 'k');
addOptional(p, 'alphaIdx', []);
addOptional(p, 'allPnts', false, @islogical);

parse(p, varargin{:})
dataMat = p.Results.dataMat;
clr = p.Results.clr;
alphaIdx = p.Results.alphaIdx;
allPnts = p.Results.allPnts;

ngrp = size(dataMat, 2);
if isempty(clr)
    clr = repmat('k', 1, ngrp);
end
clr = clr(:);
if size(clr, 2) ~= ngrp
    clr = repmat(clr, 1, ngrp);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

boxplot(dataMat, 'PlotStyle', 'traditional', 'Whisker', 2);
bh = findobj(gca, 'Tag', 'Box');
bh = flipud(bh);
if isempty(alphaIdx)
    alphaIdx = linspace(0.1, 0.3, ngrp);
end
for ibox = 1 : length(bh)
    patch(get(bh(ibox), 'XData'), get(bh(ibox), 'YData'),...
        clr(:, ibox)', 'FaceAlpha', alphaIdx(ibox))
end
hold on
if allPnts
    for ibox = 1 : length(bh)
        jitter = (0.2) .* rand(size(dataMat, 1), 1) - 0.1;
        xval = ones(size(dataMat(:, ibox), 1), 1) * ibox + jitter;
        ph = plot(xval, dataMat(:, ibox),...
            '.', 'MarkerSize', 15);
        ph.Color = 'k';
    end
else
    plot(1 : ngrp, mean(dataMat, 1, 'omitnan'), 'kd', 'markerfacecolor', 'k',...
    'MarkerSize', 10)
end

end