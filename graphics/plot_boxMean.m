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
addOptional(p, 'clr', []);
addOptional(p, 'alphaIdx', []);

parse(p, varargin{:})
dataMat = p.Results.dataMat;
clr = p.Results.clr;
alphaIdx = p.Results.alphaIdx;

ngrp = size(dataMat, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(1 : ngrp, mean(dataMat, 1, 'omitnan'), 'kd', 'markerfacecolor', 'k')
hold on
boxplot(dataMat, 'PlotStyle', 'traditional', 'Whisker', 6);
if ~isempty(clr)
    bh = findobj(gca, 'Tag', 'Box');
    bh = flipud(bh);
    if isempty(alphaIdx)
        alphaIdx = linspace(0.5, 0.9, ngrp);
    end
    for ibox = 1 : length(bh)
        patch(get(bh(ibox), 'XData'), get(bh(ibox), 'YData'),...
            clr, 'FaceAlpha', alphaIdx(ibox))
    end
end

end