function plot_silhouette(psd_cell, clr, snames, axh)

% PCA
psdMat = vertcat(psd_cell{:});
state_idx = []; clr_idx = [];
for istate = 1 : length(psd_cell)
    state_idx = [state_idx; repmat(istate, size(psd_cell{istate}, 1), 1)];
    clr_idx = [clr_idx; repmat(clr{istate}, size(psd_cell{istate}, 1), 1)];
end
[~, pc, ~, ~, ~] = pca(psdMat, 'NumComponents', 3);  % 3 PCs for visualization

% plot
silhouette(pc, state_idx, 'Euclidean');
xlim([-1 1])
yticklabels(snames)
ylabel('')
sh = get(gca, 'Children');
sh.FaceColor = 'flat';
sh.CData(~isnan(sh.YData), :) = clr_idx;
title(axh, 'Separation Quality')
yTickVals = yticks;

% add mean values to plot
silVals = silhouette(pc, state_idx, 'Euclidean');
for istate = 1 : length(psd_cell)
    meanSil = mean(silVals(state_idx == istate));
    if ~isnan(meanSil)
        text(-0.75, yTickVals(istate), num2str(meanSil, '%.3f'));
    end
end

end