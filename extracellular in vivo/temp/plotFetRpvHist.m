% -------------------------------------------------------------------------
% plot histogram of features for spks and rpvs
function plotFetRpvHist(fetclu, cluid, grpid, rpv, rmvSpkIdx, rmvflag, visible)

nfet = size(fetclu, 2);
[nsub] = numSubplots(nfet);
nbins = 100;
normmode = 'probability';
figure('units','normalized','outerposition',[0 0 1 1], 'visible', visible)
for ifet = 1 : nfet
    subplot(nsub(1), nsub(2), ifet)
    histogram(fetclu(:, ifet), nbins,...
        'FaceAlpha', 0.4, 'LineStyle', 'none', 'Normalization', normmode)
    hold on
    histogram(fetclu(rpv, ifet), nbins,...
        'FaceAlpha', 0.4, 'LineStyle', 'none', 'Normalization', normmode)
    yLim = ylim;
    if rmvflag(ifet, 1) ~= 0
        plot([rmvflag(ifet, 1) rmvflag(ifet, 1)], yLim, '--k', 'LineWidth', 3)
    end
    if rmvflag(ifet, 2) ~= 0
        plot([rmvflag(ifet, 2) rmvflag(ifet, 2)], yLim, '--k', 'LineWidth', 3)
    end
    title(sprintf('fet #%d\nremoved %d spikes', ifet, length(rmvSpkIdx)))
    if ifet == 1
        legend({'All spks', 'RPVs'})
    end
end
suptitle(sprintf('T#%d clu#%d; RPV from %.2f to %.2f; removed %d%% of spks', grpid, cluid,...
    length(rpv) / length(fetclu) * 100, rpvRatioNew, round(length(unique(rmvSpkIdx)) / length(fetclu) * 100)))
if strcmp(visible, 'off')
    mkdir(fullfile('graphics', 'clusterClean'))
    figname = fullfile('graphics', 'clusterClean', sprintf('T%d_clu%d', grpid, cluid));
    export_fig(figname)
end
end