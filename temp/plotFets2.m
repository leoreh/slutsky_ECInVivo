function plotFets(fetclu, rpv, grp, iclu, rmflag, rpvRatioNew, rmSpkIdx)

nfet = size(fetclu, 2);
[nsub] = numSubplots(nfet);
nbins = 100;
normmode = 'probability';
figure('units','normalized','outerposition',[0 0 1 1])
for ifet = 1 : nfet
    subplot(nsub(1), nsub(2), ifet)
    histogram(fetclu(:, ifet), nbins,...
        'FaceAlpha', 0.4, 'LineStyle', 'none', 'Normalization', normmode)
    hold on
    histogram(fetclu(rpv, ifet), nbins,...
        'FaceAlpha', 0.4, 'LineStyle', 'none', 'Normalization', normmode)
    yLim = ylim;
    if any(rmflag(ifet, :))
        for iflag = 1 : sum(rmflag(ifet, :) ~= 0)
            plot([rmflag(ifet, iflag) rmflag(ifet, iflag)], yLim, '--k', 'LineWidth', 3)
        end
    end
    title(sprintf('fet #%d', ifet))
    if ifet == 1
        legend({'All spks', 'RPVs'})
    end
end
suptitle(sprintf('T#%d clu#%d; RPV from %.2f to %.2f; removed %d%% of spks', grp, iclu,...
    length(rpv) / length(fetclu) * 100, rpvRatioNew, round(length(rmSpkIdx) / length(fetclu) * 100)))
mkdir(fullfile('graphics', 'clusterClean'))    
figname = fullfile('graphics', 'clusterClean', sprintf('T%d_clu%d', grp, iclu));
export_fig(figname)

end

