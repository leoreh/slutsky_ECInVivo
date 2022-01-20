isession = 1;
stateIdx = 1 : 5;
stateMfr_tmp = cellfun(@(x) mean(x, 2), fr.states.fr, 'uni', false);
stateMfr = cell2mat(stateMfr_tmp(stateIdx))';

wakeIdx = fr.states.tstamps{1} < 29729;
stateMfr(1, :) = mean(fr.states.fr{1}(:, wakeIdx), 2);
stateMfr(6, :) = mean(fr.states.fr{1}(:, ~wakeIdx), 2);

mfr = mean(fr.strd, 2);

frBoundries = [0.1, Inf; 0.1 Inf];
units = selectUnits('basepath', basepath, 'grp', [],...
    'forceA', true, 'frBoundries', frBoundries,...
    'stableFlag', false);


unitType = {'RS', 'FS'};
cfg = as_loadConfig();
nstates = cfg.nstates;
cfg.colors

fh = figure;
for iunit = 1 : 2
    
    % firing rate median
    medFr = median(mfr(units(iunit, :)));
    yLimit = ceil([0 prctile(stateMfr(:, units(iunit, :)), 98 ,'all')]);
    
    % high firing units
    medUnits = units(iunit, :) & (mfr > medFr)';
    subplot(2, 2, iunit + iunit - 1)
    dataMat = stateMfr(:, medUnits);
    plot_boxMean('dataMat', dataMat', 'clr', 'k')
    bh = findobj(gca, 'Tag', 'Box');
    bh = flipud(bh);
    for ibox = 1 : length(bh) - 1
        patch(get(bh(ibox), 'XData'), get(bh(ibox), 'YData'),...
            cfg.colors{stateIdx(ibox)}, 'FaceAlpha', 0.5)
    end
    patch(get(bh(6), 'XData'), get(bh(6), 'YData'),...
        [0.8 0.8 0.8], 'FaceAlpha', 0.5)
    ylabel('Firing Rate [Hz]')
    xticklabels([cfg.names(stateIdx); {'RUN'}])
    subtitle(sprintf('High MFR %s = %d', unitType{iunit}, sum(medUnits)))
    ylim(yLimit)
    
    % low firing units
    medUnits = units(iunit, :) & (mfr < medFr)';
    subplot(2, 2, iunit + 1 + iunit - 1)
    dataMat = stateMfr(:, medUnits);
    plot_boxMean('dataMat', dataMat', 'clr', 'k')
    bh = findobj(gca, 'Tag', 'Box');
    bh = flipud(bh);
    for ibox = 1 : length(bh) - 1
        patch(get(bh(ibox), 'XData'), get(bh(ibox), 'YData'),...
            cfg.colors{stateIdx(ibox)}, 'FaceAlpha', 0.5)
    end
    patch(get(bh(6), 'XData'), get(bh(6), 'YData'),...
        [0.8 0.8 0.8], 'FaceAlpha', 0.5)
    ylabel('Firing Rate [Hz]')
    xticklabels([cfg.names(stateIdx); {'RUN'}])
    subtitle(sprintf('Low MFR %s = %d', unitType{iunit}, sum(medUnits)))
    ylim(yLimit)
    
end