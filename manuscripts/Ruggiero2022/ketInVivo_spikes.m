
[basepaths, v, nsessions] = ketInVivo_sessions('ket');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reanalyze something
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cell_metrics = CellExplorer('basepath', basepath);

% itterate
for isession = 1 : nsessions
    basepath = basepaths{isession};
    cd(basepath)
    [mousename, basename] = fileparts(basepath);
    [~, mousename] = fileparts(mousename);

    % number of units per spike group
%     plot_nunits_session('basepath', basepath, 'frBoundries', [])

    % select specific units
    grp = v(isession).units.info.grp;
    units = selectUnits('basepath', basepath, 'grp', grp, 'saveVar', true,...
        'forceA', true, 'frBoundries', [0.05 Inf; 0.05 Inf],...
        'spikes', v(isession).spikes);
    units = v(isession).units;
%     units.clean = units.clean(1, :) | units.clean(2, :);
    
%     % plot fr vs. time
%     plot_FRtime_session('basepath', basepath,...
%         'muFlag', false, 'saveFig', false,...
%         'dataType', 'strd', 'units', units.clean)
%  
%     % state ratio according to mfr percetiles
%     frBins = catfields(v(isession).frBins, 'catdef', 'addim', 'force', false);
%     stateMfr = frBins.states.mfr(:, [1, 4], 1);
%     plot_FRstates_sextiles('stateMfr', stateMfr', 'units', units.clean,...
%         'ntiles', 2, 'saveFig', false)
% 
%     timebins = v(isession).session.general.timebins;
%     timepnt = v(isession).session.general.timepnt;
%     timebins = [15 * 60, timepnt; timepnt, timepnt + 3 * 60 * 60;...
%         timepnt + 3 * 60 * 60, timepnt + 6 * 60 * 60;
%         timepnt + 6 * 60 * 60, timepnt + 9 * 60 * 60;
%         timepnt + 9 * 60 * 60, timepnt + 12 * 60 * 60;
%         timepnt + 12 * 60 * 60, Inf];
%     winBL = [0 timepnt];
%     
%     % update session file
%     sessionfile = fullfile(basepath, [basename, '.session.mat']);
%     session = v(isession).session;
%     session.general.timebins = timebins;
%     session.general.timepnt = timepnt;
%     save(sessionfile, 'session')
% 
%     fr = firingRate(v(isession).spikes.times, 'basepath', basepath,...
%         'graphics', true, 'binsize', 60, 'saveVar', true,...
%         'smet', 'none', 'winBL', winBL, 'winCalc', [0, Inf]);
% 
%     fr_timebins('basepath', pwd,...
%         'forceA', true, 'graphics', true,...
%         'timebins', timebins, 'saveVar', true, 'sstates', [1, 4]);
end

% params
cfg = as_loadConfig();
nstates = cfg.nstates;
unitChar = {'RS', 'FS'};
unitClr = {'b', 'r'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% concate sessions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select states
sstates = [1, 4];

% cat and select units
units = catfields([v.units], 'catdef', 'long', 'force', false);
% units.clean = units.clean(1, :) | units.clean(2, :);

% concate mfr in states
stateMfr = [];
for isession = 1 : nsessions
    frBins = catfields(v(isession).frBins, 'catdef', 'addim', 'force', false);
    stateMfr = [stateMfr; frBins.states.mfr(:, sstates, :)];
end

% mfr
fr = catfields([v.fr], 'catdef', 'long', 'force', false);
mfr_bl = fr.mfr;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mfr states - box plot of gain ratio divided to sextiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_FRstates_sextiles('stateMfr', squeeze(stateMfr(:, :, 2))', 'units', units.clean,...
    'ntiles', 2, 'saveFig', false)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% percent change between bins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% params
istate = 2;
iunit = 2;
ntiles = 2;
itile = 1;

% data
mfr_bl = squeeze(stateMfr(units.clean(iunit, :), istate, 1));
mfr_chg = squeeze(stateMfr(units.clean(iunit, :), istate, 2)) ./ mfr_bl * 100;

% divide units by prctiles
tiles = prctile(mfr_bl, [1 / ntiles : 1 / ntiles : 1 - 1 / ntiles] * 100);
tiles = [0, tiles, Inf];

clear mfrMat chgMat unitsTile
for itile = 1 : ntiles
    unitsTile(:, itile) = mfr_bl >= tiles(itile) & mfr_bl < tiles(itile + 1);
    mfrMat{itile} = mfr_bl(unitsTile(:, itile));
    chgMat{itile} = mfr_chg(unitsTile(:, itile));
end
mfrMat = cell2nanmat(mfrMat, 2);
chgMat = cell2nanmat(chgMat, 2);

fh = figure;
th = tiledlayout(2, length(sstates), 'TileSpacing', 'Compact');
axh = nexttile;
plot_boxMean('dataMat', chgMat, 'clr', 'k', 'allPnts', true)





infidx = isinf(log10(mfr_bl));
mdl = fitlm(log10(mfr_bl(~infidx)), mfr_chg(~infidx),...
    'Exclude', mfr_bl(~infidx) < 0);

fh = figure;
plot(mdl)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% states mfr histogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% params
istate = 1;
ibin = 2;
iunit = 1;

fh = figure;
th = tiledlayout(2, length(sstates), 'TileSpacing', 'Compact');

for istate = 1 : 2
    data = log10(stateMfr(units.clean(iunit, :), istate, ibin));
    hh = histogram(data, 15,...
        'Normalization', 'probability');
    hh.EdgeColor = 'none';
    hh.FaceColor = cfg.colors{sstates(istate)};
    hh.FaceAlpha = 0.3;
    hold on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% box plot of mfr across states divided by median
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fh = figure;
for iunit = 1 : 2
    
    % firing rate median
    medFr = median(mfr_bl(units.clean(iunit, :)));
    yLimit = ceil([0 prctile(stateMfr(:, units.clean(iunit, :)), 98 ,'all')]);
    
    % high firing units
    subplot(2, 2, iunit + iunit - 1)
    medUnits = units.clean(iunit, :) & (mfr_bl > medFr)';
    dataMat = stateMfr(:, medUnits);
    plot_boxMean('dataMat', dataMat', 'clr', 'k')
    bh = findobj(gca, 'Tag', 'Box');
    bh = flipud(bh);
    for ibox = 1 : length(bh)
        patch(get(bh(ibox), 'XData'), get(bh(ibox), 'YData'),...
            cfg.colors{sstates(ibox)}, 'FaceAlpha', 0.5)
    end
    ylabel('Firing Rate [Hz]')
    xticklabels(cfg.names)
    title(sprintf('High MFR %s = %d', unitChar{iunit}, sum(medUnits)))
    ylim(yLimit)
    
    % low firing units
    subplot(2, 2, iunit + 1 + iunit - 1)
    medUnits = units.clean(iunit, :) & (mfr_bl < medFr)';
    dataMat = stateMfr(:, medUnits);
    plot_boxMean('dataMat', dataMat', 'clr', 'k')
    bh = findobj(gca, 'Tag', 'Box');
    bh = flipud(bh);
    for ibox = 1 : length(bh)
        patch(get(bh(ibox), 'XData'), get(bh(ibox), 'YData'),...
            cfg.colors{sstates(ibox)}, 'FaceAlpha', 0.5)
    end
    ylabel('Firing Rate [Hz]')
    xticklabels(cfg.names)
    title(sprintf('Low MFR %s = %d', unitChar{iunit}, sum(medUnits)))
    ylim(yLimit)
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mfr per session at different time points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% params
tbins = [2, 4];
iunit = 1;

chg_raw = nan(length(tbins) + 1, nsessions);
chg_norm = nan(length(tbins), nsessions);
for isession = 1 : nsessions
    unitsClean = v(isession).units.clean(iunit, :);

    chg_raw(1, isession) = mean(v(isession).frBins(1).mfr(unitsClean));
    cnt = 1;
    for ibin = 1 : length(tbins)
        chg_raw(cnt + 1, isession) = mean(v(isession).frBins(tbins(ibin)).mfr(unitsClean), 'omitnan');
        chg_norm(cnt, isession) = mean(v(isession).frBins(tbins(ibin)).mfr(unitsClean) ./...
            v(isession).frBins(1).mfr(unitsClean) * 100, 'omitnan');
        cnt = cnt + 1;
    end
end

fh = figure;
th = tiledlayout(1, 2, 'TileSpacing', 'Compact');
axh = nexttile;
plot_boxMean('dataMat', chg_raw', 'clr', 'k', 'allPnts', true)
ylabel('Mean Firing Rate [Hz]')
axh = nexttile;
plot_boxMean('dataMat', chg_norm', 'clr', 'k', 'allPnts', true)
ylabel('Mean Firing Rate [% Baseline]')

% inspect timebins
v(1).session.general.timebins / 60 / 60


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to prism
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% firing rate vs. time across all units, aligned to point of injection
units = catfields([v.units], 'catdef', 'long', 'force', false);

% firing rate vs. time across all units, aligned to point of injection. can
% extract strd or normalized data
[frMat, timeIdx] = alignFR2pnt('basepaths', basepaths, 'dataType', 'strd');

% grab fr
unitNo = 2;
prism_data = frMat(units.clean(unitNo, :), :)';
prism_tstamps = [1 : length(frMat)]' / 60;
prism_injTime = max(timeIdx) / 60;

% list of units per session
prism_rs = []; prism_fs = [];
for isession = 1 : nsessions
    units = v(isession).units;
    prism_rs = [prism_rs; ones(sum(units.clean(1, :)), 1) * isession];
    prism_fs = [prism_fs; ones(sum(units.clean(2, :)), 1) * isession];
end

% select units based on mfr
iunit = 2;
prism_data = stateMfr(units.clean(iunit, :), :, 1);

