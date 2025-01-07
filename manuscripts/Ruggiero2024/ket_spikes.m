
[basepaths, v, nfiles] = ketInVivo_sessions('ket');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reanalyze something
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cell_metrics = CellExplorer('basepaths', basepaths);

% itterate
for ifile = 1 : nfiles
    basepath = basepaths{ifile};
    cd(basepath)
    [mousename, basename] = fileparts(basepath);
    [~, mousename] = fileparts(mousename);

    %     % number of units per spike group
    %     plot_nunits_session('basepath', basepath, 'frBoundries', [])
    %
    %     % select specific units
    %     grp = v(ifile).units.info.grp;
    %     grp = [1 : 4];
    %     units = selectUnits('basepath', basepath, 'grp', grp, 'saveVar', true,...
    %         'forceA', true, 'frBoundries', [0.05 Inf; 0.05 Inf],...
    %         'spikes', v(ifile).spikes);
    %     units = v(ifile).units;
    % %     units.clean = units.clean(1, :) | units.clean(2, :);
    % %
    %     % plot fr vs. time
    %     plot_FRtime_session('basepath', basepath,...
    %         'muFlag', false, 'saveFig', false,...
    %         'dataType', 'strd', 'units', units.clean)

    % state ratio according to mfr percetiles
    %     frBins = catfields(v(ifile).frBins, 'catdef', 'addim', 'force', false);
    %     stateMfr = frBins.states.mfr(:, [1, 4], 1);
    %     plot_FRstates_sextiles('stateMfr', stateMfr', 'units', units.clean,...
    %         'ntiles', 2, 'saveFig', false)

    timebins = v(ifile).session.general.timebins;
    timepnt = v(ifile).session.general.timepnt;
    %     timebins = [10 * 60, timepnt; timepnt, timepnt + 3 * 60 * 60;...
    %         timepnt + 3 * 60 * 60, timepnt + 6 * 60 * 60;
    %         timepnt + 6 * 60 * 60, timepnt + 9 * 60 * 60;
    %         timepnt + 9 * 60 * 60, Inf];
    %     winBL = [0 timepnt];
    %
    %     % update session file
    %     sessionfile = fullfile(basepath, [basename, '.session.mat']);
    %     session = v(ifile).session;
    %     session.general.timebins = timebins;
    %     session.general.timepnt = timepnt;
    %     save(sessionfile, 'session')
    %
    %     fr = calc_fr(v(ifile).spikes.times, 'basepath', basepath,...
    %         'graphics', false, 'binsize', 60, 'saveVar', true,...
    %         'smet', 'none', 'winBL', winBL, 'winCalc', [0, Inf]);
    %
    %     fr_timebins('basepath', pwd,...
    %         'forceA', true, 'graphics', true,...
    %         'timebins', timebins, 'saveVar', true, 'sstates', [1, 4]);


    % bins w/ respec to states
    ss = v(ifile).ss.boutTimes([1, 4]);
    cnt = 1; clear bins
    for istate = 1 : 2
        bins{cnt} = ss{istate}(InIntervals(ss{istate}, timebins(1, :)), :);
        bins{cnt + 1} = ss{istate}(InIntervals(ss{istate}, timebins(2, :)), :);
        cnt = cnt + 2;
    end
    bins{5} = timebins(1, :);
    bins{6} = timebins(2, :);
    bintxt = {'AW-BSL'; 'AW-KET'; 'NREM-BSL'; 'NREM-KET'; 'BSL'; 'KET'};

    % spike timing metrics
%     st = spktimes_metrics('spikes', v(ifile).spikes, 'sunits', [],...
%         'bins', bins, 'forceA', true, 'saveVar', true, 'fullA', false);

    % brst (mea)
    brst = spktimes_meaBrst(v(ifile).spikes.times, 'binsize', [], 'isiThr', 0.05,...
        'minSpks', 2, 'saveVar', true, 'force', true, 'bins', bins);
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
for ifile = 1 : nfiles
    frBins = catfields(v(ifile).frBins, 'catdef', 'addim', 'force', false);
    stateMfr = [stateMfr; frBins.states.mfr(:, sstates, :)];
end

% mfr
fr = catfields([v.fr], 'catdef', 'long', 'force', false);
mfr_bl = fr.mfr;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% investigate burstiness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iunit = 1;

brstVar = ["rateNorm"; "rate"; "freq"; "spkprct"; "brstDur"; "ibi"; "detect"; "shortPrct"];    
brst = catfields([v(:).brst], 'catdef', 'long');
for ivar = 1 : length(brstVar)
    fh = figure;
    ydata = brst.(brstVar{ivar})(:, units.clean(iunit, :))';
    plot_boxMean('dataMat', ydata, 'clr', 'k')
    xticklabels(bintxt)
    ylabel(brstVar{ivar})
end

spkVar = ["royer"; "lidor"; "doublets"; "mizuseki"; "lvr"];    
st = catfields([v(:).st], 'catdef', 'long');
for ivar = 1 : length(spkVar)
    fh = figure;
    ydata = st.(spkVar{ivar})(:, units.clean(iunit, :))';
    plot_boxMean('dataMat', ydata, 'clr', 'k')
    xticklabels(bintxt)
    ylabel(spkVar{ivar})
end

fr.mfr(units.clean(iunit, :));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mfr states - box plot of gain ratio divided to sextiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% note that this separates according to percentiles in the input which is
% not necassarily the baseline
plot_FRstates_sextiles('stateMfr', squeeze(stateMfr(:, :, 1))', 'units', units.clean,...
    'ntiles', 2, 'saveFig', false)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% percent change and state ratio after injection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% params
istate = 1;
iunit = 1;
ntiles = 2;

% data
mfr_bl = squeeze(stateMfr(units.clean(iunit, :), 1, 1));
mfr_chg = squeeze(stateMfr(units.clean(iunit, :), istate, 2)) ./...
    squeeze(stateMfr(units.clean(iunit, :), istate, 1)) * 100;

% divide units by prctiles
tiles = prctile(mfr_bl, [1 / ntiles : 1 / ntiles : 1 - 1 / ntiles] * 100);
tiles = [0, tiles, Inf];

clear mfrMat chgMat unitsTile
for itile = 1 : ntiles
    unitsTile(:, itile) = mfr_bl >= tiles(itile) & mfr_bl < tiles(itile + 1);
    chgMat{itile} = mfr_chg(unitsTile(:, itile));
end
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
ibin = 1;
iunit = 1;
bedges = logspace(-1, 0.2, 10);

ttltxt = {'BSL'; 'KET'};

fh = figure;
th = tiledlayout(1, length(sstates), 'TileSpacing', 'Compact');

for ibin = [1, 2]
    axh = nexttile;
    for istate = [1, 2]
        %     ydata = log10(stateMfr(unitsTile(:, 1), istate, ibin));

        ydata = (stateMfr(units.clean(iunit, :), istate, ibin));
        ydata(ydata == 0) = [];
        logstats(ibin, istate, :) = lognfit(ydata);
        mstats(ibin, istate) = mean(ydata, 1, 'omitnan');
        sstats(ibin, istate) = std(ydata, 1, 'omitnan');
        ydata = log10(ydata);

        hh = histfit(ydata, 15, 'normal', 'Normalization', 'probability');
        hh(1).EdgeColor = 'none';
        hh(1).FaceColor = cfg.colors{sstates(istate)};
        hh(1).FaceAlpha = 0.3;
        hh(2).Color = cfg.colors{sstates(istate)};
        hh(2).LineWidth = 3;
        hold on
        xlim([-3 3])
        ylim([0 30])

    end
    xlabel('log MFR [Hz]')
    ylabel('Probability')
    title(axh, ttltxt{ibin})
end

fh = figure;
ydata = squeeze(stateMfr(units.clean(iunit, :), :, 1));
plot_boxMean('dataMat', ydata, 'clr', 'k')


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
% mfr at different timebins regardelss of states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% params
tbins = [1, 2, 5];
iunit = 2;

% mean of all units
chg_raw = [];
for ifile = 1 : nfiles
    unitsClean = v(ifile).units.clean(iunit, :);

    chg_tmp = nan(sum(unitsClean), length(tbins));
    for ibin = 1 : length(tbins)
        chg_tmp(:, ibin) = v(ifile).frBins(tbins(ibin)).mfr(unitsClean);   
    end
    chg_raw = [chg_raw; chg_tmp];
end
chg_raw'

% mean per session
chg_raw = nan(length(tbins) + 1, nfiles);
chg_norm = nan(length(tbins), nfiles);
for ifile = 1 : nfiles
    unitsClean = v(ifile).units.clean(iunit, :);

    chg_raw(1, ifile) = mean(v(ifile).frBins(1).mfr(unitsClean));
    cnt = 1;
    for ibin = 1 : length(tbins)
        chg_raw(cnt + 1, ifile) = mean(v(ifile).frBins(tbins(ibin)).mfr(unitsClean), 'omitnan');
        chg_norm(cnt, ifile) = mean(v(ifile).frBins(tbins(ibin)).mfr(unitsClean) ./...
            v(ifile).frBins(1).mfr(unitsClean) * 100, 'omitnan');
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

% list of all units per session
iunit = 2; 
if iunit == 1
    unitChar = 'rs';
else
    unitChar = 'fs';
end
prism_units = [];
for ifile = 1 : nfiles
    units_tmp = ones(sum(v(ifile).units.(unitChar)), 1) * ifile;
    prism_units = [prism_units; units_tmp];
end
sum(prism_units ~= 0)
prism_units(~units.clean(iunit, units.(unitChar))) = 0;

% grab fr
prism_data = frMat(units.(unitChar), :)';
prism_injTime = max(timeIdx) / 60;
prism_tstamps = [1 : length(frMat)]' / 60 - prism_injTime;

% mfr during baseline
mean(mean(prism_data(prism_tstamps < 0, prism_units > 0), 1, 'omitnan'))

% params
istate = [1, 4];
ibin = [1, 2];
iunit = 1;
ntiles = 2;

% cat units
units = catfields([v.units], 'catdef', 'long', 'force', false);

% cat mfr in states. note ratio vs. gain
stateMfr = []; gainf = []; sratio = [];
for ifile = 1 : nfiles
    frBins = catfields(v(ifile).frBins, 'catdef', 'addim', 'force', false);
    stateMfr = [stateMfr; frBins.states.mfr];
    gainf = [gainf; squeeze(frBins.states.gain(4, :, :))];      
    sratio = [sratio; squeeze(frBins.states.ratio(4, 1, :, :))];
end
stateMfr = stateMfr(units.clean(iunit, :), istate, ibin);
gainf = gainf(units.clean(iunit, :), ibin);
sratio = sratio(units.clean(iunit, :), ibin);

% data
bsl_wake = squeeze(stateMfr(:, 1, 1));
inj_wake = squeeze(stateMfr(:, 1, 2));
bsl_nrem = squeeze(stateMfr(:, 2, 1));
inj_nrem = squeeze(stateMfr(:, 2, 2));

% divide units by prctiles (median) according to wake 
tiles = prctile(bsl_wake, [1 / ntiles : 1 / ntiles : 1 - 1 / ntiles] * 100);
tiles = [0, tiles, Inf];
unitsTile = false(ntiles, length(bsl_wake));
for itile = 1 : ntiles
    unitsTile(itile, :) = bsl_wake >= tiles(itile) & bsl_wake < tiles(itile + 1);
end

% state mfr
stateMat = [bsl_wake, inj_wake, bsl_nrem, inj_nrem];

% precent change 
for itile = 1 : ntiles
    sunits = unitsTile(itile, :);
    chgCell{itile} = inj_wake(sunits) ./ bsl_wake(sunits);
    chgCell{itile + 2} = inj_nrem(sunits) ./ bsl_nrem(sunits);
end
chgMat = cell2nanmat(chgCell, 2) * 100;

% corrected gain factor
clear gaint
gaint(:, 1) = (bsl_nrem - bsl_wake)' ./ max([bsl_nrem, bsl_wake]');
gaint(:, 2) = (inj_nrem - inj_wake)' ./ max([inj_nrem, inj_wake]'); 


% gain factor
clear gainCell stateCell
for itile = 1 : ntiles
    sunits = unitsTile(itile, :);
    gainCell{itile} = gainf(sunits, 1);
    gainCell{itile + 2} = gainf(sunits, 2);
    stateCell{itile} = sratio(sunits, 1);
    stateCell{itile + 2} = sratio(sunits, 2);

    gainCell{itile} = gaint(sunits, 1);
    gainCell{itile + 2} = gaint(sunits, 2);

end
gainMat = cell2nanmat(gainCell, 2) * 100;
stateMat = cell2nanmat(stateCell, 2) * 100;


% for 2-way anova
[chgMat(:, 1)', chgMat(:, 2)'; chgMat(:, 3)', chgMat(:, 4)']

[gainMat(:, 1)', gainMat(:, 2)'; gainMat(:, 3)', gainMat(:, 4)']

[stateMat(:, 1)', stateMat(:, 2)'; stateMat(:, 3)', stateMat(:, 4)']


% kernel plot -------------------------------------------------------------
data = [bsl_wake, inj_wake ./ bsl_wake * 100];
data(31, :) = [];
fh = figure;
th = tiledlayout(1, 2, 'TileSpacing', 'compact');

% define grid for evaluating the kernel density estimate
% xgrid = linspace(-1, max(data(:, 1)));
% ygrid = linspace(-1, max(data(:, 2)));
xgrid = linspace(-1, 7);
ygrid = linspace(-50, 200);
[x1, y1] = meshgrid(xgrid, ygrid);
xi = [x1(:) y1(:)];

% estimate the kernel density
[f, ep, bw] = ksdensity(data, xi, 'function', 'pdf');

% format in matrix
X = reshape(ep(:, 1), length(xgrid), length(ygrid));
Y = reshape(ep(:, 2), length(xgrid), length(ygrid));
Z = reshape(f, length(xgrid), length(ygrid));

% plot
axh = nexttile;
ph = contourf(X, Y, Z, 8);
cr = linspace(210 / 255, 10 / 255);
map = [240 * ones(1, 100) / 255; cr; cr]';
colormap(axh, map)
% xlim([0 7])
% ylim([0 200])
xlabel('MFR BSL (Hz)')
ylabel('MFR KET / MFR BSL (%)')
title('AW')

% NREM
data = [bsl_nrem, inj_nrem ./ bsl_nrem * 100];
[f, ep, ~] = ksdensity(data, xi, 'function', 'pdf');
Z = reshape(f, length(xgrid), length(ygrid));

% plot
axh = nexttile;
ph = contourf(X, Y, Z, 10);
cr = linspace(180 / 255, 50 / 255);
map = [cr; cr; 200 * ones(1, 100) / 255]';
colormap(axh, map)
% xlim([0 7])
% ylim([0 200])
xlabel('MFR BSL (Hz)')
ylabel('MFR KET / MFR BSL (%)')
title('NREM')