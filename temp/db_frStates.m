% fr_states

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data base
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% aCSF
basepaths = [
    {'F:\Data\lh93\lh93_210811_102035'},...
    {'F:\Data\lh95\lh95_210824_083300'},...
    {'F:\Data\lh96\lh96_211201_070100'},...
    {'F:\Data\lh99\lh99_211218_090630'},...
    {'F:\Data\lh100\lh100_220405_100406'},...
    {'F:\Data\lh107\lh107_220509_095738'},...
    ];
%     {'F:\Data\Processed\lh96\lh96_211124_073800'},...
%     {'F:\Data\lh93\lh93_210811_102035'},...
%     {'F:\Data\lh99\lh99_211218_090630'},...

% local ket
basepaths = [
    {'F:\Data\lh93\lh93_210813_110609'},...
    {'F:\Data\lh95\lh95_210825_080400'},...
    {'F:\Data\lh96\lh96_211204_084200'},...
    {'F:\Data\lh100\lh100_220403_100052'},...
    {'F:\Data\lh107\lh107_220501_102641'},...
    ];
% {'F:\Data\lh99\lh99_220119_090035'},...
%     {'F:\Data\lh96\lh96_211126_072000'},...
%     {'F:\Data\lh96\lh96_211202_070500'},...

% ip ket 10 mg/kg
basepaths = [
    {'F:\Data\lh96\lh96_211206_070400'},...
    {'F:\Data\lh98\lh98_211224_084528'},...
    {'F:\Data\lh106\lh106_220512_102302'},...
    {'F:\Data\lh107\lh107_220512_102302'},...
    ];
% lh99_211224_084528
% F:\Data\lh81\lh81_210204_190000   inj 0zt

% ip ket 60 mg/kg
basepaths = [
    {'F:\Data\lh106\lh106_220513_104446'},...
    {'F:\Data\lh107\lh107_220513_104446'},...
    ];

% baclofen
basepaths = [
    {'F:\Data\Processed\lh96\lh96_211207_071500'},...
    {'K:\Data\lh99\lh99_211220_091903'},...
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reanalyze something
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load vars from each session
forceL = true;
varsFile = ["fr"; "fr_bins"; "spikes"; "datInfo"; "session"; "units"];
varsName = ["fr"; "frBins"; "spikes"; "datInfo"; "session"; "units"];
vload = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);
nsessions = length(basepaths);

cell_metrics = CellExplorer('basepaths', basepaths);

% itterate
for isession = 1 : nsessions
    basepath = basepaths{isession};
    cd(basepath)
    [mousename, basename] = fileparts(basepath);
    [~, mousename] = fileparts(mousename);

    % number of units per spike group
    plot_nunits_session('basepath', basepath, 'frBoundries', [])

    % select specific units
    units = selectUnits('basepath', basepath, 'grp', [2 : 7], 'saveVar', true,...
        'forceA', true, 'frBoundries', [0 Inf; 0 Inf],...
        'spikes', vload(isession).spikes);
    units = vload(isession).units;
    unitsClean = units.idx & units.gini' & units.stable' & units.mfrBL' &...
        units.su' & units.cnt';
    
    % plot fr vs. time
    plot_FRtime_session('basepath', basepath,...
        'muFlag', false, 'saveFig', false,...
        'dataType', 'strd', 'units', unitsClean)

    % state ratio according to mfr percetiles
    frBins = catfields(vload(isession).frBins, 'catdef', 'addim', 'force', false);
    stateMfr = frBins.states.mfr(:, [1, 4], 1);
    plot_FRstates_sextiles('stateMfr', stateMfr', 'units', unitsClean,...
        'ntiles', 2, 'saveFig', false)

%     timebins = v(isession).session.general.timebins;
%     timepnt = vload(isession).session.general.timepnt;
%     timebins = [1, timepnt; timepnt, timepnt + 3 * 60 * 60;...
%         timepnt + 3 * 60 * 60, timepnt + 6 * 60 * 60;
%         timepnt + 6 * 60 * 60, timepnt + 12 * 60 * 60;
%         timepnt + 12 * 60 * 60, Inf];
%     winBL = [0 timepnt];
% 
%     fr = firingRate(vload(isession).spikes.times, 'basepath', basepath,...
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

% load vars from each session
varsFile = ["session"; "units"; "fr_bins"; "fr"; "sleep_states"];
varsName = ["session"; "units"; "frBins"; "fr"; "ss"];
v = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);
nsessions = length(basepaths);

% select states
sstates = [1, 4];

% cat and select units
units = catfields([v.units], 'catdef', 'long', 'force', false);
unitsClean = units.idx & units.gini' & units.stable' & units.mfrBL' &...
    units.su' & units.cnt';
unitsClean = unitsClean(1, :) | unitsClean(2, :);

% concate mfr in states
stateMfr = [];
for isession = 1 : nsessions
    frBins = catfields(v(isession).frBins, 'catdef', 'addim', 'force', false);
    stateMfr = [stateMfr; frBins.states.mfr(:, sstates, :)];
end

fr = catfields([v.fr], 'catdef', 'long', 'force', false);
fr.mfr;
% squeeze(stateMfr(unitsClean, :, 2))
dataMat = fr.states.mfr(:, sstates);
plot_FRstates_sextiles('stateMfr', dataMat', 'units', unitsClean,...
    'ntiles', 2, 'saveFig', false)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics - box plot of mfr state ratio divided to sextiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_FRstates_sextiles('stateMfr', squeeze(stateMfr(:, :, 1))', 'units', unitsClean,...
    'ntiles', 2, 'saveFig', false)

istate = 1;
% normChange = (squeeze(stateMfr(unitsClean, istate, 2)) -...
%     squeeze(stateMfr(unitsClean, istate, 1))) ./...
%     (squeeze(stateMfr(unitsClean, istate, 1)) +...
%     squeeze(stateMfr(unitsClean, istate, 2)));


xval = squeeze(stateMfr(unitsClean, istate, 1));
yval = (squeeze(stateMfr(unitsClean, istate, 2)) ./...
    squeeze(stateMfr(unitsClean, istate, 1)) * 100);

infidx = isinf(log10(xval));
mdl = fitlm(log10(xval(~infidx)), yval(~infidx),...
    'Exclude', xval(~infidx) < 0.01);

% mdl = fitlm(yval(~infidx), log10(xval(~infidx)),...
%     'Exclude', xval(~infidx) < 0);

fh = figure;
plot(mdl)


fh = figure; 
plot(xval, yval, '*')
set(gca, 'XScale', 'log')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% percent change between bins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cat and select units
units = catfields([v.units], 'catdef', 'long', 'force', false);
unitsClean = units.idx & units.gini' & units.stable' & units.mfrBL' &...
    units.su' & units.cnt';
unitsClean = unitsClean(1, :) | unitsClean(2, :);

fh = figure;
th = tiledlayout(2, length(sstates), 'TileSpacing', 'Compact');
for istate = 1 : length(sstates)
    axh(istate) = nexttile;
    dataMat = squeeze(stateMfr(unitsClean, istate, :));
    plot(dataMat')
    hold on
    plot(mean(dataMat), 'LineWidth', 2, 'Color', 'k')
    set(gca, 'yscale', 'log')
    title(cfg.names{sstates(istate)})
    ylabel('MFR')
    xlabel('Bin #')
end
nexttile;
fr_diff = squeeze(stateMfr(unitsClean, :, 2)) ./...
squeeze(stateMfr(unitsClean, :, 1)) * 100;
dataMat = fr_diff;
plot_boxMean('dataMat', dataMat, 'clr', 'k', 'allPnts', false)
ylim([prctile(dataMat(:), 3), prctile(dataMat(:), 97)])
ylabel('Change in MFR [%]')
xticklabels(cfg.names(sstates))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics - psd in states according to timebins
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% analyze
for isession = 1 : nsessions
    cd(basepaths{isession});
    session = v(isession).session;
    [timebins, timepnt] = metaInfo_timebins('reqPnt', 5.5 * 60 * 60,...
        'nbins', 8);
    psdBins = psd_states_timebins('basepath', pwd,...
        'chEeg', [], 'forceA', true, 'graphics', true,...
        'timebins', timebins, 'saveVar', true, 'sstates', [1, 4, 5]);
end

% load
varsFile = ["psdBins"; "session"];
varsName = ["psdBins"; "session"];
v = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);

% 2prism
freq = psdBins.info.freq;
ibin = 2;
istate = 3;
psdState = nan(length(freq), nsessions);
for isession = 1 : nsessions
    psdtmp = squeeze(v(isession).psdBins.psdLfp(ibin, istate, :));
    psdState(:, isession) = psdtmp / sum(psdtmp);
    psdState(:, isession) = psdtmp;
end

% time in states
nbins = 4;
sstates = [1, 4];
totDur = nan(nbins, length(sstates), nsessions);
epLen = cell(length(sstates), nsessions);
for isession = 1 : nsessions
    basepath = basepaths{isession};
    cd(basepath)
    
    [tempDur, epLen(:, isession), timebins] =...
        as_plotZT('nbins', nbins, 'sstates', sstates, 'ss', ss,...
        'graphics', false);
    
    binLen = diff(timebins');
    totDur(:, :, isession) = (tempDur ./ binLen') * 100;
end   

fh = figure;
th = tiledlayout(2, 2, 'TileSpacing', 'Compact');
axh = nexttile;
dataVec = mean(totDur, 3, 'omitnan');
plot([1 : nbins], dataVec(:, istate),...
    'Color', cfg.colors{sstates(istate)}, 'LineWidth', 2)
ylabel('State duration [%]')
ylim([0 100])
ax = gca;
set(ax.YAxis, 'color', cfg.colors{sstates(istate)})


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics - box plot of mfr across states divided by median
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setMatlabGraphics(false)
stateIdx = [1, 4];

fh = figure;
for iunit = 1 : 2
    
    % firing rate median
    medFr = median(mfr(units(iunit, :)));
    yLimit = ceil([0 prctile(stateMfr(:, units(iunit, :)), 98 ,'all')]);
    
    % high firing units
    subplot(2, 2, iunit + iunit - 1)
    medUnits = units(iunit, :) & (mfr > medFr)';
    dataMat = stateMfr(:, medUnits);
    plot_boxMean('dataMat', dataMat', 'clr', 'k')
    bh = findobj(gca, 'Tag', 'Box');
    bh = flipud(bh);
    for ibox = 1 : length(bh)
        patch(get(bh(ibox), 'XData'), get(bh(ibox), 'YData'),...
            cfg.colors{stateIdx(ibox)}, 'FaceAlpha', 0.5)
    end
    ylabel('Firing Rate [Hz]')
    xticklabels(cfg.names)
    title(sprintf('High MFR %s = %d', unitChar{iunit}, sum(medUnits)))
    ylim(yLimit)
    
    % low firing units
    subplot(2, 2, iunit + 1 + iunit - 1)
    medUnits = units(iunit, :) & (mfr < medFr)';
    dataMat = stateMfr(:, medUnits);
    plot_boxMean('dataMat', dataMat', 'clr', 'k')
    bh = findobj(gca, 'Tag', 'Box');
    bh = flipud(bh);
    for ibox = 1 : length(bh)
        patch(get(bh(ibox), 'XData'), get(bh(ibox), 'YData'),...
            cfg.colors{stateIdx(ibox)}, 'FaceAlpha', 0.5)
    end
    ylabel('Firing Rate [Hz]')
    xticklabels(cfg.names)
    title(sprintf('Low MFR %s = %d', unitChar{iunit}, sum(medUnits)))
    ylim(yLimit)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics - gini ceoff vs. mfr at baseline colored by stability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load concatenated data from basepaths
frCat = catFrSessions('basepaths', basepaths, 'binidx', 3);

% assing vars
structvars(frCat);
stateMfr = frCat.stateMfr;         
stateRat = frCat.stateRat;         
units = frCat.units;               
unitsGini = frCat.unitsGini;       
unitsMfr = frCat.unitsMfr;         
mfr = frCat.mfr;                   
% giniCoeff = frCat.giniCoeff;       
stableIdx = frCat.stableIdx;       
sessionIdx = frCat.sessionIdx;     
tstamps = frCat.tstamps;     

setMatlabGraphics(false)

fh = figure;
for iunit = 1 : 2
    subplot(2, 2, iunit)
    unitIdx = unitsOrig(iunit, :) & stableIdx';
    scatter(mfr(unitIdx), giniCoeff(unitIdx), 300, '.b')
    hold on
    unitIdx = unitsOrig(iunit, :) & ~stableIdx';
    scatter(mfr(unitIdx), giniCoeff(unitIdx), 300, '.r')
    set(gca, 'XScale', 'log')
    ylim([0 1])
    plot(xlim, [0.5, 0.5], '--k')
    subtitle(unitChar{iunit})
    xlabel('Firing Rate [Hz')
    ylabel('Gini Coefficient')
    
    subplot(2, 2, iunit + 2)
    unitIdx = unitsOrig(iunit, :) & stableIdx';
    scatter(stateRat(unitIdx), giniCoeff(unitIdx), 300, '.b')
    hold on
    unitIdx = unitsOrig(iunit, :) & ~stableIdx';
    scatter(stateRat(unitIdx), giniCoeff(unitIdx), 300, '.r')
    ylim([0 1])
    plot(xlim, [0.5, 0.5], '--k')
    xlabel({sprintf('%s - %s /', cfg.names{4}, cfg.names{1}),...
        sprintf('%s + %s', cfg.names{4}, cfg.names{1})})
    ylabel('Gini Coefficient')
end

fh = figure;
unitsGini = find((unitsOrig(1, :) | unitsOrig(2, :)) & giniIdx');
for igini = 1 : length(unitsGini)
    subplot(ceil(length(unitsGini) / 2), 2, igini)
    plot([1 : length(frMat)] / 60, frMat(unitsGini(igini), :))
    subtitle(sprintf('Unit #%d, Gini = %.2f',...
        unitsGini(igini), giniCoeff(unitsGini(igini))))
end

fh = figure;
unitsGini = find((unitsOrig(1, :) | unitsOrig(2, :)) & ~giniIdx');
unitsGini = unitsGini(randperm(length(unitsGini), 16));
for igini = 1 : length(unitsGini)
    subplot(ceil(length(unitsGini) / 2), 2, igini)
    plot([1 : length(frMat)] / 60, frMat(unitsGini(igini), :))
    subtitle(sprintf('Unit #%d, Gini = %.2f',...
        unitsGini(igini), giniCoeff(unitsGini(igini))))
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to prism
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% firing rate vs. time across all units, aligned to point of injection
units = catfields([v.units], 'catdef', 'long', 'force', false);
unitsClean = units.idx & units.gini' & units.stable' & units.mfrBL' &...
    units.su' & units.cnt';

% firing rate vs. time across all units, aligned to point of injection. can
% extract strd or normalized data
[frMat, timeIdx] = alignFR2pnt('basepaths', basepaths, 'dataType', 'strd');

% grab fr
unitNo = 1;
prism_data = frMat(unitsClean(unitNo, :), :)';
prism_tstamps = [1 : length(frMat)]' / 60;
prism_injTime = max(timeIdx) / 60;

% list of units per session
prism_rs = []; prism_fs = [];
for isession = 1 : nsessions
    units = v(isession).units;
    unitsClean = units.idx & units.gini' & units.stable' & units.mfrBL' &...
    units.su' & units.cnt';
    prism_rs = [prism_rs; ones(sum(unitsClean(1, :)), 1) * isession];
    prism_fs = [prism_fs; ones(sum(unitsClean(2, :)), 1) * isession];
end

% select units based on mfr
fr = catfields([v.fr], 'catdef', 'long', 'force', false);
units_mfr = fr.mfr < 3;
unitsClean = unitsClean & units_mfr';

% single session
dataType = 'norm';
data = fr.(dataType)(units.idx(unitNo, :), :)';
data(~isfinite(data)) = nan;
excludeIdx = ~units.gini' | ~units.mfrBL';
prismIdx = excludeIdx(units.idx(unitNo, :));
tstamps = fr.tstamps / 60 / 60;