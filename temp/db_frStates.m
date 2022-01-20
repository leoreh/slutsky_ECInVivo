% fr_states

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% aCSF
basepaths = [
    {'K:\Data\lh95\lh95_210824_083300'},...
    {'F:\Data\Processed\lh96\lh96_211201_070100'},...
    {'K:\Data\lh99\lh99_211218_090630'},...
    {'G:\Data\lh93\lh93_210811_102035'},...
    {'I:\lh96\lh96_211124_073800'},...
    ];

% local ket
basepaths = [
    {'K:\Data\lh95\lh95_210825_080400'},...
    {'I:\lh96\lh96_211126_072000'},...
    {'F:\Data\Processed\lh96\lh96_211202_070500'},...
    {'K:\Data\lh99\lh99_211219_085802'},...
    ];

% baclofen
basepaths = [
    {'F:\Data\Processed\lh96\lh96_211207_071500'},...
    {'K:\Data\lh99\lh99_211220_091903'},...
    ];

nsessions = length(basepaths);

% load vars from each session
forceL = false;
if ~exist('v', 'var') || forceL
    varsFile = ["fr"; "spikes"; "cell_metrics"; "datInfo"; "session"];
    varsName = ["fr"; "spikes"; "cm"; "datInfo"; "session"];
    v = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
        'varsName', varsName);
end

% recalculate firing rate
forceA = false;
if forceA
    for isession = 1 : nsessions
        basepath = basepaths{isession};
        cd(basepath)
        
        timepnt = v(isession).session.general.timepnt;
        winBL = [0 timepnt];
        fr = firingRate(v(isession).spikes.times, 'basepath', basepath,...
            'graphics', false, 'binsize', 60, 'saveVar', true,...
            'smet', 'none', 'winBL', winBL, 'winCalc', [0, Inf]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% concate sessions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% params
stateIdx = [1, 4];
frBoundries = [0 Inf; 0 Inf];
stableFlag = false;
giniFlag = false;
cfg = as_loadConfig();
nstates = cfg.nstates;
unitType = {'RS', 'FS'};
unitClr = {'b', 'r'};

% initialize
stateMfr = [];
stateRat = [];
units = [];
mfr = [];
giniCoeff = [];
stableIdx = [];
sessionIdx = [];

for isession = 1 : nsessions
    basepath = basepaths{isession};
    cd(basepath)    
    
    % mean firing rate
    mfr = [mfr; v(isession).fr.mfr];
    
    % mean firing rate for each unit in states
    stateMfr_tmp = cellfun(@(x) mean(x, 2), v(isession).fr.states.fr, 'uni', false);
    stateMfr = [stateMfr, cell2mat(stateMfr_tmp(stateIdx))'];
    
    % state ratio (NREM - WAKE)
    stateRat = [stateRat; squeeze(v(isession).fr.states.ratio(4, 1, :))];
    
    % select specific units
    if isfield(v(isession).spikes, 'units')
        grp = v(isession).spikes.units.grp;
    else
        grp = [];
    end
    unitsSession = selectUnits('basepath', basepath, 'grp', grp,...
        'forceA', true, 'frBoundries', frBoundries,...
        'stableFlag', stableFlag, 'giniFlag', giniFlag);
    nunits = length(unitsSession);
    units = logical([units, unitsSession]);
    sessionIdx = [sessionIdx; ones(nunits, 1) * isession];
    
    % gini and stable
    giniCoeff = [giniCoeff; v(isession).fr.gini_unit];
    stableIdx = logical([stableIdx; v(isession).fr.stable]);
end

% firing rate vs. time across all units, aligned to point of injection
[frMat, timeIdx] = alignFR2pnt('basepaths', basepaths);

% manualy remove units by criteria
unitsOrig = units;
mfrIdx = mfr < 0.1 | mfr > Inf;
giniIdx = giniCoeff > 0.5;
units(1, :) = units(1, :) & ~giniIdx' & ~mfrIdx';
units(2, :) = units(2, :) & ~giniIdx' & ~mfrIdx';
unitsAny = (units(1, :) | units(2, :));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to prism
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare for prism
fr_tstamps = [1 : length(frMat)] / 60;
max(timeIdx) / 60;
unitIdx = num2str(sessionIdx);
unitIdx(mfrIdx) = repmat('m', sum(mfrIdx), 1);
unitIdx(giniIdx) = repmat('g', sum(giniIdx), 1);

unitNo = 2;
prismData = frMat(units(unitNo, :), :)';
prismIdx = string(unitIdx(units(unitNo, :)))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics - box plot of mfr state ratio divided to sextiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eqLine = 10 .^ [ceil(log10(min(stateMfr(:)))),...
    ceil(log10(max(stateMfr(:))))];
yLimit = [-max(abs(stateRat(unitsAny))), max(abs(stateRat(unitsAny)))];

fh = figure;
for iunit = 1 : 2
    
    subplot(2, 2, iunit)
    
    x = stateMfr(1, units(iunit, :));
    y = stateMfr(2, units(iunit, :));
    
    % sextiles. note calculation of tiles includes all units if they were
    % not excluded during selectUnits.
    ntiles = 6;
    tiles = prctile(x, [1 / ntiles : 1 / ntiles : 1 - 1 / ntiles] * 100);
    tiles = [0, tiles, Inf];

    clear dataMat
    for itile = 1 : ntiles
        unitsTile = units(iunit, :) & stateMfr(1, :) >= tiles(itile) &...
            stateMfr(1, :) < tiles(itile + 1);
        dataMat{itile} = stateRat(unitsTile);
    end
    dataMat = cell2nanmat(dataMat);
    plot_boxMean('dataMat', dataMat, 'clr', unitClr{iunit}, 'allPnts', true)
    hold on
    plot(xlim, [0, 0], '--k')
    ylabel({sprintf('%s - %s /', cfg.names{4}, cfg.names{1}),...
        sprintf('%s + %s', cfg.names{4}, cfg.names{1})})
    xlabel('Sixtiles (by WAKE)')
    ylim(yLimit)
    subtitle(sprintf('%s = %d', unitType{iunit}, sum(units(iunit, :))))
    
    subplot(2, 2, iunit + 2)
    fitCoeff = polyfit(log10(x), log10(y), 1);
    y2 = 10 .^ [polyval(fitCoeff, log10(x))];
    plot(x, y, '.', 'Color', unitClr{iunit}, 'MarkerSize', 10)
    hold on
    plot(eqLine, eqLine, 'k')
    plot(x, y2, '--', 'Color', unitClr{iunit})
    plot([tiles(2 : end - 1); tiles(2 : end - 1)], eqLine, '--g')
    set(gca, 'YScale', 'log', 'XScale', 'log')
    xlabel('WAKE firing rate [Hz]')
    ylabel('NREM firing rate [Hz]')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics - box plot of mfr across states divided by median
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    subtitle(sprintf('High MFR %s = %d', unitType{iunit}, sum(medUnits)))
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
    subtitle(sprintf('Low MFR %s = %d', unitType{iunit}, sum(medUnits)))
    ylim(yLimit)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics - gini ceoff vs. mfr at basline colored by stability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fh = figure;
for iunit = 1 : 2
    subplot(1, 2, iunit)
    unitIdx = unitsOrig(iunit, :) & stableIdx';
    scatter(mfr(unitIdx), giniCoeff(unitIdx), 200, '.b')
    hold on
    unitIdx = unitsOrig(iunit, :) & ~stableIdx';
    scatter(mfr(unitIdx), giniCoeff(unitIdx), 200, '.r')
    set(gca, 'XScale', 'log')
    ylim([0 1])
    subtitle(unitType{iunit})
    xlabel('Firing Rate [Hz')
    ylabel('Gini Coefficient')
end
