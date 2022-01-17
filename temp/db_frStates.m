% fr_states

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepaths = [
    {'K:\Data\lh95\lh95_210824_083300'},...
    {'F:\Data\Processed\lh96\lh96_211201_070100'},...
    {'K:\Data\lh99\lh99_211218_090630'},...
    {'G:\Data\lh93\lh93_210811_102035'},...
    ];
% {'I:\lh86\lh86_210304_070700'},...
% {'I:\lh87\lh87_210523_100607'},...
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

% params
stateIdx = [1 : 5];
frBoundries = [0.1 Inf; 0.1 Inf];
stableFlag = false;
cfg = as_loadConfig();
nstates = cfg.nstates;

% initialize
stateMfr = [];
units = [];
mfr = [];
giniCoeff = [];
stableIdx = [];
for isession = 1 : nsessions
    basepath = basepaths{isession};
    cd(basepath)    
    
    % mean firing rate
    mfr = [mfr; v(isession).fr.mfr];
    
    % mean firing rate for each unit in states
    stateMfr_tmp = cellfun(@(x) mean(x, 2), v(isession).fr.states.fr, 'uni', false);
    stateMfr = [stateMfr, cell2mat(stateMfr_tmp(stateIdx))'];
    
    % select specific units
    if isfield(v(isession).spikes, 'units')
        grp = v(isession).spikes.units.grp;
    else
        grp = [];
    end
    units = logical([units, selectUnits('basepath', basepath, 'grp', grp,...
        'forceA', true, 'frBoundries', frBoundries,...
        'stableFlag', stableFlag)]);
    
    % gini and stable
    giniCoeff = [giniCoeff; v(isession).fr.gini_unit];
    stableIdx = logical([stableIdx; v(isession).fr.stable]);
end

% firing rate vs. time across all units, aligned to point of injection
% [frMat, timeIdx] = alignFR2pnt('basepaths', basepaths);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics - mfr wake vs. nrem per unit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fh = figure;
eqLine = [min(stateMfr(:)), max(stateMfr(:))];
eqLine = [0.01, 100];

% RS
subplot(1, 2, 1)
x = stateMfr(1, units(1, :));
y = stateMfr(4, units(1, :));
fitCoeff = polyfit(log10(x), log10(y), 1);
y2 = 10 .^ [polyval(fitCoeff, log10(x))];
plot(x, y, '.b')
hold on
plot(eqLine, eqLine, 'k')
plot(x, y2, '--b')
set(gca, 'YScale', 'log', 'XScale', 'log')
xlabel('WAKE firing rate [Hz]')
ylabel('NREM firing rate [Hz]')
subtitle('RS')

% FS
subplot(1, 2, 2)
x = stateMfr(1, units(2, :));
y = stateMfr(4, units(2, :));
fitCoeff = polyfit(log10(x), log10(y), 1);
y2 = 10 .^ [polyval(fitCoeff, log10(x))];
plot(x, y, '.r')
hold on
plot(eqLine, eqLine, 'k')
plot(x, y2, '--r')
set(gca, 'YScale', 'log', 'XScale', 'log')
xlabel('WAKE firing rate [Hz]')
ylabel('NREM firing rate [Hz]')
subtitle('FS')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics - box plot of mfr across states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unitType = {'RS', 'FS'};

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
% graphics - gini ceoff vs. mfr at basline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fh = figure;
for iunit = 1 : 2
    subplot(1, 2, iunit)
    unitIdx = units(iunit, :) & stableIdx';
    scatter(mfr(unitIdx), giniCoeff(unitIdx), 200, '.b')
    hold on
    unitIdx = units(iunit, :) & ~stableIdx';
    scatter(mfr(unitIdx), giniCoeff(unitIdx), 200, '.r')
    set(gca, 'XScale', 'log')
    ylim([0 1])
    subtitle(unitType{iunit})
    xlabel('Firing Rate [Hz')
    ylabel('Gini Coefficient')
end


