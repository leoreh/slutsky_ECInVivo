
% aug 22
% investigating the state regulation of burstiness

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepaths = {
    'F:\Data\Colleagues\RA\tg4_210730_155700';...   % baseline. 24 hr rec
    'F:\Data\Colleagues\RA\hDLX_Gq_WT5\040221_0655_24hr';...    % baseline. 24 hr rec
    'F:\Data\lh52\lh52_200615_084111';...           % baseline. 9 hr rec
    'F:\Data\lh70\lh70_201015_0951';...             % baseline. 8 hr rec
    'F:\Data\lh81\lh81_210206_044300';...           % ip saline. 6 hr rec
    'F:\Data\lh86\lh86_210304_070700';...           % washout ket. 12 hr rec
    'F:\Data\lh87\lh87_210523_100607';...           % op saline
    'F:\Data\lh93\lh93_210811_102035';...           % local nacl
    'F:\Data\lh95\lh95_210824_083300';...           % local nacl
    'F:\Data\lh96\lh96_220120_090157';...           % baseline
    'F:\Data\lh98\lh98_211224_084528';...           % ket ip 10
    'F:\Data\lh100\lh100_220413_111004';...         % op acsf
    'F:\Data\lh106\lh106_220512_102302';...         % ket ip 10. too few units
    'F:\Data\lh107\lh107_220518_091200';...         % op acsf
    'F:\Data\lh111\lh111_220823_094417';...         % baseline
    'F:\Data\lh112\lh112_220828_104358';...         % baseline
    };           
nfiles = length(basepaths);
sfiles = [1 : nfiles];

clear mousenames
for ifile = 1 : nfiles
    mousepath = fileparts(basepaths{ifile});
    [~, mousenames{ifile}] = fileparts(mousepath);
end
mousenames = string(mousenames);

% load data
varsFile = ["spikes"; "datInfo"; "sleep_states"; "fr"; "units";...
    "st_metrics"; "st_brst"; "session"; "sr"];
varsName = ["spikes"; "datInfo"; "ss"; "fr"; "units";...
    "st"; "brst"; "session"; "sr"];
v = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);

sstates = [1, 4, 5];
stateNames = v(1).ss.info.names;
brstVar = ["rateNorm"; "rate"; "freq"; "spkprct"; "brstDur"; "ibi"; "detect"; "shortPrct"];    

units = catfields([v(sfiles).units], 'catdef', 'long');
unitIdx = units.rs;
sunits = find(unitIdx);
nunits = length(sunits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% single session analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% analyze
for ifile = 1 : nfiles
    
    basepath = basepaths{ifile};
    cd(basepath)

    % states
    bins = v(ifile).ss.stateEpochs(sstates);
    bins{end + 1} = [0 Inf];

    % firing rate
%     fr = calc_fr(v(ifile).spikes.times, 'basepath', basepath,...
%         'graphics', false, 'binsize', 60, 'saveVar', true,...
%         'smet', 'none', 'winBL', [0 Inf], 'winCalc', [0, Inf], 'forceA', true);

    % spike timing metrics
    st = spktimes_metrics('spikes', v(ifile).spikes, 'sunits', [],...
        'bins', bins, 'forceA', true, 'saveVar', true, 'fullA', false);
% 
%     % brst (mea)
%     brst = spktimes_meaBrst(v(ifile).spikes.times, 'binsize', [], 'isiThr', 0.1,...
%         'minSpks', 8, 'saveVar', true, 'force', true, 'bins', bins);
% 
%     % select specific units
%     units = selectUnits('basepath', pwd, 'grp', [1 : 4], 'saveVar', true,...
%         'forceA', true, 'frBoundries', [0.05 Inf; 0.05 Inf],...
%         'spikes', v(ifile).spikes);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% investigate burst criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% analyze and combine in single struct

bfile = fullfile('F:\Data', 'brst_sessions.mat');
minSpks = [2 : 6, 8, 10, 15, 20, 30];
isiThr = [0.003, 0.005, 0.008, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1];
% 
% clear b
% for ifile = 1 : nfiles
%     basepath = basepaths{ifile};
%     cd(basepath)
%     bins = v(ifile).ss.stateEpochs(sstates);
%     bins{end + 1} = [0 Inf];
% 
%     for ispk = 1 : length(minSpks)
%         for ithr = 1 : length(isiThr)
%             btmp = spktimes_meaBrst(v(ifile).spikes.times,...
%                 'binsize', [], 'isiThr', isiThr(ithr),...
%                 'minSpks', minSpks(ispk), 'saveVar', false,...
%                 'force', true, 'bins', bins);
%             
%             if isfield(btmp, 'all')
%                 btmp = rmfield(btmp, 'all');
%             end
%             b(ifile, ispk, ithr) = btmp;
%         end
%     end
% end
% 
% % organize struct and save
% binfo.runtime          = datetime(now, 'ConvertFrom', 'datenum');
% binfo.minSpks          = minSpks;
% binfo.isiThr           = isiThr;
% binfo.basepaths        = basepaths;
% save(bfile, 'b', 'binfo');
load(bfile, 'b')

% organize across sessions
clear brst
for ispk = 1 : length(minSpks)
    for ithr = 1 : length(isiThr)
        brst(ispk, ithr) = catfields([b(sfiles, ispk, ithr)], 'catdef', 'long');
    end
end

% calc gain factor
clear gfactor
for ispk = 1 : length(minSpks)
    for ithr = 1 : length(isiThr)
        for ivar = 1 : length(brstVar)
            vec1 = brst(ispk, ithr).(brstVar{ivar})(1, :)';
            vec2 = brst(ispk, ithr).(brstVar{ivar})(2, :)';
            gfactor.(brstVar{ivar})(ispk, ithr, :) = (vec2 - vec1) ./ sum([vec1, vec2]')';
        end
    end
end

% -------------------------------------------------------------------------
% plot gain factor of vars as a function of brst params

% all units combines
for ithr = 1 : length(isiThr)
    fh = figure;
    th = tiledlayout(1, length(brstVar), 'TileSpacing', 'Compact');
    title(th, num2str(isiThr(ithr)))
    
    for ivar = 1 : length(brstVar)
        axh = nexttile;
        dataMat = squeeze(gfactor.(brstVar{ivar})(:, ithr, unitIdx))';
        plot_boxMean(dataMat, 'clr', 'k', 'allPnts', false)
        title(brstVar{ivar})
        ylabel('GainFactor')
        xticklabels(split(num2str(minSpks)))
        xlabel('minSpks')
    end
end

% -------------------------------------------------------------------------
% plot mean value of brst vars as a function of burst params
clear mdata
for ivar = 1 : length(brstVar)
    for ispk = 1 : length(minSpks)
        for ithr = 1 : length(isiThr)
            mdata.(brstVar(ivar))(ispk, ithr) =...
                mean(brst(ispk, ithr).(brstVar(ivar))(4, unitIdx), 'omitnan');
        end
    end
end

ivar = 8;
fh = figure;
plot(minSpks, mdata.(brstVar(ivar)), 'LineWidth', 2)
% set(gca, 'yscale', 'log')
xlabel('Min Spikes')
ylabel(brstVar(ivar))
legend(split(num2str(isiThr)))

% -------------------------------------------------------------------------
% per-session analysis of state ratio and consistancy

% mean per session (data only, from b struct)
% creates 4D mat of nfiles x minSPks, isiThr, states
% mbr = nan(nfiles, length(minSpks), length(isiThr), length(sstates));
% mbrStat = nan(length(minSpks), length(isiThr));
% set(0, 'DefaultFigureVisible', 'off');
% for ispk = 1 : length(minSpks)
%     for ithr = 1 : length(isiThr)
%         for ifile = 1 : nfiles
%             mbr(ifile, ispk, ithr, :) =...
%                 mean(b(ifile, ispk, ithr).rate([1 : length(sstates)],...
%                 v(ifile).units.rs), 2, 'omitnan');
%         end
% 
%         % calc consistancy
%         [~, tbl, stats] = friedman(squeeze(mbr(:, ispk, ithr, :)));
%         if isempty(tbl{2, 5})
%             mbrStat(ispk, ithr) = 0;
%         else
%             mbrStat(ispk, ithr) = tbl{2, 5};
%         end
%     end
% end
% close all
% set(0,'DefaultFigureVisible','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% investigate specific criterion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fr in states by percentiles. ogranize to prism
fitLine = nan(9, nfiles);
for ifile = 1 : nfiles
    [~, fitLine(:, ifile)] = plot_frStateTiles('basepath', basepaths{ifile},...
        'ntiles', 2, 'graphics', false);

end

st.cv(1 : 3, unitIdx);

% -------------------------------------------------------------------------
% plot mfr and burst frequency in states 

% cat
fr = catfields([v(:).fr], 'catdef', 'symmetric');
st = catfields([v(:).st], 'catdef', 'long');
units = catfields([v(:).units], 'catdef', 'long');
brst = catfields([b(:, 5, 8)], 'catdef', 'long');

% plot
fh = figure;
th = tiledlayout(1, 2, 'TileSpacing', 'Compact');
axh = nexttile;
ydata = fr.states.mfr(unitIdx, [1, 4, 5]);
plot_boxMean(ydata, 'clr', 'k', 'allPnts', false)
ylabel('MFR [Hz]')
set(gca, 'yscale', 'log')
xticklabels(stateNames(sstates))

axh = nexttile;
ydata = brst.rate(1 : length(sstates), unitIdx)';
plot_boxMean(ydata, 'clr', 'k', 'allPnts', false)
set(gca, 'yscale', 'log')
ylabel('MFR [Hz]')
xticklabels(stateNames(sstates))

% to prism
squeeze(fr.states.ratio(4, 1, unitIdx))
fr.states.gain(4, unitIdx)

% select only high-firing units
unitHiFr = fr.states.mfr(:, 1) > median(fr.states.mfr(unitIdx, 1), 'omitnan');
ydata = fr.states.mfr(unitIdx & unitHiFr', [1, 4, 5]);

% select only high-mbr units
unitHiBrst = brst.rate(1, :) > median(brst.rate(1, unitIdx), 'omitnan');
ydata = fr.states.mfr(unitIdx & unitHiBrst, [1, 4, 5]);

% burst index and lvr across states
fh = figure;
th = tiledlayout(1, 2, 'TileSpacing', 'Compact');
axh = nexttile;
ydata = st.lvr(:, unitIdx)';
plot_boxMean(ydata, 'clr', 'k', 'allPnts', false)
ylabel('Burst Index (Lidor)')
xticklabels(stateNames(sstates))

axh = nexttile;
ydata = st.lvr(:, unitIdx)';
plot_boxMean(ydata, 'clr', 'k', 'allPnts', false)
ylabel('Firing Irregularity (LvR)')
xticklabels(stateNames(sstates))

% select only bursty units
unitHiBrst = st.lidor(1, :) > median(st.lidor(1, unitIdx), 'omitnan');
ydata = fr.states.mfr(unitIdx & unitHiBrst, [1, 4, 5]);


% mean per session (data only, from v struct)
ivar = 2;
ispk = 6;
ithr = 7;
mbr = nan(nfiles, length(sstates));
mfr = nan(nfiles, length(sstates));
for ifile = 1 : nfiles
    idxUnits = v(ifile).units.rs;

    mbr(ifile, :) = mean(b(ifile, ispk, ithr).(brstVar{ivar})(1 : 3, idxUnits), 2, 'omitnan');
    mfr(ifile, :) = mean(v(ifile).fr.states.mfr(idxUnits, [1, 4, 5]), 1, 'omitnan');
end

% -------------------------------------------------------------------------
% plot correlation of mfr w/ brst freq in AW and NREM
fh = figure;
th = tiledlayout(2, 2, 'TileSpacing', 'Compact');

axh = nexttile;
xdata = fr.states.mfr(unitIdx, 1);
ydata = brst.rate(1, unitIdx);
plot(xdata, ydata, '.', 'MarkerSize', 20)
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
xlabel('MFR in AW [Hz]')
ylabel('MBR in AW (Hz)')
xlim([0.01, 10])
ylim([0.000001, 1])

axh = nexttile;
xdata = fr.states.mfr(unitIdx, 4);
ydata = brst.rate(2, unitIdx);
plot(xdata, ydata, '.', 'MarkerSize', 20)
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
xlabel('MFR in NREM [Hz]')
ylabel('MBR in NREM (Hz)')
xlim([0.01, 10])
ylim([0.000001, 1])

axh = nexttile;
xdata = fr.states.mfr(unitIdx, 1);
ydata = fr.states.mfr(unitIdx, 4);
plot(xdata, ydata, '.', 'MarkerSize', 20)
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
xlabel('MFR in AW [Hz]')
ylabel('MFR in NREM (Hz)')
hold on
plot(xlim, ylim, '--k')
yval = median(fr.states.mfr(unitIdx, 4));
xval = median(fr.states.mfr(unitIdx, 1));
plot([xval xval], ylim, '--r')
plot(xlim, [yval yval], '--r')

axh = nexttile;
xdata = brst.rate(1, unitIdx);
ydata = brst.rate(2, unitIdx);
plot(xdata, ydata, '.', 'MarkerSize', 20)
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
xlabel('MBR in AW [Hz]')
ylabel('MBR in NREM (Hz)')
hold on
plot(xlim, ylim, '--k')
yval = median(brst.rate(2, unitIdx), 'omitnan');
xval = median(brst.rate(1, unitIdx), 'omitnan');
plot([xval xval], ylim, '--r')
plot(xlim, [yval yval], '--r')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyse state differences in burstiness through acg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% lidor: sum of spikes in 2-10 ms normalized to sum in 100-200
stinfo = catfields([v(:).st], 'catdef', 'symmetric');
tstamps = stinfo.info.acg_wide_tstamps(:, 1);

for iunit = 1 : length(fr.mfr)
    for ibin = 1 : length(sstates)
        t1 = find(tstamps > 0.002 & tstamps < 0.01);
        t2 = find(tstamps > 0.1 & tstamps < 0.108);
        burst_temp = sum(st.acg_wide(t1, ibin, iunit));
        bl_temp = sum(st.acg_wide(t2, ibin, iunit));
        st.slowBrst(ibin, iunit) = (burst_temp - bl_temp) ./ (burst_temp + bl_temp);
    end
end

fh = figure;
ydata = st.slowBrst(:, unitIdx)';
plot_boxMean(ydata, 'clr', 'k', 'allPnts', true)
% ylabel('Firing Irregularity (LvR)')
xticklabels(stateNames(sstates))

fh = figure;
th = tiledlayout(2, 2);
axh = nexttile;
plot(fr.mfr(unitIdx), st.lidor(1, unitIdx), '.', 'MarkerSize', 20)
set(gca, 'xscale', 'log')
% set(gca, 'yscale', 'log')

axh = nexttile;
plot(fr.mfr(unitIdx), st.slowBrst(1, unitIdx), '.', 'MarkerSize', 20)
set(gca, 'xscale', 'log')
% set(gca, 'yscale', 'log')

axh = nexttile;
xdata = fr.mfr(unitIdx);
ydata = brst(2, 6).freq(1, unitIdx);
plot(xdata, ydata, '.', 'MarkerSize', 20)
set(gca, 'xscale', 'log')
% set(gca, 'yscale', 'log')

corr(log10(xdata), ydata', 'Rows','complete')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze CV of isi 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% see also  Vyazovskiy et al (Tononi), Neuron, 2009

thr = 20; % [s]
cnt = 1;
isiCell = cell(length(unitIdx), length(sstates));
for ifile = 1 : nfiles
    bins = v(ifile).ss.stateEpochs(sstates);
    for iunit = 1 : length(v(ifile).spikes.times)
        spktimes = v(ifile).spikes.times{iunit};
        for ibin = 1 : length(sstates)

            spks = spktimes(InIntervals(spktimes, bins{ibin}));
            isi = diff(spks);
            isi(isi > thr) = [];
            isi(isi == 0) = [];
            isiCell{cnt, ibin} = (isi);
        end
        cnt = cnt + 1;
    end
end
% mdata = cellfun(@(x) mean(log10(x), 1, 'omitnan'), isiCell(unitIdx, :), 'uni', true);
% sdata = cellfun(@(x) std(log10(x), 1, 'omitnan'), isiCell(unitIdx, :), 'uni', true);
% cvdata = sdata ./ mdata;

% calculate the geometric cv
clear isiHists 
histEdges = logspace(-3, 1.5, 30);
for iunit = 1 : length(unitIdx)
    for ibin = 1 : length(sstates)
        tmp = lognfit(isiCell{iunit, ibin});
        mlog(iunit, ibin) = tmp(1);
        vlog(iunit, ibin) = tmp(2);

        isiHists{ibin}(iunit, :) =...
            histcounts(log10(isiCell{iunit, ibin}), log10(histEdges),...
            'Normalization', 'probability');
    end
end
cvdata = vlog ./ mlog;
gcv = sqrt(exp(vlog .^ 2) - 1);
% gcv = exp(1) .^ log((vlog .* log(10))) - 1;
cvdata = cvdata(unitIdx, :);
gcv = gcv(unitIdx, :);


fh = figure;
plot_boxMean('dataMat', gcv, 'clr', 'k', 'allPnts', false)
set(gca, 'yscale', 'log')

fh = figure;
plot(fr.mfr(unitIdx), gcv, '*')
set(gca, 'xscale', 'log')


% plot isi histograms in states
fh = figure;
th = tiledlayout(1, length(sstates));
ylimit = [0, max(cellfun(@(x) max(x, [], 'all'), isiHists, 'uni', true))];
xLimit = log10([histEdges(1), histEdges(end)]);
for istate = 1 : length(sstates)
    axh = nexttile;
    ydata = isiHists{istate}(unitIdx, :);
    plot(log10(histEdges(2 : end)), ydata)
    hold on
    plot(log10(histEdges(2 : end)), mean(ydata, 1),...
        'k', 'LineWidth', 2)
    title(axh, stateNames{sstates(istate)})
    xlabel('log ISI')
    ylabel('Probability')
    ylim(yLimit)
    xlim(xLimit)
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% per session analysis of metric vs. mfr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear xdata ydata
for ifile = 1 : nfiles
    ydata(ifile) = mean(squeeze(v(ifile).fr.states.ratio(4, 1, v(ifile).units.rs)), 1, 'omitnan');
    xdata(ifile) = sum(v(ifile).units.rs) / (sum(v(ifile).units.rs) + sum(v(ifile).units.fs)) * 100;
end



stVar = 'lvr';

for ifile = 1 : nfiles
    fh = figure;
    th = tiledlayout(2, 2);
    axh = nexttile;
    xdata = v(ifile).st.(stVar)(1, v(ifile).units.rs);
    ydata = v(ifile).st.(stVar)(2, v(ifile).units.rs);
    plot(xdata, ydata, '.', 'MarkerSize', 20)
    maxLimit = max([xlim, ylim]);
    ylim([0 maxLimit])
    xlim([0 maxLimit])
    hold on
    plot(xlim, ylim, '--k')
    title(basepaths{ifile})
    xlabel('CV in AW')
    ylabel('CV in NREM')

    axh = nexttile;
    xdata = v(ifile).fr.states.mfr(v(ifile).units.rs, 1);
    ydata = v(ifile).fr.states.mfr(v(ifile).units.rs, 2);
    plot(xdata, ydata, '.', 'MarkerSize', 20)
    maxLimit = max([xlim, ylim]);
    ylim([0 maxLimit])
    xlim([0 maxLimit])
    hold on
    plot(xlim, ylim, '--k')
    title(basepaths{ifile})
    xlabel('MFR in AW')
    ylabel('MFR in NREM')

    axh = nexttile;
    xdata = v(ifile).fr.mfr(v(ifile).units.rs);
    ydata = v(ifile).st.(stVar)(4, v(ifile).units.rs);
    plot(xdata, ydata, '.', 'MarkerSize', 20)
    title(basepaths{ifile})
    xlabel('MFR')
    ylabel('CV')

    axh = nexttile;
    xdata = v(ifile).fr.mfr(v(ifile).units.rs);
    ydata = squeeze(v(ifile).fr.states.ratio(4, 1, v(ifile).units.rs));
    plot(xdata, ydata, '.', 'MarkerSize', 20)
    ylim([-1 1])
    hold on
    plot(xlim, [0 0], '--k')
    title(basepaths{ifile})
    xlabel('MFR')
    ylabel('State Ratio')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze calcium data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assumes var 'ca' contains MCaR of states from a single session

% subsample few units and check stateRatio
nrep = 8;
nunitsE = 20;
nunitsCa = size(ca, 1);

caShort = nan(nrep, length(sstates));
for irep = 1 : nrep
    uidx = randperm(nunitsCa, nunitsE);
    caShort(irep, :) = mean(ca(uidx, :), 1, 'omitnan');
end


% plot histogram per state
fh = figure;
th = tiledlayout(1, 3);
axh = nexttile;
histogram(log10(ca(:, 1)), 100)
% set(gca, 'xscale', 'log')

% state ratio
(ca(:, 2) - ca(:, 1)) ./ (ca(:, 2) + ca(:, 1));


vec1 = brst.rate(1, unitIdx);
vec2 = brst.rate(2, unitIdx);

(vec2 - vec1) ./ (vec2 + vec1)

fh = figure;
plot(vec2, vec1, '.', 'MarkerSize', 20)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze spike rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% analyze
cnt = 1;
nspksTet = [];
for ifile = 1 : nfiles

    basepath = basepaths{ifile};
    cd(basepath)
    [~, basename] = fileparts(basepath);

    % states
    bins = v(ifile).ss.stateEpochs(sstates);
    bins{end + 1} = [0 Inf];

    load([basename, '.spktimes.mat'])
    fs = v(ifile).session.extracellular.sr;
    spkgrp = v(ifile).session.extracellular.spikeGroups.channels;

    spktimes = cellfun(@(x) (x / fs), spktimes, 'uni', false);

    nspksTet = [nspksTet, cellfun(@numel, spktimes, 'uni', true)];
    
    for igrp = 1 : length(spkgrp)
        grpIdx = v(ifile).spikes.shankID == igrp;
        grpUnits = v(ifile).spikes.times(grpIdx);
        nspksUnit(cnt) = sum(cellfun(@numel, grpUnits, 'uni', true));

        cnt = cnt + 1;
    end

%     sr = calc_fr(spktimes, 'basepath', basepath,...
%         'graphics', false, 'binsize', 60, 'saveVar', 'sr', 'smet', 'none',...
%         'winBL', [0 Inf]);

end


sr = catfields([v(:).sr], 'catdef', 'symmetric');
ydata = sr.states.mfr(:, sstates);

xdata = (nspksUnit ./ nspksTet) * 100;
vec1 = sr.states.mfr(:, 1);
vec2 = sr.states.mfr(:, 4);
ydata = (vec2 - vec1) ./ (vec2 + vec1);


