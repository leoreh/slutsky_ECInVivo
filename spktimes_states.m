
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepaths = {
    'F:\Data\lh86\lh86_210304_070700';...
    'F:\Data\lh87\lh87_210523_100607';...
    'F:\Data\lh95\lh95_210824_083300';...
    'F:\Data\lh96\lh96_220121_090213';...
    'F:\Data\lh98\lh98_211224_084528';...
    'F:\Data\lh100\lh100_220405_100406';...
    'F:\Data\lh106\lh106_220512_102302';...
    'F:\Data\lh107\lh107_220509_095738'};
nfiles = length(basepaths);

clear mousenames
for ifile = 1 : nfiles
    mousepath = fileparts(basepaths{ifile});
    [~, mousenames{ifile}] = fileparts(mousepath);
end
mousenames = string(mousenames);

% load data
varsFile = ["spikes"; "datInfo"; "sleep_states"; "fr"; "units";...
    "st_metrics"; "st_brst"; "session"];
varsName = ["spikes"; "datInfo"; "ss"; "fr"; "units";...
    "st"; "brst"; "session"];
v = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName);

sstates = [1, 4, 5];
stateNames = v(1).ss.info.names;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% single session analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% analyze
for ifile = 1 : nfiles
    
    basepath = basepaths{ifile};
    cd(basepath)
    
    % states
    bins = v(ifile).ss.stateEpochs(sstates);
    
    % select specific units
    units = selectUnits('basepath', pwd, 'grp', [1 : 4], 'saveVar', true,...
        'forceA', true, 'frBoundries', [0.05 Inf; 0.05 Inf],...
        'spikes', v(ifile).spikes);
    
    % spike timing metrics
    st = spktimes_metrics('spikes', v(ifile).spikes, 'sunits', [],...
        'bins', bins, 'forceA', true, 'saveVar', true, 'fullA', false);

    % brst (mea)
    brst = spktimes_meaBrst(v(ifile).spikes.times, 'binsize', [], 'isiThr', 0.1,...
        'minSpks', 8, 'saveVar', true, 'force', true, 'bins', bins);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% investigate burst criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% analyze

bfile = fullfile('F:\Data', 'brst_sessions.mat');
minSpks = [2 : 2 : 10];
isiThr = [0.02, 0.05, 0.1];

load(bfile, 'b')
% clear b
% for ifile = 1 : nfiles
%     basepath = basepaths{ifile};
%     cd(basepath)
%     bins = v(ifile).ss.stateEpochs(sstates);
%     for ispk = 1 : length(minSpks)
%         for ithr = 1 : length(isiThr)
%             b(ifile, ispk, ithr) = spktimes_meaBrst(v(ifile).spikes.times,...
%                 'binsize', [], 'isiThr', isiThr(ithr),...
%                 'minSpks', minSpks(ispk), 'saveVar', false,...
%                 'force', true, 'bins', bins);
%         end
%     end
% end
% b = rmfield(b, 'all');
% save(bfile, 'b', '-v7.3');

% -------------------------------------------------------------------------
% organize across sessions and calc gain factor

% cat units from sessions
sfiles = [1 : nfiles];
clear brst
for ispk = 1 : length(minSpks)
    for ithr = 1 : length(isiThr)
        brst(ispk, ithr) = catfields([b(sfiles, ispk, ithr)], 'catdef', 'long');
    end
end

% calc gain factor
brstVar = ["rateNorm"; "rate"; "freq"; "spkprct"; "brstDur"];    
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

units = catfields([v(sfiles).units], 'catdef', 'long');
unitIdx = units.rs;

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
% plot mfr and burst frequency in states 

% cat
fr = catfields([v(:).fr], 'catdef', 'symmetric');
st = catfields([v(:).st], 'catdef', 'long');
units = catfields([v(:).units], 'catdef', 'long');
brst = catfields([v(:).brst], 'catdef', 'long');
brst = catfields([b(:, 1, 1)], 'catdef', 'long');

% select units
unitIdx = units.rs;

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
ydata = brst.rate(:, unitIdx)';
plot_boxMean(ydata, 'clr', 'k', 'allPnts', false)
set(gca, 'yscale', 'log')
ylabel('MFR [Hz]')
xticklabels(stateNames(sstates))

% select only high-firing units
unitHiFr = fr.states.mfr(:, 4) > median(fr.states.mfr(unitIdx, 4), 'omitnan');
ydata = fr.states.mfr(unitIdx & unitHiFr', [1, 4, 5]);

% select only high-mbr units
unitHiBrst = brst.rate(1, :) > median(brst.rate(1, unitIdx), 'omitnan');
ydata = fr.states.mfr(unitIdx & unitHiBrst, [1, 4, 5]);

% burst index and lvr across states
fh = figure;
th = tiledlayout(1, 2, 'TileSpacing', 'Compact');
axh = nexttile;
ydata = st.lidor(:, unitIdx)';
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

% various vars across states
ydata = brst.spkprct(:, unitIdx)';
ydata = brst.nspks(:, unitIdx)';

% mean per session (data only, from v struct)
ivar = 2;
for ifile = 1 : nfiles
    mbr(ifile, :) = mean(v(ifile).brst.(brstVar{ivar})(:, v(ifile).units.rs), 2, 'omitnan');
    mfr(ifile, :) = mean(v(ifile).fr.states.mfr(v(ifile).units.rs, [1, 4, 5]), 1, 'omitnan');
end

% mean per session (data only, from b struct)
ivar = 2;
ispk = 2; ithr = 2;
for ifile = 1 : nfiles
    mbr(ifile, :) = mean(b(ifile, ispk, ithr).(brstVar{ivar})(:, v(ifile).units.rs), 2, 'omitnan');
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
% check percent spikes given different burst params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% analyze
for ifile = 1 : nfiles
    
    basepath = basepaths{ifile};
    cd(basepath)
    
    % states
    bins = v(ifile).ss.stateEpochs(sstates);
    
    % brst (mea)
    btmp(ifile) = spktimes_meaBrst(v(ifile).spikes.times, 'binsize', [], 'isiThr', 0.005,...
        'minSpks', 2, 'saveVar', false, 'force', true, 'bins', bins);
end

btmp = catfields([btmp], 'catdef', 'long');
ydata = btmp.spkprct(:, unitIdx)';


