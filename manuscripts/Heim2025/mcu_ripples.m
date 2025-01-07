

% see also:
% https://www.sciencedirect.com/science/article/pii/S2352289521000357?via%3Dihub


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test the analysis on single session
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% file params
ch = [1 : 4];
recWin = [0 * 60 * 60, 4 * 60 * 60];

basepath = pwd;
cd(basepath)
[~, basename] = fileparts(basepath);

% load data
load([basename, '.sleep_sig.mat'], 'emg')


% get ripples
ripp = getRipples('basepath', basepath, 'rippCh', ch,...
    'emg', emg, 'recWin', recWin, 'saveVar', true,...
    'graphics', true, 'saveVar', true);

% ripple relation to states
ripp = rippleStates(ripp, 'basepath', basepath, 'saveVar', true,...
    'graphics', true);

% ripple relation to spikes
ripp = rippleSpks(ripp, 'basepath', basepath, 'graphics', true,...
    'saveVar', true, 'fullAnalysisFlag', false)

% plot ripples
plot_ripples(ripp, 'basepath', basepath, 'saveFig', true)
plot_rippleSpks(ripp, 'basepath', basepath, 'saveFig', true)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the analysis on all sessions of a mouse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mname = 'lh132';
basepaths = [mcu_sessions(mname)];
nfiles = length(basepaths);

recWin = [0, Inf];

for ifile = 1 : nfiles

    basepath = basepaths{ifile};
    cd(basepath)
    [~, basename] = fileparts(basepath);

%     % load data
%     load([basename, '.sleep_sig.mat'], 'emg')
% 
%     % get ripples
%     ripp = getRipples('basepath', basepath, 'rippCh', ch,...
%         'emg', emg, 'recWin', recWin, 'saveVar', true,...
%         'graphics', true, 'saveVar', true);
% 
%     % ripple relation to states 
    ripp = rippleStates(ripp, 'basepath', basepath, 'saveVar', true,...
        'graphics', true)

    load([basename, '.ripp.mat'], 'ripp')

    % ripple relation to spikes
    ripp = rippleSpks(ripp, 'basepath', basepath, 'graphics', true,...
        'saveVar', true, 'fullAnalysisFlag', false)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize data from all sessions and plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load structs
vars = ["ripp"];
[basepaths, v] = mcu_sessions(mname, vars);
nfiles = length(basepaths);

% cat
ripp = catfields([v(:).ripp], 'addim', true);

% organize ripple rate
rippRate = catfields([ripp(:).rate], 2, true);
rippRate = rippRate.rate(:);
binsize = unique(ripp.info.binsizeRate);
tstamps = [1 : length(rippRate)] / 60;

% organize vars per session
vars = ["maxFreq"; "peakFreq"; "peakAmp"; "dur"]
for ivar = 1 : length(vars)
    vData.(vars(ivar)) = median(squeeze(ripp.(vars(ivar))), 1, 'omitnan')'
end


% graphics

% open figure
setMatlabGraphics(true)
fh = figure;
set(fh, 'WindowState', 'maximized');
tlayout = [ceil(length(vars) / 2) + 1, 2];
th = tiledlayout(tlayout(1), tlayout(2));
th.TileSpacing = 'tight';
th.Padding = 'none';
title(th, mname, 'interpreter', 'none', 'FontSize', 20)
set(fh, 'DefaultAxesFontSize', 16);

% rate of ripples
axh = nexttile(th, 1, [1, 2]); cla; hold on
plot(tstamps, rippRate)
xlabel('Time (hr)')
ylabel('Rate (1/min)')
axis tight

% vars
for ivar = 1 : length(vars)
    axh = nexttile; cla; hold on
    dataMat = squeeze(ripp.(vars(ivar)));
    plot_boxMean('dataMat', dataMat, 'plotType', 'bar', 'axh', axh,...
        'clr', [0.5 0.5 0.5])
    title(axh, vars(ivar))
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare ripples params during baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ifile = 1 : nfiles

    plot_ripples(v(ifile).ripp, 'basepath', basepaths{ifile}, 'saveFig', false)
    plot_rippleSpks(v(ifile).ripp, 'basepath', basepaths{ifile}, 'saveFig', false)

end



varsFile = ["ripp"; "units"];
mname = ["wt_bsl_ripp"; "mcu_bsl"];

clear vData
for igrp = 1 : 2
    
    [basepaths, v] = mcu_sessions(mname{igrp}, varsFile);
    nfiles = length(basepaths);

    % cat
    ripp = catfields([v(:).ripp], 'addim', true);

    % organize vars per session
    varsRipp = ["maxFreq"; "peakFreq"; "peakAmp"; "dur"; "peakPow"]
    for ivar = 1 : length(varsRipp)
        vData(igrp).(varsRipp(ivar)) = median(squeeze(ripp.(varsRipp(ivar))), 1, 'omitnan')';
    end
    vData(igrp).rate = squeeze(mean(ripp.rate.rate, 1, 'omitnan'));
    vData(igrp).rateNrem = cellfun(@(x) mean(x, 'omitnan'), ripp.states.rate(1, 4, :), 'uni', true);

    % calculate spike-ripple gain
    rippGain = cell(2, 1);
    rippMfr = cell(2, 1);
    randMfr = cell(2, 1);
    for ifile = 1 : nfiles

        tmp_ripp = squeeze(mean(mean(v(ifile).ripp.spks.su.rippMap, 2, 'omitnan'), 3, 'omitnan'));
        tmp_rand = squeeze(mean(mean(v(ifile).ripp.spks.su.randMap, 2, 'omitnan'), 3, 'omitnan'));
        tmp_gain = (tmp_ripp - tmp_rand) ./ (tmp_ripp + tmp_rand);

        for iunit = 1 : 2
            unitIdx = v(ifile).units.clean(iunit, :);
            rippGain{iunit} = [rippGain{iunit}; tmp_gain(unitIdx)];
            rippMfr{iunit} = [rippMfr{iunit}; tmp_ripp(unitIdx)];
            randMfr{iunit} = [randMfr{iunit}; tmp_rand(unitIdx)];
            vData(igrp).gain_m(ifile, iunit) = mean(tmp_gain(unitIdx), 'omitnan');
        end
    end
    vData(igrp).gain = cell2padmat(rippGain, 2);
    vData(igrp).rippMfr = cell2padmat(rippMfr, 2);
    vData(igrp).randMfr = cell2padmat(randMfr, 2);
end

vData = catfields([vData(:)], 'addim', true);

% graphics
setMatlabGraphics(true)
fh = figure;
set(fh, 'WindowState', 'maximized');
tlayout = [4, 4];
th = tiledlayout(tlayout(1), tlayout(2));
th.TileSpacing = 'tight';
th.Padding = 'none';
title(th, 'Ripples in Baseline', 'interpreter', 'none', 'FontSize', 20)
set(fh, 'DefaultAxesFontSize', 16);

% vars
varsRipp = ["maxFreq"; "peakFreq"; "peakAmp"; "dur"; "peakPow"; "rate"; "rateNrem"];
for ivar = 1 : length(varsRipp)
    axh = nexttile; cla; hold on
    dataMat = squeeze(vData.(varsRipp(ivar)));
    plot_boxMean('dataMat', dataMat, 'plotType', 'bar', 'axh', axh,...
        'clr', [0.5 0.5 0.5])
    title(axh, varsRipp(ivar))
end

for iunit = 1 : 2
    axh = nexttile; cla; hold on
    dataMat = squeeze(vData.gain(:, iunit, :));
    plot_boxMean('dataMat', dataMat, 'plotType', 'bar', 'axh', axh,...
        'clr', [0.5 0.5 0.5])
    title(axh, 'spike gain')
end

for iunit = 1 : 2
    axh = nexttile; cla; hold on
    dataMat = squeeze(vData.gain_m(:, iunit, :));
    plot_boxMean('dataMat', dataMat, 'plotType', 'bar', 'axh', axh,...
        'clr', [0.5 0.5 0.5])
    title(axh, 'spike gain')
end


for iunit = 1 : 2
    axh = nexttile; cla; hold on
    for igrp = 1 : 2
        xdata = squeeze(vData.randMfr(:, iunit, igrp));
        ydata = squeeze(vData.rippMfr(:, iunit, igrp));
        ph = plot(xdata, ydata, '.b', 'MarkerSize', 10);
        set(gca, 'yscale', 'log', 'xscale', 'log')
        eqLim = [min([ylim, xlim]), max([ylim, xlim])];
        plot(eqLim, eqLim, '--k', 'LineWidth', 1)
        xlim(eqLim)
        ylim(eqLim)
        ylabel('MFR during ripples')
        xlabel('MFR during non-ripple bouts')
        title(axh, 'spike gain')
    end
    ph.Color = 'r';
    legend({'WT', 'MCU-KO'})
end

