

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = as_loadConfig;
sstates = [1, 4];

mname = 'lh132';
mcu_psd(mname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load files of single mouse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mname = 'lh132';

varsFile = ["fr"; "sr"; "sleep_states";...
    "datInfo"; "session"; "units"; "psd"];
varsName = ["fr"; "sr"; "ss"; "datInfo"; "session";...
    "units"; "psd"];
xlsname = 'D:\OneDrive - Tel-Aviv University\PhD\Slutsky\Data summaries\sessionList.xlsx';
[v, basepaths] = getSessionVars('mname', mname, 'varsFile', varsFile,...
    'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
    'xlsname', xlsname);
nfiles = length(basepaths);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% average bands and psd from multiple mice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mname = {'lh96', 'lh107', 'lh122', 'lh142'};
mname = {'lh132', 'lh133', 'lh134', 'lh136', 'lh140'};
flgNormTime = false;
flgNormBand = true;

clear bands powdb
% bands(imouse, iband, istate, isession)
% powdb(imouse, istate, ifreq, isession)
for imouse = 1 : length(mname)

    [bands(imouse, :, :, :), powdb(imouse, :, :, :)] =...
        sessions_psd(mname{imouse}, 'flgNormBand', flgNormBand,...
        'flgAnalyze', false, 'flgNormTime', flgNormTime,...
        'flgEmg', true, 'idxBsl', [1], 'graphics', false, 'flgDb', false);

end

% get mouse data (transpose s.t. sessions are rows). if flgNormBand, ignore
% first column (broadband)
istate = 2;
imouse = 5;
prismData = squeeze(bands(imouse, :, istate, :))';
prismData = squeeze(powdb(imouse, istate, :, :));

% reshape powdb to 2d for grouped graph in prism
x = squeeze(powdb(:, istate, :, :));
prismData = [];
cnt = 1;
for ifile = 1 : size(x, 3)
    for imouse = 1 : length(mname)
        prismData(:, cnt) = squeeze(x(imouse, :, ifile));
        cnt = cnt + 1;
    end
end

% reshape bands to 2d for grouped graph in prism
x = squeeze(bands(:, :, istate, :));
prismData = [];
cnt = 1;
for iband = 1 : size(x, 2)
    for imouse = 1 : length(mname)
        prismData(:, cnt) = squeeze(x(imouse, iband, :));
        cnt = cnt + 1;
    end
end


% -------------------------------------------------------------------------
% for each freq point, calc diff between bsl and bac effect, and plot as
% function of freq point. note normBand must be false
istate = 2;
isessions = [1, 3];
x = v(1).psd.info.faxis;
y1 = squeeze(powdb(2, istate, :, isessions(1)));
y2 = squeeze(powdb(2, istate, :, isessions(2)));
% y1 = 10 .^ (y1 / 10);
% y2 = 10 .^ (y2 / 10);
y = (y2 - y1) ./ (y2 + y1);
fh = figure;
plot(x, y)
hold on
xLimit = xlim;
plot(xLimit, [0, 0], '--k')
legend
% set(gca, 'xscale', 'log')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize data - psd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psd = catfields([v(:).psd], 'catdef', 'cell', 'force', false);

% ORGANIZE: statePsd is a matrix of frequency (rows) x session (column)
% depicting the psd for the selected state
clear statePsd
istate = 4;
for ifile = 1 : nfiles

    statePsd(:, ifile) = squeeze(psd.psd{ifile}(1, istate, :));

end
statePsd = statePsd ./ sum(statePsd);
freq = psd.info.faxis{1};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hypnogram during baseline / washout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setMatlabGraphics(true)

% files
basepaths = [mcu_basepaths('wt_bsl')];
idx_wt = length(basepaths);
idx_mcu = idx_wt + 1;
basepaths = [basepaths; mcu_basepaths('mcu_bsl')];
varsFile = ["fr"; "sleep_states"; "datInfo"; "session"; "units"];
varsName = ["fr"; "ss"; "datInfo"; "session"; "units"];
xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
[v, basepaths] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
    'xlsname', xlsname);
nfiles = length(basepaths);

% params
nbins = 2;
sstates = [1, 4, 5];
snames = v(1).ss.info.names(sstates);

% get bout stats
clear prctDur boutLen totDur
for ifile = 1 : nfiles

    % timebins for calculating stats
    if nbins == 4
        timebins = v(ifile).session.general.timebins;
    else
        nlabels = length(v(ifile).ss.labels);
        timebins = n2chunks('n', nlabels, 'nchunks', nbins);
    end

    labels = v(ifile).ss.labels;
    stSwitch = [6, 5; 3, 4; 2, 1];

    for i_sw = 1 : size(stSwitch, 1)
        labels(labels == stSwitch(i_sw, 1)) = stSwitch(i_sw, 2);
    end
    labels = v(ifile).ss.labels;
    minDur = [20, 5, 5, 20, 10, 5];
    minDur = 0;

    cd(basepaths{ifile})
    [totDur(ifile, :, :), prctDur(ifile, :, :), boutLen(ifile, :)] = as_plotZT('nbins', nbins,...
        'sstates', sstates, 'ss', v(ifile).ss, 'labels', labels,...
        'timebins', timebins, 'graphics', false, 'minDur', minDur);
end

% graphics
for ibin = 1 : nbins

    fh = figure;
    tlayout = [3, 3];
    th = tiledlayout(tlayout(1), tlayout(2));
    th.TileSpacing = 'tight';
    th.Padding = 'none';

    nbouts = cellfun(@(x) sum(~isnan(x(:, ibin))), boutLen, 'UniformOutput', true);

    for istate = 1 : length(sstates)

        % state duration (%)
        axh = nexttile(th, istate, [1, 1]);
        hold on

        dataMat = prctDur;
        dataMat = cell2nanmat({dataMat(1 : idx_wt, ibin, istate),...
            dataMat(idx_mcu : end, ibin, istate)}, 2);
        plot_boxMean('dataMat', dataMat, 'clr', 'kr', 'allPnts', true)
        ylabel('State Duration (%)');
        legend({'WT', 'MCU-KO'}, 'Location', 'best');
        title(axh, snames(istate))

        % bout length per mouse
        axh = nexttile(th, istate + length(sstates), [1, 1]);
        hold on

        dataMat =  cellfun(@(x) mean(x(:, ibin), 'omitnan'),...
            boutLen(:, istate), 'UniformOutput', true);
        dataMat = cell2nanmat({dataMat(1 : idx_wt), dataMat(idx_mcu : end)}, 2);
        plot_boxMean('dataMat', dataMat, 'clr', 'kr', 'allPnts', true)
        ylabel('Bout Length (s)');
        title(axh, snames(istate))

        % nbouts per mouse
        axh = nexttile(th, istate + 2 * length(sstates), [1, 1]);
        hold on
        dataMat = nbouts(:, istate);
        dataMat = cell2nanmat({dataMat(1 : idx_wt), dataMat(idx_mcu : end)}, 2);
        plot_boxMean('dataMat', dataMat, 'clr', 'kr', 'allPnts', true)
        ylabel('No. Bouts');
        title(axh, snames(istate))

    end

    if ibin == 1
        txt = 'LightPhase';
    else
        txt = 'DarkPhase';
    end
    title(th, txt)

end

%%%

% manually check states
for ifile = 1 : nfiles
    cd(basepaths{ifile})
    [~, basename] = fileparts(basepaths{ifile})
    sSig = load([basename, '.sleep_sig.mat']);
    AccuSleep_viewer(sSig, v(ifile).ss.labels, [])
end


% plot hynogram
fh = figure;
for ifile = 1 : nfiles

    sb = subplot(nfiles, 1, ifile);
    plot_hypnogram('boutTimes', v(ifile).ss.boutTimes, 'axh', sb)

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare state separation across time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%
% multiple ways to show renormalization of state-dependent psd:
% (1) delta-theta ratio across time
% (2) difference between delta (and other bands) between low and high emg
% (3) accusleep recovers its ability to classify states
% (4) plot emg vs. delta theta for bouts accusleep classified as AW and
% NREM and measure distance

% -------------------------------------------------------------------------
% state duration for each day as automatically classified

mname = {'lh96', 'lh107', 'lh122', 'lh142'};
mname = {'lh132', 'lh133', 'lh134', 'lh136', 'lh140'};

nmice = length(mname);
nbins = 1;
close all

for imouse = 1 : nmice

    % load data
    varsFile = ["sleep_states"; "datInfo"; "session"];
    varsName = ["ss"; "datInfo"; "session"];
    xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
    [v, basepaths] = getSessionVars('mname', mname{imouse}, 'varsFile', varsFile,...
        'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
        'xlsname', xlsname);
    nfiles = length(basepaths);

    stateDur = nan(nfiles * nbins, 3);
    cnt = 1;
    for ifile = 1 : nfiles
        % reclassify to 3 states
        labels = v(ifile).ss.labels;
        labels(labels == 2) = 1; % QW => AW
        labels(labels == 3) = 4; % LSLEEP => NREM
        labels(labels == 6) = 5; % N/REM => REM

        % create timebins or take from sessionInfo
        if nbins == 4
            timebins = v(ifile).session.general.timebins;
        else
            timebins = n2chunks(length(labels), 'nchunks', nbins);
        end

        % count total state dur according to timebins
        tmp = [];
        for ibin = 1 : nbins
            binIdx = timebins(ibin, 1) : timebins(ibin, 2);
            tmp(1) = sum(labels(binIdx) == 1);
            tmp(2) = sum(labels(binIdx) == 4);
            tmp(3) = sum(labels(binIdx) == 5);
            stateDur(cnt, :) = tmp ./ numel(binIdx) * 100;
            cnt = cnt + 1;
        end
    end

    % plot
    fh = figure;
    b = bar(stateDur, 'stacked', 'FaceColor', 'flat');
    cdata = v(ifile).ss.info.colors([1, 4, 5]);
    for iclr = 1 : size(stateDur, 2)
        b(iclr).CData = cdata{iclr};
    end
    title(mname{imouse})
    xticklabels({'BSL', 'BAC On', 'BAC1', 'BAC2', 'BAC3', 'BAC Off', 'WASH'})
    legend({'AW', 'NREM', 'REM'})
    % grab specific bins from all mice
    % this assumes analyzed 2 bins per day
    % bins nrem from light cycle: bsl = 1; acute = 5; chronic = 9, wash = 13.
    % bins aw from dark cycle: bsl = 2; acute = 4; chronic = 10, wash = 14.

end


% -------------------------------------------------------------------------
% state dur for each day as automatically classified

% 1. calc spec in 1-second bins
% 2. calc emg rms in 1-second bins, categorize high vs. low
% 3. calc PCs from spec and project
% 4. measure distance between high and low clusters of spec w/ mahalanobis

mname = {'lh96', 'lh107', 'lh122', 'lh142'};
mname = {'lh132', 'lh133', 'lh134', 'lh136', 'lh140'};

nmice = length(mname);
nbins = 1;
close all

for imouse = 1 : nmice

    % load data
    varsFile = ["sleep_states"; "session"; 'spec'];
    varsName = ["ss"; "session"; 'spec'];
    xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
    [v, basepaths] = getSessionVars('mname', mname{imouse}, 'varsFile', varsFile,...
        'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
        'xlsname', xlsname);
    nfiles = length(basepaths);

    fh = figure;
    figAlt = 2;         % 1 - delta / theta; 2 - pca
    for ifile = 1 : nfiles

        % grab data, ALT1 - from sSig (1s bins)
        basepath = v(ifile).session.general.basePath;
        cd(basepath)
        basename = v(ifile).session.general.name;
        load([basename, '.sleep_sig.mat'], 'info');
        load([basename, '.sleep_sig.mat'], 'spec');
        load([basename, '.sleep_sig.mat'], 'spec_freq');
        load([basename, '.sleep_sig.mat'], 'emg_rms');

        % grab data, ALT2 - calculate / load (5s bins)
        % ch = 4;
        % spec = squeeze(v(ifile).spec.s(:, :, ch));
        % load([basename, '.sleep_sig.mat'], 'emg');
        % emg_rms = processEMG(emg, 1250, 5);

        % normalize spec w/ as calibration
        specNorm = log(spec);
        calData = v(ifile).ss.info.calibrationData;
        for ifreq = 1 : size(specNorm, 2)
            specNorm(:, ifreq) = (specNorm(:, ifreq) - calData(ifreq, 1)) ./ calData(ifreq, 2);
            specNorm(:, ifreq) = (specNorm(:, ifreq) + 4.5) ./ 9; % clip z scores
        end
        % clip
        specNorm(specNorm < 0) = 0;
        specNorm(specNorm > 1) = 1;

        % measure delta to theta ratio
        [~, f1idx] = min(abs(spec_freq - 1));
        [~, f4idx] = min(abs(spec_freq - 4));
        [~, f6idx] = min(abs(spec_freq - 6));
        [~, f12idx] = min(abs(spec_freq - 12));
        sDelta = sum(specNorm(:, f1idx : f4idx), 2);
        sTheta = sum(specNorm(:, f6idx : f12idx), 2);
        sRatio = sDelta ./ sTheta;

        % 2. calc emg rms and categorize high vs. low
        emgThr = 70;
        highIdx = emg_rms > prctile(emg_rms, emgThr);
        lowIdx = emg_rms < prctile(emg_rms, 100 - emgThr);

        % 3. calc PCs from spec
        [coeff, pcEmg, ~] = pca(spec, 'NumComponents', 3);

        % 4. measure distance
        sidx = highIdx | lowIdx;
        [lRat, iDist, mDist] = cluDist(pcEmg(sidx, :), highIdx(sidx));
        d = median(mDist)

        % plot
        axh = nexttile;
        sz = 20;
        falpha = 0.2;

        if figAlt == 1
            % ALT 1 - emg and delta / theta
            scatter(sRatio(highIdx), emg_rms(highIdx),...
                sz, 'r', 'filled', 'MarkerFaceAlpha', falpha)
            hold on
            scatter(sRatio(lowIdx), emg_rms(lowIdx),...
                sz, 'b', 'filled', 'MarkerFaceAlpha', falpha)
            legend('High-EMG', 'Low-EMG');
            xlabel('delta / theta'); ylabel('EMG_RMS')
            title(d)

            % ALT 2 - pca
        else
            scatter3(pcEmg(highIdx, 1), pcEmg(highIdx, 2), pcEmg(highIdx, 3),...
                sz, 'r', 'filled', 'MarkerFaceAlpha', falpha, 'MarkerEdgeColor', 'none');
            hold on;
            scatter3(pcEmg(lowIdx, 1), pcEmg(lowIdx, 2), pcEmg(lowIdx, 3),...
                sz, 'b', 'filled', 'MarkerFaceAlpha', falpha, 'MarkerEdgeColor', 'none');
            legend('High-EMG', 'Low-EMG');
            xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
            title([basename, '. mDist = ', d])
        end


    end



end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% investigate states defined by emg during baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TO DO LIST
% recalc firing rates accorindg to new bouts

% DONE:
% lh96 (t1)
% lh122 (t3)

% files
basepaths = [mcu_basepaths('wt_bsl')];
mname = 'lh122';
basepaths = [mcu_basepaths(mname)];

% load data
varsFile = ["sleep_states"; "datInfo"; "session"];
varsName = ["ss"; "datInfo"; "session"];
xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
[v, basepaths] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
    'xlsname', xlsname);
nfiles = length(basepaths);

% flags
graphics = true;
saveVar = true;
saveFig = true;

% params
sstates = [1, 4, 5];
clr = v(1).ss.info.colors(sstates);
clrEmg = {[150, 70, 55] / 255, [200, 170, 100] / 255};
snames = v(1).ss.info.names(sstates);
ftarget = [0.5 : 0.5 : 100];

for ifile = 1 : nfiles

    % file
    basepath = basepaths{ifile};
    [~, basename] = fileparts(basepath);
    cd(basepath)
    timeWin = [v(ifile).session.general.timepnt;,...
        floor(v(ifile).session.general.duration)];
    timeWin = [seconds(hours(10)), Inf];
    timeWin = [];

    % print progress
    fprintf('working on session %d of %d, %s\n',...
        ifile, nfiles, basename)

    % load emg signals
    sleepfile = fullfile(basepath, [basename, '.sleep_sig.mat']);
    sig = load(sleepfile, 'eeg');
    sig = sig.eeg;
    emg = load(sleepfile, 'emg_rms');
    emg = emg.emg_rms;
    load(sleepfile, 'fs');

    % load spectrogram
    s = load(sleepfile, 'spec');
    load(sleepfile, 'spec_freq')
    load(sleepfile, 'spec_tstamps')
    spec.s = s.spec; spec.freq = spec_freq; spec.tstamps = spec_tstamps;
    otl = get_specOutliers('basepath', basepath, 'saveVar', true,...
        'flgCalc', false, 'flgForce', false, 'graphics', true);

    % calc psd according to as state separation
    psd = psd_states('basepath', basepath, 'sstates', sstates,...
        'sig', sig, 'fs', fs, 'saveVar', saveVar,...
        'graphics', graphics, 'forceA', true, 'ftarget', ftarget,...
        'emgThr', [], 'flgEmg', false, 'timeWin', timeWin);

    % calc psd accordign to emg state separation
    psdEmg = psd_states('basepath', basepath, 'sstates', [1, 2],...
        'ch', [], 'fs', [], 'saveVar', true,...
        'graphics', graphics, 'forceA', saveVar, 'ftarget', ftarget,...
        'emgThr', [], 'flgEmg', true, 'timeWin', timeWin);

    % graphics ------------------------------------------------------------
    if graphics

        fh = figure; clear axh;
        setMatlabGraphics(true)
        set(fh, 'WindowState','maximized');
        tlayout = [3, 4];
        th = tiledlayout(tlayout(1), tlayout(2));
        th.TileSpacing = 'tight';
        th.Padding = 'none';
        title(th, basename, 'interpreter', 'none')

        % superposition of state hypnogram and emg-states
        axh(1) = nexttile(th, 1, [1, 3]); cla; hold on
        plot_hypnogram('boutTimes', psd.bouts.boutTimes, 'sstates', [1 : 3],...
            'axh', axh, 'yshift', 3, 'clr', clr)
        plot_hypnogram('boutTimes', psdEmg.bouts.boutTimes, 'sstates', [1, 2],...
            'axh', axh, 'clr', clrEmg, 'yshift', 1)
        plot([1 : length(emg)], emg, 'k', 'LineWidth', 0.5)
        axis tight
        yLimit = ylim;
        ylim([yLimit(1) - 1, yLimit(2)])
        xval = [3600 : 3600 : length(spec.s)];
        xticks(xval);
        xticklabels(string(xval / 3600))
        xlabel('Time (h)')
        title(axh(1), 'Hypnogram')

        % histogram of emg values
        axh(2) = nexttile(th, 4, [1, 1]);
        hold on
        emgThr = psdEmg.info.emgThr;
        [counts, edges] = histcounts(emg);
        bCents = edges(1 : end - 1) + diff(edges) / 2;
        bh = bar(bCents, counts, 'BarWidth', 1);
        bh.FaceColor = 'flat';
        bh.CData(bh.XData > emgThr, :) = repmat(clrEmg{1}, sum(bh.XData > emgThr), 1);
        bh.CData(bh.XData < emgThr, :) = repmat(clrEmg{2}, sum(bh.XData < emgThr), 1);
        plot([emgThr, emgThr], ylim, '--k')
        xlabel('EMG')
        ylabel('Counts')
        title(axh(2), 'EMG Distribution')

        % spectrogram
        axh(3) = nexttile(th, 5, [1, 3]); cla
        hold on
        plot_spec(spec, 'ch', 1, 'logfreq', false, 'saveFig', false,...
            'axh', axh(3), 'xtime', 1)
        axis tight
        yLimit = ylim;
        scatter(otl.idx, ones(length(otl.idx), 1) * yLimit(2), 20, 'filled', 'r')
        xval = [3600 : 3600 : length(spec.s)];
        xticks(xval);
        xticklabels(string(xval / 3600))
        linkaxes([axh(1), axh(3)], 'x')
        hold on

        % psd comparison between AS and EMG states
        axh(4) = nexttile(th, 8, [1, 1]);
        hold on
        ph = plot(psd.info.faxis, squeeze(psd.psd([1, 2], :)), 'LineWidth', 2);
        for istate = 1 : 2
            ph(istate).Color = clr{istate};
        end
        ph = plot(psd.info.faxis, squeeze(psdEmg.psd), 'LineWidth', 2);
        for istate = 1 : 2
            ph(istate).Color = clrEmg{istate};
        end
        set(gca, 'YScale', 'log', 'XScale', 'log')
        xlabel('Frequency [Hz]')
        ylabel('PSD [mV^2/Hz]')
        title(axh(4), 'Average PSD')
        legend({'WAKE', 'NREM', 'High-EMG', 'Low-EMG'})

        % AS state clusters in feature space
        axh(5) = nexttile(th, 9, [1, 1]); cla; hold on
        pcMat = psd.bouts.pc;
        clrIdx = psd.bouts.clrIdx;
        sh = scatter3(pcMat(:, 1), pcMat(:, 2), pcMat(:, 3),...
            30, clrIdx, 'filled');
        view(-25, 25)
        grid minor
        sh.MarkerFaceAlpha = 0.6;
        xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
        title(axh(5), 'AS Feature Space')

        % AS silhouette plot
        axh(6) = nexttile(th, 10, [1, 1]); cla; hold on
        stateIdx = psd.bouts.stateIdx;
        silhouette(pcMat, stateIdx, 'Euclidean');
        xlim([-1 1])
        yticklabels(snames)
        ylabel('')
        sh = get(gca, 'Children');
        sh.FaceColor = 'flat';
        sh.CData(~isnan(sh.YData), :) = clrIdx;
        title(axh(6), 'AS Separation Quality')
        yTickVals = yticks;
        for istate = 1 : 3
            meanSil = mean(psd.bouts.sil(stateIdx == istate));
            if ~isnan(meanSil)
                text(-0.75, yTickVals(istate), num2str(meanSil, '%.3f'));
            end
        end

        % EMG state clusters in feature space
        axh(7) = nexttile(th, 11, [1, 1]); cla; hold on
        pcMat = psdEmg.bouts.pc;
        clrIdx = psdEmg.bouts.clrIdx;
        sh = scatter3(pcMat(:, 1), pcMat(:, 2), pcMat(:, 3),...
            30, clrIdx, 'filled');
        view(-25, 25)
        grid minor
        sh.MarkerFaceAlpha = 0.6;
        xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
        title(axh(7), 'EMG Feature Space')

        % silhouette plot
        axh(8) = nexttile(th, 12, [1, 1]); cla; hold on
        stateIdx = psdEmg.bouts.stateIdx;
        silhouette(pcMat, stateIdx, 'Euclidean');
        xlim([-1 1])
        yticklabels(snames)
        ylabel('')
        sh = get(gca, 'Children');
        sh.FaceColor = 'flat';
        sh.CData(~isnan(sh.YData), :) = clrIdx;
        title(axh(8), 'EMG Separation Quality')
        yTickVals = yticks;
        for istate = 1 : 2
            meanSil = mean(psdEmg.bouts.sil(stateIdx == istate));
            text(-0.75, yTickVals(istate), num2str(meanSil, '%.3f'));
        end

        if saveFig
            figpath = fullfile(basepath, 'graphics', 'sleepState');
            mkdir(figpath)
            figname = fullfile(figpath, [basename, '_ASvsEMG']);
            savefig(fh, figname, 'compact')
        end
    end

end

% get psd across sessions for one mouse
mname = 'lh122';
[bands, powdb] = sessions_psd(mname, 'flgNormBand', false,...
    'flgNormTime', true, 'flgEmg', true, 'idxBsl', [1], 'saveFig', false);

% load all figures
setMatlabGraphics(true)
for ifile = 1 : nfiles
    basepath = basepaths{ifile};
    [~, basename] = fileparts(basepath);
    figpath = fullfile(basepath, 'graphics', 'sleepState');
    figname = fullfile(figpath, [basename, '_ASvsEMG']);
    openfig(figname)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot psd across experiment per mouse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% files
mname = 'lh132';
basepaths = [mcu_basepaths(mname)];

% load data
varsFile = ["fr"; "sleep_states"; "datInfo"; "session"; "psd"; "psdEmg"];
varsName = ["fr"; "ss"; "datInfo"; "session"; "psd"; "psdEmg"];
xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
[v, basepaths] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
    'xlsname', xlsname);
nfiles = length(basepaths);

% TO DO LIST psd_states
% add tstamps
% add bands to mean psd


% flags
graphics = true;
saveVar = true;
saveFig = true;

for ifile = 1 : nfiles

    basepath = basepaths{ifile};
    [~, basename] = fileparts(basepath);
    cd(basepath)

    for istate = 1 : length(sstates)
        tstamps{istate} = diff(v(ifile).psdEmg.info.boutTimes{istate}')' +...
            v(ifile).psdEmg.info.boutTimes{istate}(:, 1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot psd during baseline per mouse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% files
mname = 'wt_wsh';
vars = ["session"; "psd"];
mNames = [unique(get_mname(mcu_basepaths('eeg')))];
basepaths = mcu_basepaths(mname);
v = basepaths2vars('basepaths', basepaths, 'vars', vars);
nfiles = length(basepaths);

psd = catfields([v(:).psd], 'addim', true);
istate = 1;
prismData = squeeze(psd.psd(istate, :, :));
prismData = mean(prismData ./ sum(prismData), 2);

psd.info{1}.faxis



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% broadband power during baseline and washout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grps
clear mname
mname{1} = 'wt_bsl';
mname{2} = 'wt_wsh';
mname{3} = 'mcu_bsl';
mname{4} = 'mcu_wsh';
ngrp = length(mname)

bband = nan(ngrp, 6, 3);
for igrp = 1 : ngrp

    % load data
    vars = ["psd"];
    basepaths = mcu_basepaths(mname{igrp});
    v = basepaths2vars('basepaths', basepaths, 'vars', vars);
    nfiles = length(basepaths);

    % nunits
    for ifile = 1 : nfiles
        bband(igrp, ifile, :) = v(ifile).psd.bands.mean(:, 1);
    end
end

% for grp plot
istate = 2;
sz = size(bband);
prismData = squeeze(bband([1 : 2], :, istate));
prismData(:, sz(2) + 1 : sz(2) * 2) = squeeze(bband([3, 4], :, istate));

% for individual lines plot
istate = 2;
prismData = squeeze(bband(:, :, istate));

