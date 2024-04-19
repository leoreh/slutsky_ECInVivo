
function mcu_psd(mname)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% investigate states defined by emg during baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TO DO LIST
% recalc firing rates accorindg to new epochs
% split psd struct to light and dark phase
% use new psd epochs for epoch stats analysis

% DONE:
% lh96 (t1)
% lh122 (t3)
% lh142 (t2)

% files
basepaths = [mcu_sessions(mname)];

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

    % print progress
    fprintf('working on session %d of %d, %s\n',...
        ifile, nfiles, basename)

    % load emg signal
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
        'flgCalc', false, 'flgForce', false, 'graphics', false);

    % calc psd according to as state separation
    psd = psd_states('basepath', basepath, 'sstates', sstates,...
        'sig', sig, 'fs', fs, 'saveVar', saveVar,...
        'graphics', true, 'forceA', true, 'ftarget', ftarget,...
        'emgThr', [], 'flgEmg', false);

    % calc psd accordign to emg state separation
    psdEmg = psd_states('basepath', basepath, 'sstates', [1, 2],...
        'sig', sig, 'fs', fs, 'saveVar', saveVar,...
        'graphics', false, 'forceA', true, 'ftarget', ftarget,...
        'emgThr', [], 'flgEmg', true);

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
        plot_hypnogram('stateEpochs', psd.epochs.stateEpochs, 'sstates', [1 : 3],...
            'axh', axh, 'yshift', 3, 'clr', clr)
        plot_hypnogram('stateEpochs', psdEmg.epochs.stateEpochs, 'sstates', [1, 2],...
            'axh', axh, 'clr', clrEmg, 'yshift', 1.5)
        plot([1 : length(emg)], emg, 'k', 'LineWidth', 0.5)
        axis tight
        yLimit = ylim;
        ylim([yLimit(1) - 1, yLimit(2) + 1])
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
        pcMat = psd.epochs.pc;
        clrIdx = psd.epochs.clrIdx;
        sh = scatter3(pcMat(:, 1), pcMat(:, 2), pcMat(:, 3),...
            30, clrIdx, 'filled');
        view(-25, 25)
        grid minor
        sh.MarkerFaceAlpha = 0.6;
        xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
        title(axh(5), 'AS Feature Space')

        % AS silhouette plot
        axh(6) = nexttile(th, 10, [1, 1]); cla; hold on
        nonEmptyStates = cellfun(@(x) ~isempty(x), psd.epochs.stateEpochs, 'uni', true);
        stateIdx = psd.epochs.stateIdx;
        silhouette(pcMat, stateIdx, 'Euclidean');
        xlim([-1 1])
        yticklabels(snames(nonEmptyStates))
        ylabel('')
        sh = get(gca, 'Children');
        sh.FaceColor = 'flat';
        sh.CData(~isnan(sh.YData), :) = clrIdx;
        title(axh(6), 'AS Separation Quality')
        yTickVals = yticks;
        cnt = 1;
        for istate = 1 : 3
            meanSil = mean(psd.epochs.sil(stateIdx == istate));
            if ~isnan(meanSil)
                text(-0.75, yTickVals(cnt), num2str(meanSil, '%.3f'));
                cnt = cnt + 1;
            end
        end

        % EMG state clusters in feature space
        axh(7) = nexttile(th, 11, [1, 1]); cla; hold on
        pcMat = psdEmg.epochs.pc;
        clrIdx = psdEmg.epochs.clrIdx;
        sh = scatter3(pcMat(:, 1), pcMat(:, 2), pcMat(:, 3),...
            30, clrIdx, 'filled');
        view(-25, 25)
        grid minor
        sh.MarkerFaceAlpha = 0.6;
        xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
        title(axh(7), 'EMG Feature Space')

        % silhouette plot
        axh(8) = nexttile(th, 12, [1, 1]); cla; hold on
        stateIdx = psdEmg.epochs.stateIdx;
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
            meanSil = mean(psdEmg.epochs.sil(stateIdx == istate));
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
[bands, powdb] = sessions_psd('mname', mname, 'flgNormBand', true, 'flgDb', false,...
    'flgNormTime', false, 'flgEmg', true, 'idxBsl', [1 : 2], 'saveFig', false);

% load all figures
setMatlabGraphics(true)
for ifile = 1 : nfiles
    basepath = basepaths{ifile};
    [~, basename] = fileparts(basepath);
    figpath = fullfile(basepath, 'graphics', 'sleepState');
    figname = fullfile(figpath, [basename, '_ASvsEMG']);
    openfig(figname)
end

end
