function psd = psd_states(varargin)

% wrapper for calc_psd dedicated to sleep states. uses the eeg signal from
% sSig by default, but can load any specified ch from a binary file.
% if sleep_states.mat is not found, will separate the recording to "AW" and
% "NREM" according to high- and low-emg activity.
%
% INPUT:
%   basepath        string. path to recording folder {pwd}
%   sig             vector of lfp signal. if empty will grab from sSig or
%                   sigfile (see below)
%   timeWin         2 element numeric vector depicting the time (s) over
%                   which to calculate btimes, i.e. the time window for
%                   calculating the mean psd. used as indices to labels.
%   ch              numeric. channels to load from the sigfile. if a vector
%                   is specified, the signal will be averaged across
%                   channels
%   nchans          numeric. no. channels in sigfile
%   fs              numeric. sampling frequency of lfp data. if empty will
%                   be extracted from session struct. also determines the
%                   new sampling frequency of the eeg signal.
%   sigfile         char. name of file to load signal from. if empty but
%                   channel is specified, will load from [basename.lfp]
%   btimes          cell of n x 2 mats. if empty will load from sleep states 
%   sstates         numeric. index of selected states to calculate psd
%   ftarget         numeric. requested frequencies for calculating the psd
%   emgThr          numeric. percent by which to seprate high- and low-emg.
%                   e.g., prct = 50 is the median
%   flgEmg          logical. calc psd in high- and low-emg even if states
%                   file exists
%   saveVar         logical. save ss var {true}
%   forceA          logical. reanalyze recordings even if ss struct
%                   exists (false)
%   graphics        logical. plot confusion chart and state separation {true}
%
% DEPENDENCIES:
%   calc_psd
%
% TO DO LIST:
%
% 07 sep 22 LH  updates:
% 20 mar 24 LH      cleanup, emg params, band analysis
% 26 mar 24 LH      removed win functionality since psd is now calculated
%                   per bout. much simpler this way. bkup exists in temp
%                   folder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p, 'basepath', pwd);
addOptional(p, 'sig', [], @isnumeric);
addOptional(p, 'ch', [], @isnumeric);
addOptional(p, 'timeWin', [], @isnumeric);
addOptional(p, 'nchans', [], @isnumeric);
addOptional(p, 'fs', 1, @isnumeric);
addOptional(p, 'sigfile', [], @ischar);
addOptional(p, 'sstates', [1, 4, 5], @isnumeric);
addOptional(p, 'btimes', []);
addOptional(p, 'ftarget', 0.5 : 0.5 : 100, @isnumeric);
addOptional(p, 'emgThr', 50, @isnumeric);
addOptional(p, 'flgEmg', false, @islogical);
addOptional(p, 'saveVar', true, @islogical);
addOptional(p, 'forceA', false, @islogical);
addOptional(p, 'graphics', true, @islogical);

parse(p, varargin{:})
basepath        = p.Results.basepath;
sig             = p.Results.sig;
timeWin         = p.Results.timeWin;
ch              = p.Results.ch;
nchans          = p.Results.nchans;
fs              = p.Results.fs;
sigfile         = p.Results.sigfile;
sstates         = p.Results.sstates;
btimes          = p.Results.btimes;
ftarget         = p.Results.ftarget;
emgThr          = p.Results.emgThr;
flgEmg          = p.Results.flgEmg;
saveVar         = p.Results.saveVar;
forceA          = p.Results.forceA;
graphics        = p.Results.graphics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% analysis params
nfet = 3;

% graphics flags
flgSaveFig = true;

% file
cd(basepath)
[~, basename] = fileparts(basepath);
sleepfile = fullfile(basepath, [basename, '.sleep_sig.mat']);

% load session vars
varsFile = ["sleep_states"; "datInfo"; "session"];
varsName = ["ss"; "datInfo"; "session"];
v = getSessionVars('basepaths', {basepath}, 'varsFile', varsFile,...
    'varsName', varsName);
recLen = length(v.ss.labels);
spkgrp = v.session.extracellular.spikeGroups.channels;

% state params
if flgEmg
    flgEmg = true;
    namePrefix = 'psdEmg';
    sstates = [1 : 2];
    cfg = as_loadConfig('flgEmg', true);
    statesfile = fullfile(basepath, [basename, '.sleep_statesEmg.mat']);
else
    flgEmg = false;
    namePrefix = 'psd';
    cfg = as_loadConfig();
    if isempty(sstates)
        sstates = 1 : cfg.nstates;
    end
    statesfile = fullfile(basepath, [basename, '.sleep_states.mat']);
end
clr = cfg.colors(sstates);
snames = cfg.names(sstates);
nstates = length(sstates);

% check if already analyzed
psdfile = fullfile(basepath, [basename, '.', namePrefix, '.mat']);
if exist(psdfile, 'file') && ~forceA
    load(psdfile)
    return
end

% assert time window
if isempty(timeWin)
    timeWin = [1, recLen];
else
    if timeWin(1) < 1
        timeWin(1) = 1;
    end
    if timeWin(2) > recLen
        timeWin(2) = recLen;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load signal from sSig or from binary if ch specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(sig)
    if ~isempty(ch)         % from binary

        % if no file specified, load from .lfp binary
        if isempty(sigfile)
            sigfile = fullfile(basepath, [basename, '.lfp']);
            fs = v.session.extracellular.srLfp;
            nchans = v.session.extracellular.nChannels;
        end

        sig = double(bz_LoadBinary(sigfile,...
            'duration', Inf,...
            'frequency', fs, 'nchannels', nchans, 'start', 0,...
            'channels', ch, 'downsample', 1));

        % average tetrode
        if length(ch) > 1
            sig = mean(sig, 2);
        end

    else                    % from sSig

        sig = load(sleepfile, 'eeg');
        sig = sig.eeg;
        load(sleepfile, 'fs');
        load(sleepfile, 'info');
        ch = info.eegCh;
    end
else                        % assumes input signal is from sSig
    load(sleepfile, 'info');
    ch = info.eegCh;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate state bouts if not provided
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get btimes
if isempty(btimes)
    ss = load(statesfile);
    varName = fieldnames(ss);   
    ss = ss.(varName{1});
    bouts = ss.bouts;
    btimes = bouts.times;

    % add timeWin(1) to each bout
    btimes = cellfun(@(x) round(x + timeWin(1) - 1),...
        btimes, 'uni', false);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc psd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, faxis, psd_bouts] = calc_psd('sig',...
    sig, 'bins', btimes,...
    'fs', fs, 'ftarget', ftarget, 'graphics', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove outlier bouts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

otl = get_otl(data, 'graphics', true);


% outliers are detected in several iterations, first by pc1 and then by
% mahalnobis distance. the pca is recalculated between iterations. this
% was determined by empiric examination of the data (march 2024).

% initialize vars
clear otl
psd.bouts.original = psd_bouts;
tmp_bouts = psd_bouts;
otlIdx = cell(nstates, 1);
origIdx = cell(nstates, 1);
for istate = 1 : nstates
    origIdx{istate} = 1 : size(btimes{istate}, 1);
end

% project bout PSDs on feature space
[~, pc, ~, ~, expl] = pca(vertcat(tmp_bouts{:}), 'NumComponents', nfet);

% create state and color indices for each bout when bouts are
% concatenated (e.g., for PCs)
stateIdx = []; clrIdx = [];
for istate = 1 : nstates
    stateIdx = [stateIdx; istate * ones(bouts.nbouts(istate), 1)];
    clrIdx = [clrIdx; repmat(clr{istate}, bouts.nbouts(istate), 1)];
end
startIdx = [0, cumsum(bouts.nbouts(1 : end - 1))];

% detect outliers
nrep = 1;
for irep = 1 : nrep

    % get pca outliers per state
    badIdx = cell(nstates, 1);
    for istate = 1 : nstates
        tmpIdx = stateIdx == istate;
        if sum(tmpIdx) > 20
            if irep < 3
                stdThr = 2;
                dists = pc(tmpIdx, 1);
            else
                stdThr = 2 * (irep - 1);
                dists = mahal(pc(tmpIdx, :), pc(tmpIdx, :));
            end

            badIdx{istate} = dists > mean(dists) + std(dists) * stdThr;
        else
            badIdx{istate} = false(sum(tmpIdx), 1);
        end

        % track indices to original bouts
        otlIdx{istate} = [otlIdx{istate}; origIdx{istate}(badIdx{istate})'];

    end

    % organize outliers struct as an array (of iterations)
    otl(irep).pc = pc;
    otl(irep).expl = expl;
    otl(irep).stateIdx = stateIdx;
    otl(irep).clrIdx = clrIdx;
    otl(irep).badIdx = badIdx;
    otl(irep).btimes = btimes;
    otl(irep).bouts = bouts;
    otl(irep).psd_bouts = tmp_bouts;
    otl(irep).sil = silhouette(pc, stateIdx, 'Euclidean');

    % remove bad indices from bouts and recalculate bout stats
    for istate = 1 : nstates
        if any(badIdx{istate})
            btimes{istate}(badIdx{istate}, :) = [];
            tmp_bouts{istate}(badIdx{istate}, :) = [];
            origIdx{istate}(badIdx{istate}) = [];
        end
    end

    % calculate bouts
    bouts.boutLen = cellfun(@(x) (diff(x')'), btimes, 'UniformOutput', false);
    bouts.nbouts = cellfun(@length, bouts.boutLen);
    bouts.totDur = cellfun(@sum, bouts.boutLen);

    % recalculate pca
    [~, pc, ~, ~, expl] = pca(vertcat(tmp_bouts{:}), 'NumComponents', nfet);

    % create state and color indices
    stateIdx = []; clrIdx = [];
    for istate = 1 : nstates
        stateIdx = [stateIdx; istate * ones(bouts.nbouts(istate), 1)];
        clrIdx = [clrIdx; repmat(clr{istate}, bouts.nbouts(istate), 1)];
    end

end

% organize global indices of outliers (to pc mat)
clear otlIdx_glbl
for istate = 1 : nstates
    otlIdx_glbl{istate} = otlIdx{istate} + startIdx(istate);
end
otlIdx_glbl = sort(cat(1, otlIdx_glbl{:}));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize psd struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% organize bouts
psd.bouts.clean = tmp_bouts;
psd.bouts.btimes = btimes;
psd.bouts.bouts = bouts;
psd.bouts.otl = otl;
psd.bouts.pc = pc;
psd.bouts.expl = expl;
psd.bouts.sil = silhouette(pc, stateIdx, 'Euclidean');
psd.bouts.stateIdx = stateIdx;
psd.bouts.clrIdx = clrIdx;
psd.bouts.otlIdx = otlIdx;
psd.bouts.otlIdx_glbl = otlIdx_glbl;

% info
psd.info.sstates = sstates;
psd.info.snames = snames;
psd.info.clr = clr; 
psd.info.faxis = faxis;
psd.info.emgThr = emgThr;
psd.info.timeWin = timeWin;
psd.info.ch = ch;
psd.info.spkgrpIdx = find(cellfun(@(x) all(ismember(ch, x)), spkgrp));
psd.info.runtime = datetime("now");

% calculate vars using cleaned bouts
for istate = 1 : length(sstates)
    if ~isempty(psd.bouts.clean{istate})
        psd.psd(istate, :) = mean(psd.bouts.clean{istate}, 1);

        % calc power in specific frequency bands
        % per bout
        [psd.bands.bouts{istate}, psd.bands.info] = calc_bands('psdData',...
            psd.bouts.clean{istate}, 'freq', faxis, 'flgNormBand', false);
        % across bouts
        [psd.bands.mean(istate, :), ~] = calc_bands('psdData',...
            psd.psd(istate, :), 'freq', faxis, 'flgNormBand', false);

    end
end

if saveVar
    save(psdfile, 'psd')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% plot spec, emg_rms, and hypnogram with outlier bouts in black.
% plot only the original and final steps of cleaning (perhaps keep the
% option to plot all steps) 


% load spectrogram for plot
if graphics
    if ~exist('emg', 'var')
        emg = load(sleepfile, 'emg_rms');
        emg = emg.emg_rms;
    end

    load(fullfile(basepath, [basename, '.spec.mat']))
end


if graphics
      
    % open figure
    setMatlabGraphics(true)
    fh = figure;
    set(fh, 'WindowState', 'maximized');
    tlayout = [4, nstates + 3];
    th = tiledlayout(tlayout(1), tlayout(2));
    th.TileSpacing = 'tight';
    th.Padding = 'none';
    title(th, basename, 'interpreter', 'none', 'FontSize', 20)
    set(fh, 'DefaultAxesFontSize', 16);

    % spectrogram
    axh1 = nexttile(th, 1, [1, tlayout(2)]); cla; hold on
    plot_spec(spec, 'ch', 1, 'logfreq', false, 'saveFig', false,...
        'axh', axh1, 'xtime', 1)
    axis tight
    yLimit = ylim;
    xticks([]);
    xlabel('')

    % add spectrogram outliers
    % if ~isempty(specOtl)
    %     scatter(specOtl.idx, ones(length(specOtl.idx), 1) * yLimit(2), 30, 'filled', 'r')
    % end

    % hypnogram and emg
    axh2 = nexttile(th, tlayout(2) + 1, [1, tlayout(2)]); cla; hold on
    plot([1 : length(emg)], emg, 'k', 'LineWidth', 0.5)
    yLimit = ylim;
    plot_hypnogram('btimes', psd.bouts.btimes,...
        'clr', clr, 'axh', axh2, 'sstates', [1 : length(sstates)],...
        'yshift', 1)
    for istate = 1 : nstates
        badBouts{istate} = psd.bouts.otl(1).btimes{istate}(psd.bouts.otlIdx{istate}, :);
    end
    yshift = 1.05;
    plot_hypnogram('btimes', badBouts,...
        'clr', repmat({[0 0 0]}, length(sstates), 1),...
        'axh', axh2, 'sstates', [1 : length(sstates)], 'yshift', yshift)
    yLimit = ylim;
    ylim([yLimit(1), yLimit(2) * yshift])
    xval = [3600 : 3600 : length(emg)];
    xticks(xval);
    xticklabels(string(xval / 3600))
    xlabel('Time [Hr]')
    set(axh2, 'YTickMode', 'auto')
    set(axh2, 'YColor', 'k')
    yticks(axh2, [])
    linkaxes([axh1, axh2], 'x')
    axis tight
    ylabel(axh2, 'EMG')

    tilebias = tlayout(2) * 2;
    for istate = 1 : nstates
        axh = nexttile(th, istate + tilebias, [1, 1]); cla; hold on

        psdMat = psd.bouts.original{istate};
        badIdx = psd.bouts.otlIdx{istate};
        if ~isempty(psdMat)
            if any(badIdx)
                ph = plot(faxis, psdMat(badIdx, :), 'LineWidth', 0.5,...
                    'Color', [0.7 0.7 0.7]);
                ph = ph(1);
            end
            ph(2) = plot(faxis, mean(psdMat, 1), 'LineWidth', 3,...
                'Color', clr{istate});
            set(gca, 'YScale', 'log', 'XScale', 'log')
            xlabel('Frequency [Hz]')
            ylabel('PSD [mV^2/Hz]')
            title(axh, snames{istate})
            lgdTxt = sprintf('Outliers n=%d', length(badIdx));
            legend(ph, {lgdTxt, 'Mean'}, 'Location', 'southwest')
        end
    end

    % scree plot
    axh = nexttile(th, nstates + 1 + tilebias, [1, 1]);
    expl = psd.bouts.otl(1).expl;
    plot(cumsum(expl), 'LineWidth', 2);
    xlabel('PC');
    ylabel('Cum. Var. (%)');
    title('Scree Plot');
    hold on
    xlim([0 10])
    bh = bar(expl);
    bh.FaceColor = 'flat';
    bh.CData(1 : nfet, :) = repmat([0 0 0], nfet, 1);
    bh.CData(nfet + 1 : end, :) = repmat([1 0 0], size(bh.CData, 1) - nfet, 1);

    % state clusters in feature space
    axh = nexttile(th, nstates + 2 + tilebias, [1, 1]); cla; hold on
    pcMat = psd.bouts.otl(1).pc;
    badIdx = psd.bouts.otlIdx_glbl;
    clrIdx = psd.bouts.otl(1).clrIdx;
    sh = scatter3(axh, pcMat(:, 1), pcMat(:, 2), pcMat(:, 3),...
        30, clrIdx, 'filled');
    scatter3(pcMat(badIdx, 1), pcMat(badIdx, 2), pcMat(badIdx, 3),...
        50, 'd', 'LineWidth', 1, 'MarkerEdgeColor', 'k');
    view(-25, 25)
    grid minor
    sh.MarkerFaceAlpha = 0.6;
    xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
    title(axh, 'Bout PSD projection')

    % silhouette plot
    axh = nexttile(th, nstates + 3 + tilebias, [1, 1]); cla; hold on
    stateIdx = psd.bouts.otl(1).stateIdx;
    silhouette(pcMat, stateIdx, 'Euclidean');
    xlim([-1 1])
    yticklabels(snames)
    ylabel('')
    sh = get(gca, 'Children');
    sh.FaceColor = 'flat';
    sh.CData(~isnan(sh.YData), :) = clrIdx;
    title(axh, 'Separation Quality')
    yTickVals = yticks;
    for istate = 1 : nstates
        meanSil = mean(psd.bouts.otl(1).sil(stateIdx == istate));
        if ~isnan(meanSil)
            text(-0.75, yTickVals(istate), num2str(meanSil, '%.3f'));
        end
    end

    % final psd after cleaning
    tilebias = tlayout(2) * 3;
    
    for istate = 1 : length(sstates)
        axh = nexttile(th, istate + tilebias, [1, 1]); cla; hold on

        psdMat = psd.bouts.clean{istate};
        if ~isempty(psdMat)
            ph = plot(faxis, psdMat, 'LineWidth', 0.5,...
                'Color', [0.7 0.7 0.7]);
            ph = ph(1);
            ph(2) = plot(faxis, mean(psdMat, 1), 'LineWidth', 3,...
                'Color', clr{istate});
        end
        set(gca, 'YScale', 'log', 'XScale', 'log')
        xlabel('Frequency [Hz]')
        ylabel('PSD [mV^2/Hz]')
        title(axh, snames{istate})
        lgdTxt = sprintf('Clean n=%d', size(psdMat, 1));
        legend(ph, {lgdTxt, 'Mean'}, 'Location', 'southwest')
    end    

    % scree plot
    axh = nexttile(th, nstates + 1 + tilebias, [1, 1]); cla; hold on
    expl = psd.bouts.expl;
    plot(cumsum(expl), 'LineWidth', 2);
    xlabel('PC');
    ylabel('Cumulative Variance (%)');
    title('Scree Plot');
    hold on
    xlim([0 10])
    bh = bar(expl);
    bh.FaceColor = 'flat';
    bh.CData(1 : nfet, :) = repmat([0 0 0], nfet, 1);

    % state clusters in feature space
    axh = nexttile(th, nstates + 2 + tilebias, [1, 1]); cla; hold on
    pcMat = psd.bouts.pc;
    clrIdx = psd.bouts.clrIdx;
    sh = scatter3(axh, pcMat(:, 1), pcMat(:, 2), pcMat(:, 3),...
        30, clrIdx, 'filled');
    view(-25, 25)
    grid minor
    sh.MarkerFaceAlpha = 0.6;
    xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
    title(axh, 'Bout PSD projection')

    % silhouette plot
    axh = nexttile(th, nstates + 3 + tilebias, [1, 1]); cla; hold on
    stateIdx = psd.bouts.stateIdx;
    silhouette(pcMat, stateIdx, 'Euclidean');
    xlim([-1 1])
    yticklabels(snames)
    ylabel('')
    sh = get(gca, 'Children');
    sh.FaceColor = 'flat';
    sh.CData(~isnan(sh.YData), :) = clrIdx;
    title(axh, 'Silhouette of PCA Clusters')
    yTickVals = yticks;
    for istate = 1 : nstates
        meanSil = mean(psd.bouts.sil(stateIdx == istate));
        if ~isnan(meanSil)
            text(-0.75, yTickVals(istate), num2str(meanSil, '%.3f'));
        end
    end

    if flgSaveFig
        figpath = fullfile(basepath, 'graphics', 'sleepState');
        mkdir(figpath)
        figname = fullfile(figpath, [basename, '_', namePrefix]);
        savefig(fh, figname, 'compact')
    end

end


%%% ALT 2 - all iterations
% if graphics
% 
%     setMatlabGraphics(true)
%     fh = figure;
%     set(fh, 'WindowState','maximized');
%     tlayout = [nrep + 1, nstates + 3];
%     th = tiledlayout(tlayout(1), tlayout(2));
%     th.TileSpacing = 'tight';
%     th.Padding = 'none';
%     title(th, basename, 'interpreter', 'none', 'FontSize', 20)
%     set(fh, 'DefaultAxesFontSize', 16);
% 
%     for irep = 1 : nrep
%         tilebias = (irep - 1) * (nstates + 3);
% 
%         for istate = 1 : nstates
%             axh = nexttile(th, istate + tilebias, [1, 1]); cla; hold on
% 
%             psdMat = psd.bouts.otl(irep).psd_bouts{istate};
%             badIdx = psd.bouts.otl(irep).badIdx{istate};
%             if ~isempty(psdMat)
%                 if any(badIdx)
%                     ph = plot(faxis, psdMat(badIdx, :), 'LineWidth', 0.5,...
%                         'Color', [0.7 0.7 0.7]);
%                 end
%                 ph = plot(faxis, mean(psdMat, 1), 'LineWidth', 3,...
%                     'Color', clr{istate});
%                 set(gca, 'YScale', 'log', 'XScale', 'log')
%                 xlabel('Frequency [Hz]')
%                 ylabel('PSD [mV^2/Hz]')
%                 title(axh, snames{istate})
%                 lgdTxt = sprintf('Ourliers n=%d', sum(badIdx));
%                 legend({'Mean', lgdTxt}, 'Location', 'southwest')
%             end
%         end
% 
%         % scree plot
%         axh = nexttile(th, nstates + 1 + tilebias, [1, 1]);
%         expl = psd.bouts.otl(irep).expl;
%         plot(cumsum(expl), 'LineWidth', 2);
%         xlabel('PC');
%         ylabel('Cum. Var. (%)');
%         title('Scree Plot');
%         hold on
%         xlim([0 10])
%         bh = bar(expl);
%         bh.FaceColor = 'flat';
%         bh.CData(1 : nfet, :) = repmat([0 0 0], nfet, 1);
%         bh.CData(nfet + 1 : end, :) = repmat([1 0 0], size(bh.CData, 1) - nfet, 1);
% 
%         % state clusters in feature space
%         axh = nexttile(th, nstates + 2 + tilebias, [1, 1]); cla; hold on
%         pcMat = psd.bouts.otl(irep).pc;
%         badIdx = vertcat(psd.bouts.otl(irep).badIdx{:});
%         clrIdx = psd.bouts.otl(irep).clrIdx;      
%         sh = scatter3(axh, pcMat(:, 1), pcMat(:, 2), pcMat(:, 3),...
%             30, clrIdx, 'filled');
%         scatter3(pcMat(badIdx, 1), pcMat(badIdx, 2), pcMat(badIdx, 3),...
%             50, 'd', 'LineWidth', 1, 'MarkerEdgeColor', 'k');
%         view(-25, 25)
%         grid minor
%         sh.MarkerFaceAlpha = 0.6;
%         xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
%         title(axh, 'Bout PSD projection')
% 
%         % silhouette plot
%         axh = nexttile(th, nstates + 3 + tilebias, [1, 1]); cla; hold on
%         stateIdx = psd.bouts.otl(irep).stateIdx;
%         silhouette(pcMat, stateIdx, 'Euclidean');
%         xlim([-1 1])
%         yticklabels(snames)
%         ylabel('')
%         sh = get(gca, 'Children');
%         sh.FaceColor = 'flat';
%         sh.CData(~isnan(sh.YData), :) = clrIdx;
%         title(axh, 'Separation Quality')
%         yTickVals = yticks;
%         for istate = 1 : nstates
%             meanSil = mean(psd.bouts.otl(irep).sil(stateIdx == istate));
%             if ~isnan(meanSil)
%                 text(-0.75, yTickVals(istate), num2str(meanSil, '%.3f'));
%             end
%         end
%     end
%     
%     % final psd after cleaning
%     tilebias = nrep * (nstates + 3);
%     for istate = 1 : length(sstates)
%         axh = nexttile(th, istate + tilebias, [1, 1]); cla; hold on
% 
%         psdMat = psd.bouts.clean{istate};
%         if ~isempty(psdMat)
%             ph = plot(faxis, psdMat, 'LineWidth', 0.5,...
%                 'Color', [0.7 0.7 0.7]);
%             ph = plot(faxis, mean(psdMat, 1), 'LineWidth', 3,...
%                 'Color', clr{istate});
%         end
%         set(gca, 'YScale', 'log', 'XScale', 'log')
%         xlabel('Frequency [Hz]')
%         ylabel('PSD [mV^2/Hz]')
%         title(axh, snames{istate})
%         lgdTxt = sprintf('All n=%d', size(psdMat, 1));
%         legend({'Mean', lgdTxt}, 'Location', 'southwest')
%     end    
% 
%     % scree plot
%     axh = nexttile(th, nstates + 1 + tilebias, [1, 1]); cla; hold on
%     expl = psd.bouts.expl;
%     plot(cumsum(expl), 'LineWidth', 2);
%     xlabel('PC');
%     ylabel('Cumulative Variance (%)');
%     title('Scree Plot');
%     hold on
%     xlim([0 10])
%     bh = bar(expl);
%     bh.FaceColor = 'flat';
%     bh.CData(1 : nfet, :) = repmat([0 0 0], nfet, 1);
% 
%     % state clusters in feature space
%     axh = nexttile(th, nstates + 2 + tilebias, [1, 1]); cla; hold on
%     pcMat = psd.bouts.pc;
%     clrIdx = psd.bouts.clrIdx;
%     sh = scatter3(axh, pcMat(:, 1), pcMat(:, 2), pcMat(:, 3),...
%         30, clrIdx, 'filled');
%     view(-25, 25)
%     grid minor
%     sh.MarkerFaceAlpha = 0.6;
%     xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
%     title(axh, 'Bout PSD projection')
% 
%     % silhouette plot
%     axh = nexttile(th, nstates + 3 + tilebias, [1, 1]); cla; hold on
%     stateIdx = psd.bouts.stateIdx;
%     silhouette(pcMat, stateIdx, 'Euclidean');
%     xlim([-1 1])
%     yticklabels(snames)
%     ylabel('')
%     sh = get(gca, 'Children');
%     sh.FaceColor = 'flat';
%     sh.CData(~isnan(sh.YData), :) = clrIdx;
%     title(axh, 'Silhouette of PCA Clusters')
%     yTickVals = yticks;
%     for istate = 1 : nstates
%         meanSil = mean(psd.bouts.sil(stateIdx == istate));
%         if ~isnan(meanSil)
%             text(-0.75, yTickVals(istate), num2str(meanSil, '%.3f'));
%         end
%     end
% 
%     if flgSaveFig
%         figpath = fullfile(basepath, 'graphics', 'sleepState');
%         mkdir(figpath)
%         figname = fullfile(figpath, [basename, '_', namePrefix]);
%         savefig(fh, figname, 'compact')
%     end
% 
% end
% 


end


% EOF