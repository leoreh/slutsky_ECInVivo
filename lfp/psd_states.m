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
%                   which to calculate stateEpochs, i.e. the time window for
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
%   stateEpochs     cell of n x 2 mats. if empty will calculate from
%                   ss.labels or from emg_labels
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
%                   per epoch. much simpler this way. bkup exists in temp
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
addOptional(p, 'stateEpochs', []);
addOptional(p, 'ftarget', [0.5 : 0.5 : 100], @isnumeric);
addOptional(p, 'emgThr', [50], @isnumeric);
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
stateEpochs     = p.Results.stateEpochs;
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

% flg ss epochs or emg epochs
if flgEmg
    flgEmg = true;
    namePrefix = 'psdEmg';

    % state params
    clr = {[150, 70, 55] / 255, [200, 170, 100] / 255};
    snames = {'High-EMG', 'Low-EMG'};
    sstates = [1, 2];

    % emg data
    emg = load(sleepfile, 'emg_rms');
    emg = emg.emg_rms;

    % state epoch duration limits
    minDur = [20];
    interDur = 0;

else
    flgEmg = false;
    namePrefix = 'psd';

    % state params
    cfg = as_loadConfig();
    if isempty(sstates)
        sstates = 1 : cfg.nstates;
    end
    clr = cfg.colors(sstates);
    snames = cfg.names(sstates);

    % state epoch duration limits
    minDur = [20, 5, 5, 20, 10, 5];;
    interDur = 0;

end

% check if already analyzed
psdfile = fullfile(basepath, [basename, '.', namePrefix, '.mat']);
if exist(psdfile, 'file') && ~forceA
    load(psdfile)
    return
end

% update params
nstates = length(sstates);

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
% calculate state epochs if not provided
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get stateEpochs
if isempty(stateEpochs)
    if ~flgEmg

        % get AS labels
        labels = v.ss.labels;
    else

        % find threshold to separate the bimodal distribution of emg
        if isempty(emgThr)
            [~, cents] = kmeans(emg(:), 2);
            emgThr = mean(cents);
        end

        % create EMG labels
        labels = double(emg > emgThr);
        labels(emg < emgThr) = 2;
    end

    % re-calc state epochs from labels
    [stateEpochs, ~] = as_epochs('labels', labels(timeWin(1) : timeWin(2)),...
        'minDur', minDur, 'interDur', interDur, 'rmOtl', true);
    stateEpochs = stateEpochs(sstates);

    % add timeWin(1) to each epoch
    stateEpochs = cellfun(@(x) round(x + timeWin(1) - 1),...
        stateEpochs, 'uni', false);

end

% calc epoch stats
epochStats.epLen = cellfun(@(x) (diff(x')'), stateEpochs, 'UniformOutput', false);
epochStats.nepochs = cellfun(@length, epochStats.epLen);
epochStats.totDur = cellfun(@sum, epochStats.epLen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc psd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, faxis, psd_epochs] = calc_psd('sig',...
    sig, 'bins', stateEpochs,...
    'fs', fs, 'ftarget', ftarget, 'graphics', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove outlier epochs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% outliers are detected in two iterations, first by pc1 and then by
% mahalnobis distance. the pca is recalculated between iterations. this
% was determined by empiric examination of the data (march 2024).
clear otl
psd.epochs.original = psd_epochs;
tmp_epochs = psd_epochs;

% project epoch PSDs on feature space
[~, pc, ~, ~, expl] = pca(vertcat(tmp_epochs{:}), 'NumComponents', nfet);

% create state and color indices for each epoch when epochs are
% concatenated (e.g., for PCs)
stateIdx = []; clrIdx = [];
for istate = 1 : nstates
    stateIdx = [stateIdx; istate * ones(epochStats.nepochs(istate), 1)];
    clrIdx = [clrIdx; repmat(clr{istate}, epochStats.nepochs(istate), 1)];
end

% detect outliers
nrep = 3;
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
    end

    % organize outliers struct as an array (of iterations)
    otl(irep).pc = pc;
    otl(irep).expl = expl;
    otl(irep).stateIdx = stateIdx;
    otl(irep).clrIdx = clrIdx;
    otl(irep).badIdx = badIdx;
    otl(irep).stateEpochs = stateEpochs;
    otl(irep).epochStats = epochStats;
    otl(irep).psd_epochs = tmp_epochs;
    otl(irep).sil = silhouette(pc, stateIdx, 'Euclidean');

    % remove bad indices from epochs and recalculate epoch stats
    for istate = 1 : nstates
        if any(badIdx{istate})
            stateEpochs{istate}(badIdx{istate}, :) = [];
            tmp_epochs{istate}(badIdx{istate}, :) = [];
        end
    end
    epochStats.epLen = cellfun(@(x) (diff(x')'), stateEpochs, 'UniformOutput', false);
    epochStats.nepochs = cellfun(@length, epochStats.epLen);
    epochStats.totDur = cellfun(@sum, epochStats.epLen);

    % recalculate pca
    [~, pc, ~, ~, expl] = pca(vertcat(tmp_epochs{:}), 'NumComponents', nfet);

    % create state and color indices
    stateIdx = []; clrIdx = [];
    for istate = 1 : nstates
        stateIdx = [stateIdx; istate * ones(epochStats.nepochs(istate), 1)];
        clrIdx = [clrIdx; repmat(clr{istate}, epochStats.nepochs(istate), 1)];
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize psd struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% organize epochs
psd.epochs.clean = tmp_epochs;
psd.epochs.stateEpochs = stateEpochs;
psd.epochs.epochStats = epochStats;
psd.epochs.otl = otl;
psd.epochs.pc = pc;
psd.epochs.expl = expl;
psd.epochs.sil = silhouette(pc, stateIdx, 'Euclidean');
psd.epochs.stateIdx = stateIdx;
psd.epochs.clrIdx = clrIdx;

% info
psd.info.sstates = sstates;
psd.info.snames = snames;
psd.info.clr = clr; 
psd.info.faxis = faxis;
psd.info.emgThr = emgThr;
psd.info.timeWin = timeWin;
psd.info.ch = ch;
psd.info.spkgrpIdx = find(cellfun(@(x) all(ismember(ch, x)), spkgrp));
psd.info.runtime = datetime(now, 'ConvertFrom', 'datenum');

% calculate vars using cleaned epochs
for istate = 1 : length(sstates)
    if ~isempty(psd.epochs.clean{istate})
        psd.psd(istate, :) = mean(psd.epochs.clean{istate}, 1);

        % calc power in specific frequency bands
        % per epoch
        [psd.bands.epochs{istate}, psd.bands.info] = calc_bands('psdData',...
            psd.epochs.clean{istate}, 'freq', faxis, 'flgNormBand', false);
        % across epochs
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
if graphics

    setMatlabGraphics(true)
    fh = figure;
    set(fh, 'WindowState','maximized');
    tlayout = [nrep + 1, nstates + 3];
    th = tiledlayout(tlayout(1), tlayout(2));
    th.TileSpacing = 'tight';
    th.Padding = 'none';
    title(th, basename, 'interpreter', 'none')

    for irep = 1 : nrep
        tilebias = (irep - 1) * (nstates + 3);

        for istate = 1 : nstates
            axh = nexttile(th, istate + tilebias, [1, 1]); cla; hold on

            psdMat = psd.epochs.otl(irep).psd_epochs{istate};
            badIdx = psd.epochs.otl(irep).badIdx{istate};
            if ~isempty(psdMat)
                if any(badIdx)
                    ph = plot(faxis, psdMat(badIdx, :), 'LineWidth', 0.5,...
                        'Color', [0.7 0.7 0.7]);
                end
                ph = plot(faxis, mean(psdMat, 1), 'LineWidth', 3,...
                    'Color', clr{istate});
                set(gca, 'YScale', 'log', 'XScale', 'log')
                xlabel('Frequency [Hz]')
                ylabel('PSD [mV^2/Hz]')
                title(axh, snames{istate})
                lgdTxt = sprintf('Ourliers n=%d', sum(badIdx));
                legend({'Mean', lgdTxt})
            end
        end

        % scree plot
        axh = nexttile(th, nstates + 1 + tilebias, [1, 1]);
        expl = psd.epochs.otl(irep).expl;
        plot(cumsum(expl), 'LineWidth', 2);
        xlabel('PC');
        ylabel('Cumulative Variance (%)');
        title('Scree Plot');
        hold on
        xlim([0 10])
        bh = bar(expl);
        bh.FaceColor = 'flat';
        bh.CData(1 : nfet, :) = repmat([0 0 0], nfet, 1);
        bh.CData(nfet + 1 : end, :) = repmat([1 0 0], size(bh.CData, 1) - nfet, 1);

        % state clusters in feature space
        axh = nexttile(th, nstates + 2 + tilebias, [1, 1]); cla; hold on
        pcMat = psd.epochs.otl(irep).pc;
        badIdx = vertcat(psd.epochs.otl(irep).badIdx{:});
        clrIdx = psd.epochs.otl(irep).clrIdx;      
        sh = scatter3(axh, pcMat(:, 1), pcMat(:, 2), pcMat(:, 3),...
            30, clrIdx, 'filled');
        scatter3(pcMat(badIdx, 1), pcMat(badIdx, 2), pcMat(badIdx, 3),...
            50, 'd', 'LineWidth', 1, 'MarkerEdgeColor', 'k');
        view(-25, 25)
        grid minor
        sh.MarkerFaceAlpha = 0.6;
        xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
        title(axh, 'Epoch PSD projection')

        % silhouette plot
        axh = nexttile(th, nstates + 3 + tilebias, [1, 1]); cla; hold on
        stateIdx = psd.epochs.otl(irep).stateIdx;
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
            meanSil = mean(psd.epochs.otl(irep).sil(stateIdx == istate));
            if ~isnan(meanSil)
                text(-0.75, yTickVals(istate), num2str(meanSil, '%.3f'));
            end
        end
    end
    
    % final psd after cleaning
    tilebias = nrep * (nstates + 3);
    for istate = 1 : length(sstates)
        axh = nexttile(th, istate + tilebias, [1, 1]); cla; hold on

        psdMat = psd.epochs.clean{istate};
        if ~isempty(psdMat)
            ph = plot(faxis, psdMat, 'LineWidth', 0.5,...
                'Color', [0.7 0.7 0.7]);
            ph = plot(faxis, mean(psdMat, 1), 'LineWidth', 3,...
                'Color', clr{istate});
        end
        set(gca, 'YScale', 'log', 'XScale', 'log')
        xlabel('Frequency [Hz]')
        ylabel('PSD [mV^2/Hz]')
        title(axh, snames{istate})
        lgdTxt = sprintf('All n=%d', size(psdMat, 1));
        legend({'Mean', lgdTxt})
    end    

    % scree plot
    axh = nexttile(th, nstates + 1 + tilebias, [1, 1]); cla; hold on
    expl = psd.epochs.expl;
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
    pcMat = psd.epochs.pc;
    clrIdx = psd.epochs.clrIdx;
    sh = scatter3(axh, pcMat(:, 1), pcMat(:, 2), pcMat(:, 3),...
        30, clrIdx, 'filled');
    view(-25, 25)
    grid minor
    sh.MarkerFaceAlpha = 0.6;
    xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
    title(axh, 'Epoch PSD projection')

    % silhouette plot
    axh = nexttile(th, nstates + 3 + tilebias, [1, 1]); cla; hold on
    stateIdx = psd.epochs.stateIdx;
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
        meanSil = mean(psd.epochs.sil(stateIdx == istate));
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

end


% EOF