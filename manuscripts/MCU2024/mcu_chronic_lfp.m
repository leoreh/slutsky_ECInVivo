

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = as_loadConfig;
sstates = [1, 4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load files of single mouse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mname = 'lh134';

varsFile = ["fr"; "sr"; "sleep_states";...
    "datInfo"; "session"; "units"; "psd"];
varsName = ["fr"; "sr"; "ss"; "datInfo"; "session";...
    "units"; "psd"];
xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
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
% hypnogram during baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mcu baseline basepaths
basepaths = {...
    'F:\Data\lh132\lh132_230413_094013',...
    'F:\Data\lh133\lh133_230413_094013',...
    'F:\Data\lh134\lh134_230504_091744',...
    'F:\Data\lh136\lh136_230519_090043',...
    'F:\Data\lh140\lh140_230619_090023'};

% wt baseline basepaths
basepaths = {...
    'F:\Data\lh96\lh96_220120_090157',...
    'F:\Data\lh107\lh107_220518_091200',...
    'F:\Data\lh122\lh122_221223_092656',...
    'F:\Data\lh142\lh142_231005_091832'};

% load data
varsFile = ["sleep_states"; "datInfo"; "session"];
varsName = ["ss"; "datInfo"; "session"];
xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
[v, basepaths] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
    'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
    'xlsname', xlsname);

% params
nfiles = length(basepaths);

% % manually check states
% for ifile = 1 : nfiles
%     cd(basepaths{ifile})
%     [~, basename] = fileparts(basepaths{ifile})
%     sSig = load([basename, '.sleep_sig.mat']);
%     AccuSleep_viewer(sSig, v(ifile).ss.labels, [])
% end

% % add timebins to session
% for ifile = 1 : nfiles
%     cd(basepaths{ifile})
%     [timebins, timepnt] = metaInfo_timebins('reqPnt', 6 * 60 * 60, 'nbins', 4);
%     timebins / 60 / 60
% end

% % plot hynogram
% fh = figure;
% for ifile = 1 : nfiles
% 
%     sb = subplot(nfiles, 1, ifile);
%     plot_hypnogram('stateEpochs', v(ifile).ss.stateEpochs, 'axh', sb)
% 
% end

% plot state duration in timebins
nbins = 2;
sstates = [1, 4, 5];
clear prctDur epLen
for ifile = 1 : nfiles
    
    % timebins for calculating stats
    if nbins == 4
        timebins = v(ifile).session.general.timebins;
    else
        nlabels = length(v(ifile).ss.labels);
        timebins = n2chunks('n', nlabels, 'nchunks', nbins);
    end

    cd(basepaths{ifile})
    [totDur, prctDur(ifile, :, :), epLen(ifile, :)] = as_plotZT('nwin', 4,...
        'sstates', sstates, 'ss', v(ifile).ss,...
        'timebins', timebins, 'graphics', false);
end

% ---------- reorganize for prism
% state duration
istate = 1;
prismData = prctDur(:, :, istate)';

% epoch length
% concatenated across mice. each column corresponds to timebine.
istate = 1;
prismData = vertcat(epLen{:, istate});

% nepochs 
ibin = 2;
prismData = cellfun(@(x) sum(~isnan(x(:, ibin))), epLen, 'UniformOutput', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare state separation across time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%
% multiple ways to show renormalization of state-dependent psd:
% (1) delta-theta ratio across time
% (2) difference between delta (and other bands) between low and high emg
% (3) accusleep recovers its ability to classify states
% (4) plot emg vs. delta theta for epochs accusleep classified as AW and
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
        prct = 70;
        highIdx = emg_rms > prctile(emg_rms, prct);
        lowIdx = emg_rms < prctile(emg_rms, 100 - prct);

        % 3. calc PCs from spec
        [coeff, score, ~] = pca(spec, 'NumComponents', 3);

        % 4. measure distance
        sidx = highIdx | lowIdx;
        [lRat, iDist, mDist] = cluDist(score(sidx, :), highIdx(sidx));
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
            scatter3(score(highIdx, 1), score(highIdx, 2), score(highIdx, 3),...
                sz, 'r', 'filled', 'MarkerFaceAlpha', falpha, 'MarkerEdgeColor', 'none');
            hold on;
            scatter3(score(lowIdx, 1), score(lowIdx, 2), score(lowIdx, 3),...
                sz, 'b', 'filled', 'MarkerFaceAlpha', falpha, 'MarkerEdgeColor', 'none');
            legend('High-EMG', 'Low-EMG');
            xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
            title([basename, '. mDist = ', d])
        end


    end



end

