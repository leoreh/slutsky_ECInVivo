

mnames = ["lh96"; "lh107"; "lh122"; "lh142"; "lh100"; "lh99"];
mnames = ["lh132"; "lh133"; "lh134"; "lh136"; "lh140"];

for iname = 1 : length(mnames)
    mname = mnames{iname};
    mcu_psd(mnames{iname})
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

queryStr = 'lh142';
basepaths = mcu_basepaths(queryStr);

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

    % load lfp signal
    sleepfile = fullfile(basepath, [basename, '.sleep_sig.mat']);
    sig = load(sleepfile, 'eeg');
    sig = sig.eeg;
    load(sleepfile, 'fs');

    % get spectrogram outliers
    otl = get_specOutliers('basepath', basepath, 'saveVar', true,...
        'flgCalc', false, 'flgForce', false, 'graphics', false,...
        'stdThr', 5);

    % calc psd according to as state separation
    psd = psd_states('basepath', basepath, 'sstates', sstates,...
        'sig', sig, 'fs', fs, 'saveVar', saveVar,...
        'graphics', true, 'forceA', true, 'ftarget', ftarget,...
        'emgThr', [], 'flgEmg', false);

    % calc psd according to emg state separation
    psd = psd_states('basepath', basepath, 'sstates', [1, 2],...
        'sig', sig, 'fs', fs, 'saveVar', saveVar,...
        'graphics', true, 'forceA', true, 'ftarget', ftarget,...
        'emgThr', [], 'flgEmg', true);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bout stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grpnames{1} = ["lh96"; "lh107"; "lh122"; "lh142"; "lh100"; "lh99"];
grpnames{2} = ["lh132"; "lh133"; "lh134"; "lh136"; "lh140"; "lh137"];

% ALT 2 - bouts from labels w/o processing
sstates = [1, 4, 5];
minDur = [10, 5, 5, 10, 5, 5];
minDur = minDur(sstates);
interDur = 2;
nbins = 1;
clear grpStats
for igrp = 1 : 2
    mnames = grpnames{igrp};

    clear boutStats tmpStats mStats
    for imouse = 1 : length(mnames)
        queryStr = mnames{imouse};
        basepaths = mcu_basepaths(queryStr);
        if ~strcmp(queryStr, 'lh137')
            basepaths = basepaths([1, length(basepaths)]);
        end

        % load data
        varsFile = ["sleep_states"; "datInfo"; "session"];
        varsName = ["ss"; "datInfo"; "session"];
        xlsname = 'D:\Google Drive\PhD\Slutsky\Data Summaries\sessionList.xlsx';
        [v, basepaths] = getSessionVars('basepaths', basepaths, 'varsFile', varsFile,...
            'varsName', varsName, 'pcond', ["tempflag"], 'ncond', [""],...
            'xlsname', xlsname);
        nfiles = length(basepaths);

        % get bouts stats
        for ifile = 1 : nfiles

            % process labels
            labels = v(ifile).ss.labels;
            %             labels(labels == 3) = 4;
            %             labels(labels == 6) = 5;
            %             labels(labels == 2) = 1;

            [~, tmpStats(imouse, ifile)] = as_bouts('labels', labels,...
                'minDur', minDur, 'interDur', interDur,...
                'sstates', sstates, 'nbins', nbins, 'graphics', false);
        end
        if isempty(tmpStats(imouse, 2).boutLen)
            tmpStats(imouse, 2).boutLen = {};
        end
        mStats(imouse) = catfields(tmpStats(imouse, :), 'addim', true);
    end

    % organize
    grpStats(igrp) = catfields(mStats, 'addim', true);

end

boutStats = catfields(grpStats, 3, true);

% to prism
istate = 3;
prismData = squeeze(boutStats.prctDur(:, istate, :, :));
prismData = reshape(permute(prismData, [1, 3, 2]), nbins, []);

ibin = 1;
istate = 1;
prismData = squeeze(boutStats.prctDur(ibin, istate, :, :));
x = prismData([1, 2], :);
x([1, 2], 7 : 12) = prismData([3, 4], :);



% state params
cfg = as_loadConfig;
nstates = length(sstates);
snames = cfg.names(sstates);
clr = v(1).psd.info.clr;
clr = cfg.colors(sstates);

% open figure
fh = figure;
setMatlabGraphics(true)
set(fh, 'WindowState','maximized');
tlayout = [3, nstates];
th = tiledlayout(tlayout(1), tlayout(2));
th.TileSpacing = 'tight';
th.Padding = 'none';
% title(th, mname, 'interpreter', 'none')

tilebias = 0;
for istate = 1 : nstates
    axh = nexttile(th, istate + tilebias, [1, 1]); cla
    dataMat = squeeze(boutStats.prctDur(:, istate, :));
    plot_boxMean('dataMat', dataMat, 'clr', clr{istate},...
        'plotType', 'bar', 'axh', axh, 'allpnts', true)
    ylabel('Duration (%)')
end

tilebias = nstates;
for istate = 1 : nstates
    axh = nexttile(th, istate + tilebias, [1, 1]); cla
    dataMat = squeeze(boutStats.nbouts(:, istate, :));
    plot_boxMean('dataMat', dataMat, 'clr', clr{istate},...
        'plotType', 'bar', 'axh', axh, 'allpnts', true)
    ylabel('No. Bouts')
end

tilebias = nstates * 2;
for istate = 1 : nstates
    axh = nexttile(th, istate + tilebias, [1, 1]); cla
    dataMat = cellfun(@mean, squeeze(boutStats.boutLen(:, istate, :)), 'uni', true);
    plot_boxMean('dataMat', dataMat, 'clr', clr{istate},...
        'plotType', 'bar', 'axh', axh, 'allpnts', true)
    ylabel('Mean Bout Length (s)')
end





