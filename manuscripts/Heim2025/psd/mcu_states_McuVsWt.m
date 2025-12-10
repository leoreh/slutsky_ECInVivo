% mcu_states_McuVsWt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize and load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepaths = [mcu_basepaths('wt_bsl'); mcu_basepaths('mcu_bsl')];
nfiles = length(basepaths);
mnames = get_mname(basepaths);

% load data
v = basepaths2vars('basepaths', basepaths, 'vars',...
    ["sleep_labelsMan"; "sleep_states"]);

% state params
cfg = as_loadConfig();
sstates = [1, 4, 5];
nstates = length(sstates);
clr = cfg.colors(sstates);
snames = cfg.names(sstates);

% bout params
interDur = 4;
minDur = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update sleep_states struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% go over files
for ifile = 1 : nfiles

    % file
    basepath = basepaths{ifile};
    cd(basepath)
    [~, basename] = fileparts(basepath);

    % insert manual labels to auto labels
    manIdx = v(ifile).labels < cfg.nstates + 2;
    ss = v(ifile).ss;
    labels = ss.labels;
    labels(manIdx) = v(ifile).labels(manIdx);

    % calculate state bouts
    [bouts] = as_bouts('labels', labels,...
        'minDur', minDur, 'interDur', interDur);

    % save bckup of original ss struct
    ssOrig = v(ifile).ss;
    mkdir(fullfile(basepath, 'ss'))
    ctime = datetime('now', 'Format', 'yyMMdd_HHmmss');
    bkupName = [basename, '.sleep_states', '_', char(ctime), '.mat'];
    bkupFile = fullfile(basepath, 'ss', bkupName);
    save(bkupFile, 'ss')

    % save updated ss struct
    ssFile = fullfile(basepath, [basename, '.sleep_states.mat']);
    ss.labels = labels;
    ss.bouts = bouts;
    save(ssFile, 'ss')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reload data
clear basepaths
basepaths{1} = mcu_basepaths('wt_bsl');
basepaths{2} = mcu_basepaths('mcu_bsl');

% initialize
boutLen = nan(2, nfiles, length(sstates));
prctDur = nan(2, nfiles, length(sstates));

for igrp = 1 : 2
    v = basepaths2vars('basepaths', basepaths{igrp}, 'vars',...
        ["sleep_states"]);
    nfiles = length(basepaths{igrp});
    ss = catfields([v(:).ss], 1, true);

    for ifile = 1 : nfiles
        boutLen(igrp, ifile, :) = mean(cell2padmat(ss.bouts.boutLen(ifile, sstates), 2), 'omitnan');
        prctDur(igrp, ifile, :) = ss.bouts.prctDur(ifile, sstates);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open figure
fh = figure;
setMatlabGraphics(true);
set(fh, 'WindowState','maximized');
tlayout = [2, length(sstates)];
th = tiledlayout(tlayout(1), tlayout(2));
th.TileSpacing = 'tight';
th.Padding = 'none';

% Plot bout length for each state
tilebias = 0;
for istate = 1 : length(sstates)
    axh = nexttile(th, tilebias + istate); hold on; cla;
    stateIdx = sstates(istate);
    dataMat = boutLen(:, :, istate)';
    plot_boxMean('axh', axh, 'dataMat', dataMat, 'allPnt', false,...
        'plotType', 'bar')
    ylabel('Bout Length (s)');
    title([snames{istate}]);
    ylim([0, ceil(max(dataMat(:)))]);
    xlim([0.5 2.5])
    xticks([1, 2])
    xticklabels({'WT', 'MCU-KO'});
end


% Plot bout length for each state
tilebias = 3;
for istate = 1 : length(sstates)
    axh = nexttile(th, tilebias + istate); hold on; cla;
    stateIdx = sstates(istate);
    dataMat = prctDur(:, :, istate)';
    plot_boxMean('axh', axh, 'dataMat', dataMat, 'allPnt', false,...
        'plotType', 'bar')
    ylabel('State Duration (%)');
    title([snames{istate}]);
    ylim([0, ceil(max(dataMat(:)))]);
    xlim([0.5 2.5])
    xticks([1, 2])
    xticklabels({'WT', 'MCU-KO'});
end


