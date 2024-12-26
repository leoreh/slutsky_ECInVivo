% mcu_states_McuVsWt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize and load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basepaths = [mcu_sessions('mcu_bsl')];
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
    ss.labels(manIdx) = v(ifile).labels(manIdx);
    
    % calculate state epochs

    % auto
    [stateEpochs, epochStats] = as_epochs('labels', labelsAuto(curIdx),...
        'minDur', minDur, 'interDur', interDur);

    % manual
    [~, epochStats(ifile, 2)] = as_epochs('labels', labelsMan(curIdx),...
        'minDur', minDur, 'interDur', interDur);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

epLen = nan(2, nfiles, length(sstates));
for ifile = 1 : nfiles
    for icur = 1 : 2
        epLen(icur, ifile, :) = mean(cell2padmat(epochStats(ifile, icur).epLen(sstates), 2), 'omitnan');
        prctDur(icur, ifile, :) = epochStats(ifile, icur).prctDur(sstates);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inspect data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sSig = load([basename, '.sleep_sig.mat']);
% AccuSleep_viewer(sSig, v(ifile).labels, [])
% AccuSleep_viewer(sSig, labelsAuto, [])

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
title(th, basename, 'interpreter', 'none');

% Plot epoch length for each state
tilebias = 0;
for istate = 1 : length(sstates)
    axh = nexttile(th, tilebias + istate); hold on; cla;
    stateIdx = sstates(istate);
    dataMat = epLen(:, :, istate)';
    ph = plot([1, 2], dataMat);
    ylabel('Epoch Length (s)');
    title([snames{istate}]);
    ylim([0, ceil(max(dataMat(:)))]);
    xlim([0.5 2.5])
    xticks([1, 2])
    xticklabels({'Auto', 'Man'});
end
legend(mnames)

% Plot epoch length for each state
tilebias = 3;
for istate = 1 : length(sstates)
    axh = nexttile(th, tilebias + istate); hold on; cla;
    stateIdx = sstates(istate);
    dataMat = prctDur(:, :, istate)';
    ph = plot([1, 2], dataMat);
    ylabel('State Duration (%)');
    title([snames{istate}]);
    ylim([0, ceil(max(dataMat(:)))]);
    xlim([0.5 2.5])
    xticks([1, 2])
    xticklabels({'Auto', 'Man'});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot confusion matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % confusion matrix (relative to manual labels)
% [ss.netPrecision, ss.netRecall] = as_cm(labelsMan, labelsAuto,...
%     'saveFig', false);