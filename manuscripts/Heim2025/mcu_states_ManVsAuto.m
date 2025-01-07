% mcu_states_ManVsAuto


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize and load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% files with manual curation
basepaths = {...
    'E:\Data\lh132\lh132_230413_094013',...
    'E:\Data\lh133\lh133_230413_094013',...
    'E:\Data\lh96\lh96_220120_090157',...
    'E:\Data\lh142\lh142_231005_091832',...
    'E:\Data\lh107\lh107_220518_091200',...
    };
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
% find indices of manual curation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize 
manIdx = nan(nfiles, 2);
clear boutStats

% go over files
for ifile = 1 : nfiles
    
    % file
    basepath = basepaths{ifile};
    cd(basepath)
    [~, basename] = fileparts(basepath);

    % get labels
    labelsAuto = v(ifile).ss.labels_net;
    labelsMan = v(ifile).labels;

    % ---------------------------------------------------------------------
    % detect undefined labels  
    diffs = diff([0; labelsMan ~= cfg.nstates + 2; 0]);  % Add padding to detect edges
    startIdx = find(diffs == 1);
    endIdx = find(diffs == -1) - 1;
    
    % manual override for lh142
    endIdx = min([endIdx, 6 * 60 * 60]);

    % filter out undefined segments 
    validSegments = arrayfun(@(s, e) any(labelsMan(s:e) ~= cfg.nstates + 2),...
        startIdx, endIdx);

    % extract start and end indices of valid segments
    manIdx(ifile, :) = [startIdx(validSegments), endIdx(validSegments)];
    curIdx = manIdx(ifile, 1) : manIdx(ifile, 2);     % manually curated indices 
    
    % ---------------------------------------------------------------------
    % calculate boutStats for Man and Auto

    % auto
    bouts(ifile, 1) = as_bouts('labels', labelsAuto(curIdx),...
        'minDur', minDur, 'interDur', interDur);

    % manual
    bouts(ifile, 2) = as_bouts('labels', labelsMan(curIdx),...
        'minDur', minDur, 'interDur', interDur);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

boutLen = nan(2, nfiles, length(sstates));
for ifile = 1 : nfiles
    for icur = 1 : 2
        boutLen(icur, ifile, :) = mean(cell2padmat(bouts(ifile, icur).boutLen(sstates), 2), 'omitnan');
        prctDur(icur, ifile, :) = bouts(ifile, icur).prctDur(sstates);
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

% Plot bout length for each state
tilebias = 0;
for istate = 1 : length(sstates)
    axh = nexttile(th, tilebias + istate); hold on; cla;
    stateIdx = sstates(istate);
    dataMat = boutLen(:, :, istate)';
    ph = plot([1, 2], dataMat);
    ylabel('Bout Length (s)');
    title([snames{istate}]);
    ylim([0, ceil(max(dataMat(:)))]);
    xlim([0.5 2.5])
    xticks([1, 2])
    xticklabels({'Auto', 'Man'});
end
legend(mnames)

% Plot bout length for each state
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