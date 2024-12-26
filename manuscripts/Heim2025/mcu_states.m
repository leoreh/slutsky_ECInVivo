
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% manually curate state epochs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

queryStr = 'wt_bsl';
basepaths = [mcu_sessions(queryStr)];
nfiles = length(basepaths);

% select session
ifile = 3;
basepath = basepaths{ifile};
[~, basename] = fileparts(basepath);
cd(basepath)

% create bkup of manual labels
labelsrevfile = [basename, '.sleep_labelsRev.mat'];
if exist(labelsrevfile, "file")
    load(labelsrevfile)
    mkdir(fullfile(basepath, 'ss'))
    ctime = datetime('now', 'Format', 'yyMMdd_HHmmss');
    bkupName = [basename, '.labelsRev', '_', char(ctime), '.mat'];
    bkupFile = fullfile(basepath, 'ss', bkupName);
    save(bkupFile, 'labels')
else
    load([basename, '.sleep_states.mat'])
    labels = ss.labels;
end

% load sleep sig 
sSig = load([basename, '.sleep_sig.mat']);

% manually revise network labels
AccuSleep_viewer(sSig, labels, labelsrevfile)

% create bkup of labelsMan
labelsmanfile = [basename, '.sleep_labelsMan.mat'];
ctime = datetime('now', 'Format', 'yyMMdd_HHmmss');
bkupName = [basename, '.labelsMan', '_', char(ctime), '.mat'];
bkupFile = fullfile(basepath, 'ss', bkupName);
copyfile(labelsmanfile, bkupFile)

% replace labelsMan with revised labels
idx = 1 : 6 * 60 * 60;
load(labelsrevfile)
labelsRev = labels;
labels = ones(length(labelsRev), 1) * 8;
labels(idx) = labelsRev(idx);
AccuSleep_viewer(sSig, labels, labelsmanfile)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reclassify states automatically
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % classify with a network
% netfile = 'D:\Code\slutsky_ECInVivo\lfp\SleepStates\AccuSleep\trainedNetworks\net_230212_103132.mat';
% calData = [];
% % calData = ss.info.calibrationData;
% ss = as_classify(sSig, 'basepath', basepath, 'inspectLabels', false,...
%     'saveVar', true, 'forceA', true, 'netfile', netfile,...
%     'graphics', true, 'calData', calData);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% states by emg across mice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mice = [mcu_sessions('wt')];

for imouse = 1 : length(mice)
    ssEmg = mcu_StatesEmgMouse(mice{imouse}, true);
end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % compare manual and automatic classification
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % -------------------------------------------------------------------------
% % organize and load data
% 
% % files with manual curation
% basepaths = {...
%     'E:\Data\lh132\lh132_230413_094013',...
%     'E:\Data\lh133\lh133_230413_094013',...
%     'E:\Data\lh96\lh96_220120_090157',...
%     'E:\Data\lh142\lh142_231005_091832',...
%     'E:\Data\lh107\lh107_220518_091200',...
%     };
% nfiles = length(basepaths);
% 
% % load data
% % basepaths = xls2basepaths('mname', 'lh132', 'pcond', "tempflag");
% v = basepaths2vars('basepaths', basepaths, 'vars',...
%     ["sleep_labelsMan"; "sleep_states"]);
% 
% % -------------------------------------------------------------------------
% % find indices of manual curation
% 
% % file
% ifile = 5; 
% basepath = basepaths{ifile};
% cd(basepath)
% [~, basename] = fileparts(basepath);
% 
% % load automatic labels
% labelsAuto = v(ifile).ss.labels_net;
% 
% % load sleep signal and inspect manLabels
% % sSig = load([basename, '.sleep_sig.mat']);
% % AccuSleep_viewer(sSig, v(ifile).labels, [])
% % AccuSleep_viewer(sSig, labelsAuto, [])
% 
% % -------------------------------------------------------------------------
% % get cirated indices from labelMan.m file
% 
% labelsMan = v(ifile).labels;
% 
% % Find indices where labels are not 8
% curated = labelsMan ~= 8;
% 
% % Find transitions (start and end of scored segments)
% diffs = diff([0; curated; 0]);  % Add padding to detect edges
% startIdx = find(diffs == 1);
% endIdx = find(diffs == -1) - 1;
% 
% % Filter out segments that are entirely label 8
% validSegments = arrayfun(@(s, e) any(labelsMan(s:e) ~= 8), startIdx, endIdx);
% 
% % Extract start and end indices of valid segments
% startIdx = startIdx(validSegments);
% endIdx = endIdx(validSegments);
% 
% % Display results
% manIdx = [startIdx : endIdx];
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot confusion matrix
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % confusion matrix (relative to manual labels)
% [ss.netPrecision, ss.netRecall] = as_cm(labelsMan, labelsAuto,...
%     'saveFig', false);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % calculate stateEpochs for Man and Auto
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % state params
% sstates = [1, 4, 5];
% nstates = length(sstates);
% 
% % auto
% interDur = 5;
% minDur = 5;
% [stateEpochs, epochStats] = as_epochs('labels', labelsAuto(manIdx),...
%     'minDur', minDur, 'interDur', interDur);
% 
% % manual
% interDur = 5;
% minDur = 5;
% [stateEpochs_man, epochStats_man] = as_epochs('labels', labelsMan(manIdx),...
%     'minDur', minDur, 'interDur', interDur);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plot Comparison
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Colors and state names
% cfg = as_loadConfig();
% clr = cfg.colors(sstates);
% snames = cfg.names(sstates);
% 
% % Open figure
% fh = figure;
% setMatlabGraphics(true);
% set(fh, 'WindowState','maximized');
% tlayout = [3, length(sstates)];
% th = tiledlayout(tlayout(1), tlayout(2));
% th.TileSpacing = 'tight';
% th.Padding = 'none';
% title(th, basename, 'interpreter', 'none');
% 
% % Plot epoch length for each state
% tilebias = 0;
% for istate = 1 : length(sstates)
%     axh = nexttile(th, tilebias + istate); hold on; cla;
%     stateIdx = sstates(istate);
%     epLen = cell2padmat({epochStats.epLen{stateIdx}, epochStats_man.epLen{stateIdx}}, 2);
%     plot_boxMean('axh', axh, 'dataMat', epLen, 'plotType', 'bar', 'clr', clr{(istate)});
%     ylabel('Epoch Length (s)');
%     title([snames{istate}]);
%     xticklabels({'Auto', 'Man'});
% end
% 
% % Plot number of epochs
% tilebias = 3;
% nepochs = [epochStats.nepochs(sstates); epochStats_man.nepochs(sstates)]';
% for istate = 1 : length(sstates)
%     axh = nexttile(th, tilebias + istate); hold on; cla;
%     stateIdx = sstates(istate);
%     plot_boxMean('axh', axh, 'dataMat', nepochs(istate, :), 'plotType', 'bar', 'clr', clr{(istate)});
%     ylabel('Number Epochs');
%     title([snames{istate}]);
%     xticklabels({'Auto', 'Man'});
% end
% 
% % Plot number of epochs
% tilebias = 6;
% prctDur = [epochStats.prctDur(sstates); epochStats_man.prctDur(sstates)]';
% for istate = 1 : length(sstates)
%     axh = nexttile(th, tilebias + istate); hold on; cla;
%     stateIdx = sstates(istate);
%     plot_boxMean('axh', axh, 'dataMat', prctDur(istate, :), 'plotType', 'bar', 'clr', clr{(istate)});
%     ylabel('State Duration (%)');
%     title([snames{istate}]);
%     xticklabels({'Auto', 'Man'});
% end










% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % compare manual and automatic classification (RA data)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% queryStr = 'mcu_bsl';
% vars = ["sleep_states"];
% [basepaths, v] = mcu_sessions(queryStr, vars);
% nfiles = length(basepaths);
% 
% 
% % params
% interDur = 0;
% minDur = 0;
% sstates = [1, 4, 5];
% nstates = length(sstates);
% 
% % file
% ifile = 2;
% basepath = basepaths{ifile};
% cd(basepath)
% [~, basename] = fileparts(basepath);
% ss = v(ifile).ss;
% 
% idx = 1 : 6 * 60 * 60;
% idx = 1 : length(ss.labels);
% [stateEpochs, epochStats] = as_epochs('labels', ss.labels_man(idx),...
%     'minDur', minDur, 'interDur', interDur);
% [stateEpochs2, epochStats2] = as_epochs('labels', ss.labels_net(idx),...
%     'minDur', minDur, 'interDur', interDur);
% 
% epLen = cell(nstates, 1);
% for istate = 1 : nstates
%     tmp{1} = epochStats.epLen{sstates(istate)};
%     tmp{2} = epochStats2.epLen{sstates(istate)};
%     epLen{istate} = cell2nanmat(tmp, 2);
% end
% 
% prctDur = epochStats.prctDur(sstates);
% prctDur(2, :) = epochStats2.prctDur(sstates);
% 
% 
% % -------------------------------------------------------------------------
% % compare state dependent firing rate
% load([basename, '.spikes.cellinfo.mat'])
% load([basename, '.units.mat'])
% 
% fr = calc_fr(spikes.times, 'basepath', basepath,...
%     'graphics', false, 'binsize', 60, 'saveVar', false, 'forceA', true,...
%     'smet', 'none', 'winBL', [0, Inf], 'winCalc', [0, Inf], 'stateEpochs', stateEpochs);
% 
% fr2 = calc_fr(spikes.times, 'basepath', basepath,...
%     'graphics', false, 'binsize', 60, 'saveVar', false, 'forceA', true,...
%     'smet', 'none', 'winBL', [0, Inf], 'winCalc', [0, Inf], 'stateEpochs', stateEpochs2);
% 
% iunit = 1;
% unitIdx = units.clean(iunit, :);
% clear mfr
% for istate = 1 : nstates
%     tmp{1} = fr.states.mfr(unitIdx, sstates(istate));
%     tmp{2} = fr2.states.mfr(unitIdx, sstates(istate));
%     mfr{istate} = cell2nanmat(tmp, 2);
% end
% 
% % calculate % change
% clear prctdiff
% for istate = 1 : nstates
%     prctdiff(:, istate) = (mfr{istate}(:, 1) - mfr{istate}(:, 2)) ./...
%         ((mfr{istate}(:, 1) + mfr{istate}(:, 2)) / 2) * 100;
% end