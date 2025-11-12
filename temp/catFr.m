
%% correct fr

% basepaths = mcu_sessions('wt_bac_off');
% v = basepaths2vars('basepaths', basepaths, 'vars', vars);
% 
% for iPath = 1 : length(basepaths)
% 
%     % get mfr per bout
%     fr = calc_fr(v(iPath).spikes.times, 'basepath', basepaths{iPath},...
%         'graphics', false, 'binsize', 60, 'saveVar', 'fr', 'forceA', true,...
%         'smet', 'none', 'winBL', [0, Inf], 'winCalc', [0, Inf]);
% 
% end


%%

% concatenate var from different sessions
basepaths = mcu_sessions('lh140');


  
close all
hFig = figure;

tiledlayout('flow', 'Padding', 'compact', 'TileSpacing', 'compact');
hAx = gca;

uNames = {'pPYR', 'pINT'};
clr = [0, 0, 1; 1, 0, 0];
lgdTxt = cell(1, 2);


for iUnit = 1 : 2
    
    [frMat, info] = sessions_catVarTime('basepaths', basepaths,...
    'dataPreset', 'fr', 'graphics', false, 'xTicksBinsize', 12,...
    'markRecTrans', true, 'saveFig', false, 'dataAlt', iUnit);
    tstamps = info.x_data / 60 / 60;

    frMat = mea_frDenoise(frMat', info.x_data, 'flgPlot', false);

    plot_stdShade('dataMat', frMat, 'xVal', tstamps,...
        'hAx', hAx, 'clr', clr(iUnit, :));

    nUnits = sum(~all(isnan(frMat), 2));
    lgdTxt{iUnit} = sprintf('%s (n=%d)', uNames{iUnit}, nUnits);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PERTURBATION DETECTION FROM MFR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detect perturbation onset time. Find the N largest MFR drops (most
% negative derivative) and select the one that occurs at the lowest MFR,
% which is characteristic of a major network state transition.

% Crude outlier detection, only for cleaning the MFR before perturbation
% detection
uOtl = isoutlier(range(frMat, 2), 'median', 'ThresholdFactor', 7);
mfr = mean(frMat(~uOtl, :), 1, 'omitnan');
dmfr = diff(mfr);

% Limit the search to a pre-defined window
binsize = 60;
pertWin = round([24, 48] * 60 * 60 / binsize);
pertWin = pertWin(1) : pertWin(2);
% dmfr is one element shorter than mfr, so limit pertWin to valid dmfr indices
pertWin = pertWin(pertWin <= length(dmfr));
dmfrWin = dmfr(pertWin);
mfrWin = mfr(pertWin);

% Find the top N most negative derivatives
nCandidates = 10;
[sortedDrops, sortIdx] = sort(dmfrWin, 'ascend');
nTop = min(nCandidates, length(sortedDrops));
topIdx = sortIdx(1:nTop);

% From these candidates, find the one with the lowest MFR value
candidateMFRs = mfrWin(topIdx);
[~, minMfrIdx] = min(candidateMFRs);

% The perturbation onset is the index of this best candidate
bestCandidateIdx = topIdx(minMfrIdx);
idxPert = pertWin(bestCandidateIdx);
pertTime = idxPert / 60;

plot([pertTime, pertTime], ylim, '--k')