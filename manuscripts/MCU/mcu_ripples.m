






%% ========================================================================
%  ANALYZE (staged: detect -> curate -> analyze)
%  ========================================================================
% The ripple pipeline runs in three separable stages so each session can be
% manually curated between detection and the heavy spike/phase/map analysis.
% Run the loops in order; loop 2 is manual, one mouse at a time.

basepaths = [mcu_basepaths('wt_bsl_ripp'), mcu_basepaths('mcu_bsl'), ...
    mcu_basepaths('ra')];
nFiles = numel(basepaths);
met = ripp_methods('default');          % detection + default QA filter (met.qa)

% Loop 1 - DETECT 
for iFile = 1 : nFiles
    ripp_wrapper('basepath', basepaths{iFile}, 'met', met, 'win', [0 Inf], ...
        'flgSave', true, 'flgForce', true, 'flgDetectOnly', true, ...
        'rippCh', []);
end

% Loop 2 - CURATE + INSPECT 
iFile = 1;
ripp_curate(basepaths{iFile}, 'met', met); % bulk curation GUI

% first open (slow) 
[~, vm, gm] = guiPath(basepaths{iFile}, 'preset', 'ripp');

% reopen (fast)
vm.ripp.data = []; % the ONLY entry re-read
guiPath(basepaths{iFile}, 'varMap', vm, 'guiMap', gm);


% Loop 3 - ANALYZE 
for iFile = 2 : nFiles
    ripp_wrapper('basepath', basepaths{iFile}, 'met', met, 'win', [0 Inf], ...
        'flgSave', true, 'flgForce', true, 'flgDetectOnly', true, ...
        'rippCh', []);
    
    ripp_curate(basepaths{iFile}, 'met', met, 'flgGui', false); % bulk curation GUI

    ripp_analyze(basepaths{iFile}, 'flgPlot', false);
end



%% ========================================================================
%  RATE & DENSITY (STATE-DEPENDENT)
%  ========================================================================

basepaths = [mcu_basepaths('wt_bsl_ripp'), mcu_basepaths('mcu_bsl'), ...
    mcu_basepaths('ra')];
nFiles = length(basepaths);

% RIPPLE STATES
presets = {'rippStates'};
tblStates = mcu_tblVivo('basepaths', basepaths, 'presets', presets);

% NREM Only
tblPlot = tblStates(tblStates.State == 'NREM', :);
tblPlot = tblStates;

guiTbl_bar(tblPlot, 'xVar', 'genotype', 'yVar', 'Density');
guiTbl_scatHist(tblPlot, 'xVar', 'Density', 'yVar', 'Rate', 'grpVar', 'genotype');

% Run LME
frml = 'Density ~ (Duration + Rate) * genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblPlot, frml);



%% ========================================================================
%  RIPP SPIKES
%  ========================================================================

presets = {'rippSpks', 'burst'};
[tbl, ~, ~, xVec] = mcu_tblVivo('basepaths', basepaths, 'presets', presets);

% Select
tblPlot = tbl;
tblPlot = tbl(tbl.unitType == 'RS', :);
% tblPlot(tblPlot.sbjID == 'lh137', :) = [];
% tblPlot.sbjID = removecats(tblPlot.sbjID, {'lh137'});

% Add logit pBurst
tblTrans = tbl_trans(tblPlot, 'varsInc', {'pBurst'}, 'logBase', 'logit');
tblPlot.pBurst_trans = tblTrans.pBurst;
tblTrans = tbl_trans(tblPlot, 'varsInc', {'bRoy'}, 'logBase', 10);
tblPlot.bRoy_trans = tblTrans.bRoy;

% Plot
guiTbl_bar(tblPlot, 'xVar', 'genotype', 'yVar', 'frZ');
guiTbl_scatHist(tblPlot, 'xVar', 'asym', 'yVar', 'bRoy', 'grpVar', 'genotype');
guiTbl_xy(xVec, tbl, 'grpVar', 'genotype');

tblPlot.burstClu = tblPlot.pBurst > 0.25;
tblPlot.pethNorm = normalize(tblPlot.peth, 2, "norm");
tblPlot.pethCumSum = normalize(cumsum(tblPlot.peth, 2), 2, "range");
tblPlot.pethCumSum = cumsum(tblPlot.peth, 2) ./ sum(tblPlot.peth, 2);

guiTbl_xy(xVec, tblPlot, 'grpVar', 'genotype', 'yVar', 'pethCumSum');
xlim([-0.05, 0.05])

% LME
xVar = 'pBurst';
frml = sprintf('com ~ (fr + %s) + genotype + (1|sbjID)', xVar);
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblPlot, frml, 'dist', 'normal');

% Partial Dependence
hFig = figure;
hAx = nexttile;
vars = {xVar, 'genotype'};
[pdRes, hFig] = lme_lsmeans(lmeMdl, vars, 'transParams', lmeInfo.transParams, ...
    'hAx', hAx, 'xLims', {[0, 1], []});

hAx = nexttile;
xVar = 'fr';
vars = {xVar, 'genotype'};
[pdRes, hFig] = lme_lsmeans(lmeMdl, vars, 'transParams', lmeInfo.transParams, ...
    'hAx', hAx);
set(hAx, 'XScale', 'log')

% To Prism
grpIdx = pdRes.genotype == "MCU-KO";
[pdRes(grpIdx, {xVar}), ...
    pdRes(grpIdx, {'com_pred', 'com_upper', 'com_lower'})]

ylim([-2.5, 0.5])
set(gca,'XScale','log')


% Summary
tblSum = groupsummary(tblPlot, {'genotype', 'sbjID'}, 'mean', ...
    vartype("numeric"));

% To Prism (Metrics)
prismMat = tbl2prism(tblPlot, 'yVar', 'com', 'grpVar', 'genotype');
mean(prismMat, 1, 'omitnan');

% To prism (Time)
yVar = 'pethNorm';
grpIdx = tblPlot.genotype == 'MCU-KO';
prismIdx = grpIdx;
nUnits = sum(prismIdx);
prismMat = [mean(tblPlot{prismIdx, yVar}, 1, 'omitnan')', ...
    std(tblPlot{prismIdx, yVar}, [], 1, 'omitnan')', ...
    repmat(nUnits, length(xVec), 1)];


%% ========================================================================
%  RIPPLE PARAMS
%  ========================================================================

presets = {'ripp'};
tblRipp = mcu_tblVivo('basepaths', basepaths, 'presets', presets);
tblPlot = tblRipp(tblRipp.state == 'NREM', :);

% Plot
guiTbl_bar(tblPlot, 'xVar', 'genotype', 'yVar', 'dur');
guiTbl_scatHist(tblRipp, 'xVar', 'dur', 'yVar', 'amp', 'grpVar', 'genotype');

% Summary
% tblSum = groupsummary(tblRipp, {'genotype', 'sbjID'}, 'mean', ...
%     vartype("numeric"))
% tblBsl = tblVivo(tblVivo.day == 'BSL', :);
% tblSum = groupsummary(tblBsl, {'genotype', 'sbjID'}, 'mean', ...
%     vartype("numeric"))

% LME
frml = 'dur ~ (freq + amp + com) * genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblRipp, frml);

frml = 'dur ~ genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblRipp, frml);

% To Prism (Metrics)
prismMat = tbl2prism(tblRipp, 'yVar', 'freq', 'grpVar', 'genotype');
mean(prismMat, 1, 'omitnan');


%% ========================================================================
%  RIPPLE MAPS
%  ========================================================================

presets = {'rippMaps'};
[tblMaps, ~, ~, xVec] = mcu_tblVivo('basepaths', basepaths, 'presets', presets);

% Plot
guiTbl_xy(xVec, tblMaps, 'yVar', 't_lfp', 'grpVar', 'genotype');

% To prism
yVar = 't_freq';
grpIdx = tblMaps.genotype == 'Control';
nRipp = height(tblMaps(grpIdx, :));
prismMat = [mean(tblMaps{grpIdx, yVar}, 1, 'omitnan')', ...
    std(tblMaps{grpIdx, yVar}, [], 1, 'omitnan')', ...
    repmat(nRipp, length(xVec), 1)];








%% ========================================================================
%  POLAR PLOT
%  ========================================================================

% Figure Parameters
hFig = figure;
fntSize = 16; FntName = 'Arial';
txtUnit = cfg.lbl.unit;
txtGrp = cfg.lbl.grp;

% Plot each group
nGrp = length(grps);
iUnit = 1;
for iGrp = 1 : nGrp
    % Get specific data from table
    idxUnit = tblLme.UnitType == categorical(txtUnit(iUnit));
    idxGrp = tblLme.genotype == categorical(txtGrp(iGrp));
    idxSgn = tblLme.pVal < 0.05;
    idxTbl = idxUnit & idxGrp & idxSgn;
    grpTbl = tblLme(idxTbl, :);

    % Plot units, colored by type if population info is available.
    hPlt = polarscatter(grpTbl.Theta, grpTbl.MRL, 50, ...
        cfg.clr.grp(iGrp, :), 'filled', ...
        'MarkerFaceAlpha', 0.3);
    hold on;
end
rlim([0 0.6])
rticks(0 : 0.3 : 1)
thetaticks(0:90:270)
hAx = gca;
hAx.ThetaAxisUnits = 'degrees';
hAx.GridAlpha = 0.2;
legend(txtGrp, 'Location', 'northwest', 'Interpreter', 'none');
set(hAx, 'FontName', 'Arial', 'FontSize', fntSize);
set(hFig, 'Color', 'w');

% Assert Size
plot_axSize('hFig', hFig, 'szOnly', true, 'axShape', 'square', 'axHeight', 300);

% Save
fname = ['Ripp~SpkPolar_', txtUnit{iUnit}];
lme_save('hFig', hFig, 'fname', fname, 'frmt', {'svg', 'mat'});




%% ========================================================================
%  RATE-PHASE MAP
%  ========================================================================

% Select
flgCbar = false;
iGrp = 1;
iUnit = 1;

% get map data
nSgn = nan(2, 2);
prctSgn = nan(2, 2);
nMice = length(v{iGrp});
mapData = cell(nMice, 1);
for iMouse = 1 : nMice
    ripp = v{iGrp}(iMouse).ripp;
    mapData{iMouse} = ripp.spkLfp.rateMap.rate;
end
mapData = cell2padmat(mapData, 3);
rateMap = ripp.spkLfp.rateMap;

% get unit indices
idxGrp = tblLme.genotype == categorical(txtGrp(iGrp));
grpTbl = tblLme(idxGrp, :);
idxUnit = grpTbl.UnitType == categorical(txtUnit(iUnit));
idxSgn = grpTbl.pVal < 0.05;
idxMap = idxUnit & idxSgn;

% store number of significant units
nSgn(iGrp, iUnit) = sum(idxUnit & idxSgn);
prctSgn(iGrp, iUnit) = sum(idxUnit & idxSgn) / sum(idxUnit) * 100;

% Plot Mean Power-Phase Rate Map (averaged across cells).
% This 2D heatmap shows the average firing rate of neurons as a function of
% LFP phase (x-axis) and LFP power (y-axis). The phase axis is duplicated
% (0 to 4*pi) to visualize cyclic nature. A cosine wave is overlaid as a phase reference.


[hFig, hAx] = plot_axSize('szOnly', false);

mapAvg = mean(mapData(:, :, idxMap), 3, 'omitnan'); % Average rate map across units.
imagesc(hAx, rateMap.phaseBins, rateMap.powBins, mapAvg);
hold on
% Overlay cosine wave for phase reference.
plot(hAx, linspace(0, 2*pi, 100), ...
    cos(linspace(0, 2*pi, 100)) * (range(rateMap.powBins)/4) + mean(rateMap.powBins), ...
    'k--', 'LineWidth', 0.5);
axis xy
if flgCbar
    hCb = colorbar;
    hCb.Label.String = 'Firing Rate (Hz)';
end
colormap(hAx, "pink")
clim([0 12])
ylim(hAx, [min(rateMap.powBins), max(rateMap.powBins)])
xlim([0 2 * pi])
xticks(0:pi/2:2*pi)
hAx.XTickLabel = {'0', '90', '180', '270', '360'};
xlabel('Phase (°)')
ylabel('LFP Power (z-score)');
title([cfg.lbl.grp{iGrp}]);
hTtl = get(hAx, 'Title');
set(hTtl, 'FontSize', 18, 'FontWeight', 'bold');

% Assert Size
plot_axSize('hFig', hFig, 'szOnly', false, 'axWidth', 232, 'axHeight', 300);

% Save
fname = ['Ripp~SpkPhaseMap_', cfg.lbl.grp{iGrp}, '_', cfg.lbl.unit{iUnit}];
lme_save('hFig', hFig, 'fname', fname, 'frmt', {'svg', 'mat'});









