

%% ========================================================================
%  RE-ANALYSE SPECIFIC STEPS
%  ========================================================================

basepaths = [mcu_basepaths('wt_bsl_ripp'), mcu_basepaths('mcu_bsl')];
nFiles = length(basepaths);

% Load Session & Data
vars = {'session', 'spikes', 'brst', 'ripp', 'units'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

clear rippSpks
for iFile = 1 : nFiles

    basepath = basepaths{iFile};
    cd(basepath)
    [~, basename] = fileparts(basepath);

    ripp = v(iFile).ripp;
    uType = v(iFile).units.type;

    % Prepare Single-Unit (SU) Spike Times
    spkTimes = v(iFile).spikes.times;
    nUnits = length(spkTimes);

    % Firing Metrics
    rippSpks = ripp_spks(spkTimes, ...
        ripp.times, ...
        ripp.ctrlTimes, ...
        ripp.peakTime, ...
        'unitType', uType, ...
        'flgSave', false);

    % Add per ripple spike metrics to ripp struct
    ripp.spks = rippSpks.all.events;
    rippSpks = rmfield(rippSpks.all, 'events');

    % Save
    rippSpksFile = fullfile(basepath, [basename, '.rippSpks.mat']);
    rippFile = fullfile(basepath, [basename, '.ripp.mat']);
    save(rippSpksFile, 'rippSpks', '-v7.3');
    save(rippFile, 'ripp', '-v7.3');

end




%
% funcon = [];
% for iFile = 1 : nFiles
%     basepath = basepaths{iFile};
%     cd(basepath)
%     [~, basename] = fileparts(basepath);
%
%     load([basename, '.rippSpks.mat'])
%     cc = fr_corr(rippSpks.su.rippRates, 'nShuffles', 50, 'flgPlot', false);
%     funcon = [funcon; cc.shuffle.funcon];
% end




%% ========================================================================
%  ANALYZE
%  ========================================================================

basepaths = [mcu_basepaths('wt_bsl_ripp'), mcu_basepaths('mcu_bsl')];
nFiles = length(basepaths);

for iFile = 1 : nFiles
    basepath = basepaths{iFile};
    cd(basepath)
    [~, basename] = fileparts(basepath);
    tic
    ripp = ripp_wrapper('basepath', pwd, ...
        'win', [0 12] * 3600, ...
        'rippCh', [], ...
        'flgPlot', false, ...
        'flgSave', true, ...
        'steps', 'spks');
    toc
end




%% ========================================================================
%  RATE & DENSITY (STATE-DEPENDENT)
%  ========================================================================

basepaths = [mcu_basepaths('wt_bsl_ripp'), mcu_basepaths('mcu_bsl')];
nFiles = length(basepaths);

% RIPPLE STATES
presets = {'rippStates'};
tblStates = mcu_tblVivo('basepaths', basepaths, 'presets', presets);

% NREM Only
tblPlot = tblStates(tblStates.State == 'NREM', :);
tblPlot = tblStates;

tblGUI_bar(tblPlot, 'xVar', 'Group', 'yVar', 'Density');
tblGUI_scatHist(tblPlot, 'xVar', 'Density', 'yVar', 'Rate', 'grpVar', 'Group');

% Run LME
frml = 'Density ~ (Duration + Rate) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblPlot, frml);



%% ========================================================================
%  RIPP SPIKES
%  ========================================================================

presets = {'rippSpks'};
[tbl, ~, ~, xVec] = mcu_tblVivo('basepaths', basepaths, 'presets', presets);

% Select
tblPlot = tbl;
tblPlot = tbl(tbl.UnitType == 'RS', :);
% tblPlot(tblPlot.Name == 'lh137', :) = [];
% tblPlot.Name = removecats(tblPlot.Name, {'lh137'});

% Plot
tblGUI_bar(tblPlot, 'xVar', 'Group', 'yVar', 'frZ');
tblGUI_scatHist(tblPlot, 'xVar', 'asym', 'yVar', 'bRoy', 'grpVar', 'Group');
tblGUI_xy(xVec, tbl);

% LME
frml = 'com ~ Group * bRoy + (1|Name)';
[lmeMdl, lmeStates, lmeInfo] = lme_analyse(tblPlot, frml);

% Summary
tblSum = groupsummary(tblPlot, {'Group', 'Name'}, 'mean', ...
    vartype("numeric"));




%% ========================================================================
%  RIPPLE PARAMS
%  ========================================================================

presets = {'ripp'};
tblRipp = mcu_tblVivo('basepaths', basepaths, 'presets', presets);

% Plot
tblGUI_bar(tblRipp, 'xVar', 'Group', 'yVar', 'dur');
tblGUI_scatHist(tblRipp, 'xVar', 'dur', 'yVar', 'amp', 'grpVar', 'Group');

% Summary
tblSum = groupsummary(tblRipp, {'Group', 'Name'}, 'mean', ...
    vartype("numeric"));

% Prism
tblSum = groupsummary(tblMap(:, {'Group', 't_lfp'}), {'Group'}, {'mean', 'std'}, ...
    vartype("numeric"));
tblSum.mean_t_lfp';
repmat(tblSum.GroupCount(2), 127, 1)

% LME
frml = 'dur ~ (freq + amp) * Group + (1|Name)';
[lmeMdl, lmeStates, lmeInfo] = lme_analyse(tblRipp, frml);




%% ========================================================================
%  RIPPLE MAPS
%  ========================================================================

presets = {'rippMaps'};
[tblMaps, ~, ~, xVec] = mcu_tblVivo('basepaths', basepaths, 'presets', presets);

% Plot
tblGUI_xy(xVec, tblMaps, 'yVar', 't_z', 'grpVar', 'states');







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
    idxGrp = tblLme.Group == categorical(txtGrp(iGrp));
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
idxGrp = tblLme.Group == categorical(txtGrp(iGrp));
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









