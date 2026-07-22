% MCU_RIPPLES  SWR figures for the MCU manuscript (three genotypes).
%
% Each figure gets its own section. A section (1) builds the table with an
% mcu_tblVivo preset, (2) opens it in a table GUI for inspection, (3) writes
% a Prism-ready block to the clipboard - raw replicates via tbl2prism, or
% Mean/SD/N via tbl2prismSum - and (4) runs the matching LME, as in
% mcu_lme2xls. Sections run independently.
%
% NREM only. The curation state scope is NREM (ripp_methods qa.states), so
% .accepted is NREM-confined, and the 'ripp' / 'rippMaps' / 'rippStates'
% presets return the accepted subset - every metric below is already NREM.
%
% Cohort. mcu_basepaths('bsl3_ripp') = Control (ripple-curated WT sessions),
% MCU-KO (germline), CAG-MCU-KO (viral). One baseline session per mouse.


%% ========================================================================
%  PIPELINE (detect -> curate -> analyze)
%  ========================================================================
% Three separable stages so each session can be curated by hand between
% detection and the heavy spike/phase analysis. Loop 2 is manual, one mouse
% at a time; loops 1 and 3 are batch.

basepaths = mcu_basepaths('bsl3_ripp');
nFiles = numel(basepaths);
met = ripp_methods('default');          % detection + default QA filter (met.qa)

% Loop 1 - DETECT
for iFile = 1 : nFiles
    ripp_wrapper('basepath', basepaths{iFile}, 'met', met, 'win', [0 Inf], ...
        'flgSave', true, 'flgForce', true, 'flgDetectOnly', true, ...
        'rippCh', []);
end

% Loop 2 - CURATE + INSPECT
iFile = 16;
ripp_curate(basepaths{iFile}, 'met', met);          % bulk curation GUI

[~, vm, gm] = guiPath(basepaths{iFile}, 'preset', 'ripp');     % first open (slow)
vm.ripp.data = [];                                             % the ONLY entry re-read
guiPath(basepaths{iFile}, 'varMap', vm, 'guiMap', gm);         % reopen (fast)

% Loop 3 - ANALYZE
for iFile = 13 : nFiles
    ripp_analyze(basepaths{iFile}, 'flgPlot', false);
end


%% ========================================================================
%  SWR PROPERTIES (freq, amp, dur)
%  ========================================================================
% Per-event frequency, amplitude and duration across all NREM SWRs. Two
% views: every event (bars + Mean/SD/N), and one average per mouse (points).

basepaths = mcu_basepaths('bsl3_ripp');
tblRipp = mcu_tblVivo('basepaths', basepaths, 'presets', {'ripp'});

% All events - GUI + Prism (Mean/SD/N)
guiTbl_bar(tblRipp, 'xVar', 'genotype', 'yVar', 'freq', 'grpVar', 'sbjID');
prismFreq = tbl2prismSum(tblRipp, 'yVar', 'freq');
prismAmp  = tbl2prismSum(tblRipp, 'yVar', 'amp');
prismDur  = tbl2prismSum(tblRipp, 'yVar', 'dur');

tblRipp.amp(tblRipp.genotype == 'CAG-MCU-KO');

% One average per mouse - GUI (points) + Prism (raw)
tblMouse = groupsummary(tblRipp, {'genotype', 'sbjID'}, 'mean', ...
    {'freq', 'amp', 'dur'});
guiTbl_bar(tblMouse, 'xVar', 'genotype', 'yVar', 'mean_freq', 'mode', 'points');
prismFreqMouse = tbl2prism(tblMouse, 'yVar', 'mean_freq', 'grpVar', 'genotype');
prismAmpMouse  = tbl2prism(tblMouse, 'yVar', 'mean_amp',  'grpVar', 'genotype');
prismDurMouse  = tbl2prism(tblMouse, 'yVar', 'mean_dur',  'grpVar', 'genotype');

% LME (events, mouse as random effect) - matches Table S5
frml = 'freq ~ genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblRipp, frml, 'dist', 'normal');

frml = 'amp ~ genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblRipp, frml, 'dist', 'log-normal');

frml = 'dur ~ genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblRipp, frml, 'dist', 'log-normal');


%% ========================================================================
%  SWR RATE (per mouse)
%  ========================================================================
% Rate = number of NREM SWRs / total NREM duration, one value per mouse.
% rippStates is per-bout (Rate = count / boutDuration), so pool bouts:
% total events = sum(Rate .* Duration), total time = sum(Duration).

basepaths = mcu_basepaths('bsl3_ripp');
tblStates = mcu_tblVivo('basepaths', basepaths, 'presets', {'rippStates'});
nrem = tblStates(tblStates.State == 'NREM', :);

nrem.nEvt = nrem.Rate .* nrem.Duration;             % events per bout
tblRate = groupsummary(nrem, {'genotype', 'sbjID'}, 'sum', {'nEvt', 'Duration'});
tblRate.rate = tblRate.sum_nEvt ./ tblRate.sum_Duration;    % Hz
tblRate.genotype = removecats(tblRate.genotype);

% GUI (points) + Prism (raw)
guiTbl_bar(tblRate, 'xVar', 'genotype', 'yVar', 'rate', 'mode', 'points');
prismRate = tbl2prism(tblRate, 'yVar', 'rate', 'grpVar', 'genotype');

% Stat: one value per mouse -> one-way genotype comparison (OLS = ANOVA)
mdlRate = fitlm(tblRate, 'rate ~ genotype');
anova(mdlRate)                                      % omnibus genotype effect
mdlRate.Coefficients                                % pairwise vs Control


%% ========================================================================
%  SWR WAVEFORM & FIRING
%  ========================================================================
% LFP waveform and the normalised peri-SWR firing rate of RS units. Two Prism
% exports: a group sheet averaged across events (Mean/SD/N over time), and one
% XY block per genotype where each mouse carries its own Mean/SD/N - overlay the
% per-mouse traces to show the between-mouse spread (<20 mice, mcu_ed style).

basepaths = mcu_basepaths('bsl3_ripp');

% LFP waveform (per-event maps)
[tblMaps, ~, ~, xMaps] = mcu_tblVivo('basepaths', basepaths, 'presets', {'rippMaps'});
guiTbl_xy(xMaps, tblMaps, 'yVar', 't_lfp', 'grpVar', 'genotype');
prismLfp = tbl2prismSum(tblMaps, 'yVar', 't_lfp', 'xVec', xMaps);
prismLfp(:, 1) = prismLfp(:, 1) * 1000;

% Firing rate, normalised per unit (per-unit peri-SWR PETH)
[tblSpk, ~, ~, xSpk] = mcu_tblVivo('basepaths', basepaths, ...
    'presets', {'rippSpks'}, 'flgClean', true);
tblSpk.pethNorm = normalize(tblSpk.peth, 2, 'norm');
guiTbl_xy(xSpk, tblSpk, 'yVar', 'pethNorm', 'grpVar', 'genotype');
prismFr = tbl2prismSum(tblSpk, 'yVar', 'pethNorm', 'xVec', xSpk);
prismFr(:, 1) = prismFr(:, 1) * 1000;

% Overlaid per mouse, one Prism XY block per genotype (each mouse a Mean/SD/N
% triple across its own events) - <20 mice, so show the spread directly rather
% than average it away, as in mcu_ed. Tile by genotype and colour by mouse to
% inspect the overlay; wv2prism copies the 'copy' level, re-run with the next
% genotype (or paste wvLfp(iGrp).str).
guiTbl_xy(xMaps * 1000, tblMaps, 'yVar', 't_lfp', 'tileVar', 'genotype', ...
    'grpVar', 'sbjID', 'xLbl', 'time (ms)');
wvLfp = wv2prism(tblMaps, xMaps * 1000, 'yVar', 't_lfp', 'grpVar', 'sbjID', ...
    'splitVar', 'genotype', 'xLbl', 'time (ms)', 'copy', 'CAG-MCU-KO');

% FR is already normalised (unitless), so scale 1 - not wv2prism's uV->mV default.
guiTbl_xy(xSpk * 1000, tblSpk, 'yVar', 'pethNorm', 'tileVar', 'genotype', ...
    'grpVar', 'sbjID', 'xLbl', 'time (ms)');
    wvFr = wv2prism(tblSpk, xSpk * 1000, 'yVar', 'pethNorm', 'grpVar', 'sbjID', ...
        'splitVar', 'genotype', 'xLbl', 'time (ms)', 'scale', 1, 'copy', 'CAG-MCU-KO');


%% ========================================================================
%  SWR DISTRIBUTIONS
%  ========================================================================
% Per-event distributions across all NREM SWRs. The scatter GUI shows the
% joint plot with marginal histograms; switch X/Y in the dropdowns to view
% amp, freq (instantaneous Hilbert) and freqPeak (1/f-corrected). Prism gets
% the raw replicates per genotype for its own histograms.

basepaths = mcu_basepaths('bsl3_ripp');
tblRipp = mcu_tblVivo('basepaths', basepaths, 'presets', {'ripp'});

guiTbl_scatHist(tblRipp, 'xVar', 'freq', 'yVar', 'amp', 'grpVar', 'genotype');

prismDistAmp      = tbl2prism(tblRipp, 'yVar', 'amp',      'grpVar', 'genotype');
prismDistFreq     = tbl2prism(tblRipp, 'yVar', 'freq',     'grpVar', 'genotype');
prismDistFreqPeak = tbl2prism(tblRipp, 'yVar', 'freqPeak', 'grpVar', 'genotype');


%% ========================================================================
%  SPIKE ORGANIZATION DURING SWR
%  ========================================================================
% Timing of RS spikes within the SWR (center of mass, com), and whether it
% tracks firing rate and burstiness. Matches Table S4; the LS-means panel is
% the Spike-CoM vs P_burst figure.

basepaths = mcu_basepaths('bsl3_ripp');
[tblRipp, ~, ~, xVec] = mcu_tblVivo('basepaths', basepaths, ...
    'presets', {'rippSpks', 'burst'}, 'flgClean', true);
tblTrans = tbl_trans(tblRipp, 'varsInc', {'pBurst'}, 'logBase', 'logit');
tblRipp.pBurst_trans = tblTrans.pBurst;

% GUI
guiTbl_scatHist(tblRipp, 'xVar', 'pBurst', 'yVar', 'com', 'grpVar', 'genotype');
guiTbl_xy(xVec, tblRipp, 'yVar', 'peth', 'grpVar', 'genotype');

% LME (genotype only, then adjusted for fr + burstiness)
frml = 'com ~ genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblRipp, frml, 'dist', 'normal', 'flgStnd', false);
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);

frml = 'com ~ (fr + pBurst) + genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblRipp, frml, 'dist', 'normal', 'flgStnd', false);
lmeTbls = [lmeTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];

% Partial dependence: predicted com vs P_burst, then vs fr
hFig = figure;
hAx = nexttile;
pdBurst = lme_lsmeans(lmeMdl, {'pBurst', 'genotype'}, ...
    'transParams', lmeInfo.transParams, 'hAx', hAx, 'xLims', {[0, 1], []});
hAx = nexttile;
pdFr = lme_lsmeans(lmeMdl, {'fr', 'genotype'}, ...
    'transParams', lmeInfo.transParams, 'hAx', hAx);
set(hAx, 'XScale', 'log')

% Prism: the predicted curves live in pdBurst / pdFr (genotype, grid,
% com_pred, com_lower, com_upper). Raw com per mouse for the points plot.
tblComMouse = groupsummary(tblRipp, {'genotype', 'sbjID'}, 'mean', 'com');
prismCom = tbl2prism(tblComMouse, 'yVar', 'mean_com', 'grpVar', 'genotype');

prismDistCom      = tbl2prism(tblRipp, 'yVar', 'com',      'grpVar', 'genotype');


%% ========================================================================
%  STATE DEPENDENCE (exploratory)
%  ========================================================================
% Per-bout SWR rate and density by vigilance state. Kept for reference; the
% short-bout bias here is why the rate figure above pools to one value per
% mouse instead.

basepaths = mcu_basepaths('bsl3_ripp');
tblStates = mcu_tblVivo('basepaths', basepaths, 'presets', {'rippStates'});
tblPlot = tblStates(tblStates.State == 'NREM', :);

guiTbl_bar(tblPlot, 'xVar', 'genotype', 'yVar', 'Rate');
guiTbl_scatHist(tblPlot, 'xVar', 'Duration', 'yVar', 'Rate', 'grpVar', 'genotype');

frml = 'Rate ~ genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblPlot, frml);


%% ========================================================================
%  SPIKE-PHASE POLAR (legacy)
%  ========================================================================
% Spike-LFP phase coupling per unit, coloured by genotype. Legacy view -
% expects a preloaded tblLme (from ripp_screen) with Theta/MRL/pVal and cfg.

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
%  RATE-PHASE MAP (legacy)
%  ========================================================================
% Mean firing rate as a function of LFP phase and power, averaged across
% significant units. Legacy view - expects preloaded v{iGrp} and tblLme.

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
