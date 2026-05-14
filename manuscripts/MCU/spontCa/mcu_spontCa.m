%% mcu_spontCa.m  Post-detection analyses for spontaneous Ca imaging.
%
% PURPOSE
%   Load pre-built tables produced by spontCa_detectWrapper, derive
%   analysis metrics via spontCa2_metrics, apply selective filtering,
%   and run summary plots + focused analyses + interactive exploration.
%   This script does NOT touch the source Excel file - all detection,
%   curation, and assembly happen upstream in spontCa_detectWrapper.
%
% PIPELINE
%   spontCa_detectWrapper  loadXls -> detect -> manCur -> assemble ->
%                          mandatory filter -> save cache/spontCa_tbl.mat.
%   spontCa2_metrics       enrich tblEvent and tblCell with amp/dur/flux/
%                          pairing/aggregates. Re-callable on filtered
%                          subsets.
%   spontCa_explore        interactive two-panel viewer.
%   spontCa_summary        canonical 3x4 distribution + scatter figures
%                          and the LME/OLS battery (genotype contrasts,
%                          scaling slopes).
%
% SECTIONS
%   1. LOAD
%   2. METRICS
%   3. FILTERING (selective)
%   4. SUMMARY PLOTS
%   5. FOCUSED ANALYSES (triggering, compensation, compensation model)
%   6. INTERACTIVE EXPLORATION


%% ========================================================================
%  LOAD
%  ========================================================================
% tblCell carries per-session traces + genotype + compartment. tblEvent
% holds the minimal per-event schema {sbjID, genotype, compartment,
% start, stop} after the mandatory-filter pass in detectWrapper.

spDir = fileparts(which('spontCa_detect'));
load(fullfile(spDir, 'cache', 'spontCa_tbl.mat'), 'tblCell', 'tblEvent', 'fs');


%% ========================================================================
%  METRICS
%  ========================================================================
% Recompute amp / dur / flux from traces, run cyto-mito pairing, attach
% cell-level aggregates. aggFcn drives mean vs median collapsing of
% per-event quantities to per-cell.

[tblCell, tblEvent] = spontCa2_metrics(tblCell, tblEvent, fs, ...
    'aggFcn', 'mean');


%% ========================================================================
%  FILTERING (SELECTIVE)
%  ========================================================================
% Parameter screen. Each criterion is a clearly named scalar; flip them
% to rerun analyses on a different subset. Mandatory filtering (NaN
% events, cells with zero events in either compartment after assembly)
% already happened in spontCa_detectWrapper - this section is purely
% optional analysis-time refinement.

minAmpCyto   = 0.025;   % drop cyto events below this dF/F (detector noise floor)
minNEvents   = 0;       % drop cells with fewer than this many events per compartment

% --- Build event-level keep mask ---
keepEv = true(height(tblEvent), 1);
isCytoEv = tblEvent.compartment == 'Cyto';
keepEv(isCytoEv & tblEvent.amp < minAmpCyto) = false;
fprintf('FILTER: dropped %d cyto events with amp < %g\n', ...
    sum(~keepEv), minAmpCyto);
tblEvent = tblEvent(keepEv, :);

% --- Build cell-level keep mask (cells must retain >= minNEvents per compartment) ---
nCper = arrayfun(@(s) sum(tblEvent.sbjID == s & ...
    tblEvent.compartment == 'Cyto'), tblCell.sbjID);
nMper = arrayfun(@(s) sum(tblEvent.sbjID == s & ...
    tblEvent.compartment == 'Mito'), tblCell.sbjID);
keepCell = nCper > minNEvents & nMper > minNEvents;
fprintf('FILTER: dropped %d cells below minNEvents=%d\n', ...
    sum(~keepCell) / 2, minNEvents);
tblCell  = tblCell(keepCell, :);
tblEvent = tblEvent(ismember(tblEvent.sbjID, tblCell.sbjID), :);

% --- Re-run metrics on the filtered tables so aggregates are consistent ---
[tblCell, tblEvent] = spontCa2_metrics(tblCell, tblEvent, fs, ...
    'aggFcn', 'mean');

% NOTE on minAmpCyto threshold (decision aid). Triggering cytos give the
% empirical "real cyto" amp distribution; below its lower tail is
% detector noise. To revisit, plot histograms of all-cyto vs triggering
% cyto amps in Control cells and pick from the lower tail:
%   mask    = tblEvent.compartment == 'Cyto' & tblEvent.genotype == 'Control';
%   ampAll  = tblEvent.amp(mask);
%   ampTrig = tblEvent.amp(mask & ~isnan(tblEvent.pairIdx));
%   thrSugg = prctile(ampTrig, 5);


%% ========================================================================
%  SUMMARY PLOTS
%  ========================================================================
% Canonical 3x4 figure (event distributions + scatters) and the LME /
% OLS battery. spontCa_summary auto-switches modes on table identity.

[hFigEv, sEv] = spontCa_summary(tblEvent);
[hFigCl, sCl] = spontCa_summary(tblCell);

% Legacy interactive triage (kept as a 1-liner reference):
% tblGUI_bar(tblCell, 'yVar', 'fluxRate', 'xVar', 'compartment', 'grpVar', 'genotype');


%% ========================================================================
%  FOCUSED ANALYSES
%  ========================================================================
% Three questions raised after the canonical summary.
%
% (1) TRIGGERING CYTOS. Are cyto events that triggered a mito event
%     larger than non-triggering cyto events, and does that triggering-
%     related shift differ by genotype?
% (2) COMPENSATION (per-cell, log-log). Do cells with weaker per-event
%     mito uptake compensate with larger cyto events? Tested across
%     three aggregation metrics (amp, flux, fluxRate).
% (3) COMPENSATION MODEL. Do Ctrl and KO cells trace a single continuous
%     curve T = load_mito / load_cyto vs load_cyto? Tested as joint
%     genotype null in a log-log OLS fit (per MCU compensation model).

cfg = mcu_cfg();
clr = cfg.clr.grp;


% --- (1) TRIGGERING CYTOS ----------------------------------------------
% tblEvent.triggered is set by spontCa2_metrics. LME on cyto events:
% amp ~ triggered * genotype + (1|sbjID), log-normal.

cytoEv = tblEvent(tblEvent.compartment == 'Cyto', :);
[mdlTrig, stTrig, infoTrig] = lme_analyse(cytoEv, ...
    'amp ~ triggered * genotype + (1|sbjID)', ...
    'dist', 'Log-Normal', 'flgPlot', false, 'verbose', true, 'flgStnd', false);

hFigTrig = figure('Name', 'Triggering cytos', 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.08 0.20 0.62 0.55]);
tlT = tiledlayout(hFigTrig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
title(tlT, 'Cyto amp split by triggering status', 'FontWeight', 'bold');

axT1 = nexttile(tlT, 1);
sub = cytoEv(cytoEv.genotype == 'Control', :);
plot_hist(sub, 'amp', 'g', 'triggered', 'hAx', axT1, ...
    'c', [0.45 0.45 0.45; 0.05 0.05 0.05], ...
    'flgKDE', true, 'flgStat', true, 'scale', 'log');
set(axT1, 'XScale', 'log');
xlabel(axT1, 'cyto amp (dF/F)'); ylabel(axT1, 'pdf');
title(axT1, sprintf('Control (n=%d)', height(sub)));
legend(axT1, 'Location', 'best', 'Box', 'off');

axT2 = nexttile(tlT, 2);
sub = cytoEv(cytoEv.genotype == 'MCU-KO', :);
plot_hist(sub, 'amp', 'g', 'triggered', 'hAx', axT2, ...
    'c', [0.90 0.75 0.55; 0.55 0.40 0.20], ...
    'flgKDE', true, 'flgStat', true, 'scale', 'log');
set(axT2, 'XScale', 'log');
xlabel(axT2, 'cyto amp (dF/F)'); ylabel(axT2, 'pdf');
title(axT2, sprintf('MCU-KO (n=%d)', height(sub)));
legend(axT2, 'Location', 'best', 'Box', 'off');


% --- (2) COMPENSATION (per-cell log-log) -------------------------------
% Reshape to one row per cell with cyto/mito columns per metric, fit
% log-log OLS with genotype interaction.

isC = tblCell.compartment == 'Cyto';
isM = tblCell.compartment == 'Mito';
tblComp = tblCell(isC, {'sbjID', 'genotype'});
tblComp.cytoAmp      = tblCell.amp(isC);
tblComp.mitoAmp      = tblCell.amp(isM);
tblComp.cytoFlux     = tblCell.flux(isC);
tblComp.mitoFlux     = tblCell.flux(isM);
tblComp.cytoFluxRate = tblCell.fluxRate(isC);
tblComp.mitoFluxRate = tblCell.fluxRate(isM);

metrics   = {'amp', 'flux', 'fluxRate'};
metricLbl = {'mean per-event amp (dF/F)', ...
             'mean per-event flux (dF/F\cdots)', ...
             'fluxRate (dF/F/s)'};
compRows = {};
for iM = 1:numel(metrics)
    m = metrics{iM};
    xCol = ['cyto' upper(m(1)) m(2:end)];
    yCol = ['mito' upper(m(1)) m(2:end)];
    tblFit = tblComp;
    tblFit.(xCol) = log(tblComp.(xCol));
    tblFit.(yCol) = log(tblComp.(yCol));
    mdl   = fitlm(tblFit, sprintf('%s ~ %s * genotype', yCol, xCol));
    coefs = mdl.Coefficients;
    nm    = string(coefs.Properties.RowNames);
    bMain  = coefs.Estimate(nm == xCol);
    pMain  = coefs.pValue(nm == xCol);
    isInt  = contains(nm, xCol) & contains(nm, "genotype");
    bInter = coefs.Estimate(isInt);
    pInter = coefs.pValue(isInt);
    compRows(end+1, :) = {m, bMain, pMain, bInter, pInter, ...
        bMain + bInter, mdl.Rsquared.Ordinary};                  %#ok<*AGROW>
end
tblCompStats = cell2table(compRows, 'VariableNames', ...
    {'metric', 'b_Ctrl', 'p_main', 'b_inter', 'p_inter', 'b_KO', 'R2'});
disp(tblCompStats);

hFigComp = figure('Name', 'Compensation: mito ~ cyto per cell', ...
    'Color', 'w', 'Units', 'normalized', 'Position', [0.06 0.20 0.86 0.55]);
tlC = tiledlayout(hFigComp, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
title(tlC, 'Compensation: mito ~ cyto per cell (log-log)', 'FontWeight', 'bold');
for iM = 1:numel(metrics)
    m = metrics{iM};
    xCol = ['cyto' upper(m(1)) m(2:end)];
    yCol = ['mito' upper(m(1)) m(2:end)];
    ax = nexttile(tlC, iM);
    set(ax, 'XScale', 'log', 'YScale', 'log');
    plot_scat(tblComp, xCol, yCol, 'g', 'genotype', ...
        'hAx', ax, 'c', clr, 'fitType', 'Linear', 'flgStats', true, ...
        'sz', 40, 'alpha', 0.7);
    xlabel(ax, sprintf('cyto %s', metricLbl{iM}));
    ylabel(ax, sprintf('mito %s', metricLbl{iM}));
    rowS = tblCompStats(strcmp(tblCompStats.metric, m), :);
    title(ax, sprintf('%s | Ctrl %+.2f, KO %+.2f, p(int)=%.3f', ...
        m, rowS.b_Ctrl, rowS.b_KO, rowS.p_inter));
end


% --- (3) COMPENSATION MODEL TEST ---------------------------------------
% Single continuous T(load_cyto) curve across genotypes? Compensation
% supported iff adding genotype to a pooled fit of log T ~ log load_cyto
% does not improve it (joint F-test on all genotype-related coefficients).
% Uses the detection-free `load` from spontCa2_metrics.

isC = tblCell.compartment == 'Cyto';
isM = tblCell.compartment == 'Mito';
tblT = table(tblCell.sbjID(isC), tblCell.genotype(isC), ...
    tblCell.load(isC), tblCell.load(isM), ...
    'VariableNames', {'sbjID', 'genotype', 'loadCyto', 'loadMito'});
ok = isfinite(tblT.loadCyto) & tblT.loadCyto > 0 & ...
     isfinite(tblT.loadMito) & tblT.loadMito > 0;
tblT = tblT(ok, :);
tblT.T          = tblT.loadMito ./ tblT.loadCyto;
tblT.logLcyto   = log(tblT.loadCyto);
tblT.logT       = log(tblT.T);
tblT.genotype   = setcats(tblT.genotype, {'Control', 'MCU-KO'});

mdlPool = fitlm(tblT, 'logT ~ logLcyto');
mdlGeno = fitlm(tblT, 'logT ~ logLcyto * genotype');

% Build the H matrix from coefficient names so the joint test does not
% depend on MATLAB's ordering convention (categorical before continuous).
coefNames = mdlGeno.CoefficientNames;
genoIdx   = find(contains(coefNames, 'genotype'));
H = zeros(numel(genoIdx), numel(coefNames));
for k = 1:numel(genoIdx)
    H(k, genoIdx(k)) = 1;
end
[p_geno, F_geno, df1] = coefTest(mdlGeno, H);
df2 = mdlGeno.DFE;

fprintf('\n=== Compensation model test (per-cell T vs load_cyto) ===\n');
fprintf('n = %d cells (%d Ctrl, %d KO).\n', height(tblT), ...
    sum(tblT.genotype == 'Control'), sum(tblT.genotype == 'MCU-KO'));
fprintf('\nPooled fit:  log T ~ log load_cyto\n');
disp(mdlPool.Coefficients);
fprintf('R^2 = %.3f\n', mdlPool.Rsquared.Ordinary);
fprintf('\nWith-genotype fit:  log T ~ log load_cyto * genotype\n');
disp(mdlGeno.Coefficients);
fprintf('R^2 = %.3f\n', mdlGeno.Rsquared.Ordinary);
fprintf('\nJoint test of genotype contribution:\n');
fprintf('  F(%d, %d) = %.2f, p = %.4f\n', df1, df2, F_geno, p_geno);
if p_geno >= 0.05
    fprintf('  -> genotype does NOT improve the fit. Compensation supported.\n');
else
    fprintf('  -> genotype DOES improve the fit. Compensation rejected.\n');
end
dAIC = mdlPool.ModelCriterion.AIC - mdlGeno.ModelCriterion.AIC;
dBIC = mdlPool.ModelCriterion.BIC - mdlGeno.ModelCriterion.BIC;
fprintf('Delta AIC (pooled - genotype) = %+.2f (positive favors genotype)\n', dAIC);
fprintf('Delta BIC (pooled - genotype) = %+.2f (positive favors genotype)\n', dBIC);

hFigT = figure('Name', 'Compensation model: T vs load_cyto', 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.08 0.22 0.72 0.55]);
tlT = tiledlayout(hFigT, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
title(tlT, sprintf('Compensation model: T vs load_{cyto} (joint geno p = %.3f)', p_geno), ...
    'FontWeight', 'bold');
axT1 = nexttile(tlT, 1);
set(axT1, 'XScale', 'log', 'YScale', 'log');
plot_scat(tblT, 'loadCyto', 'T', 'g', 'genotype', ...
    'hAx', axT1, 'c', clr, 'fitType', 'Linear', 'flgStats', true, ...
    'sz', 50, 'alpha', 0.85);
xlabel(axT1, 'load_{cyto} (\DeltaF/F)');
ylabel(axT1, 'T = load_{mito} / load_{cyto}');
title(axT1, sprintf('T vs load_{cyto} (n=%d cells)', height(tblT)));

axT2 = nexttile(tlT, 2);
set(axT2, 'XScale', 'log', 'YScale', 'log');
plot_scat(tblT, 'loadCyto', 'loadMito', 'g', 'genotype', ...
    'hAx', axT2, 'c', clr, 'fitType', 'Linear', 'flgStats', true, ...
    'sz', 50, 'alpha', 0.85);
xlabel(axT2, 'load_{cyto} (\DeltaF/F)');
ylabel(axT2, 'load_{mito} (\DeltaF/F)');
title(axT2, 'load_{mito} vs load_{cyto}');


%% ========================================================================
%  INTERACTIVE EXPLORATION
%  ========================================================================
% Two-panel viewer. Panel 1: within-compartment scatter+hist. Panel 2:
% cross-compartment scatter+hist (one row per paired cyto-mito cell).
% Top-bar controls switch Level (Event/Cell), Compartment (Cyto/Mito),
% and the cross-compartment Metric.

spontCa_explore(tblEvent, tblCell);
