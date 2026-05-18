%% mcu_spontCa.m  Post-detection analyses for spontaneous Ca imaging.
%
% PURPOSE
%   Load pre-built tables produced by spontCa_detectWrapper, derive
%   analysis metrics via spontCa2_metrics, apply selective filtering,
%   run the predefined stats tables, and open the interactive explorer.
%   This script does NOT touch the source Excel file - all detection,
%   curation, and assembly happen upstream in spontCa_detectWrapper.
%
% PIPELINE
%   spontCa_detectWrapper  loadXls -> detect -> manCur -> assemble ->
%                          mandatory filter -> save cache/spontCa_tbl.mat.
%   spontCa2_metrics       enrich tblEvent and tblCell with amp/dur/flux/
%                          pairing/aggregates. Re-callable on filtered
%                          subsets.
%   spontCa_stats          four predefined tables: genotype contrasts and
%                          relationships at cell + event level.
%   spontCa_explore        interactive two-panel viewer.
%
% SECTIONS
%   1. LOAD
%   2. METRICS
%   3. FILTERING (selective)
%   4. STATS
%   5. INTERACTIVE EXPLORATION


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


% out = spontCa_sweep(tblCell, tblEvent, fs, 'var', 'amp', 'flgPair', 'all');           



%% ========================================================================
%  FILTERING (SELECTIVE)
%  ========================================================================
% Parameter screen for analysis-time refinement. Mandatory filtering
% (NaN events, cells with zero events in either compartment) already
% happened in spontCa_detectWrapper. Pass 1 metrics (amp/dur/flux/snr/
% fluxOther) are stable under filtering; only Pass 2 (pairing + cell
% aggregates) is refreshed afterwards via mode='cellOnly'.

[tblCell, tblEvent] = spontCa_filter(tblCell, tblEvent, ...
    'minAmp',    [0; 0.03], ...
    'minSNR',    [0; 0], ...
    'minEvents', [0; 0], ...
    'flgPair',   'all');

[tblCell, tblEvent] = spontCa2_metrics(tblCell, tblEvent, fs, ...
    'mode', 'cellOnly', 'aggFcn', 'mean');



%% ========================================================================
%  INTERACTIVE EXPLORATION
%  ========================================================================
% Two-panel viewer. Panel 1: within-compartment scatter+hist. Panel 2:
% cross-compartment scatter+hist (one row per paired cyto-mito cell).
% Top-bar controls switch Level (Event/Cell), Compartment (Cyto/Mito),
% and the cross-compartment Metric.

spontCa_explore(tblEvent, tblCell, fs);

% tblGUI_bar(tblEvent, 'yVar', 'fluxRate', 'xVar', 'compartment', 'grpVar', 'genotype');



%% ========================================================================
%  STATS
%  ========================================================================
% Four predefined tables: T1 cell-level genotype contrasts, T2 event-level
% genotype contrasts, T3 cell-level relationships (main + interaction),
% T4 event-level relationships.

% stats = spontCa_stats(tblCell, tblEvent);


frml = 'ampRate ~ compartment * genotype + (1 | sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblCell, frml, 'dist', 'log-normal', ...
    'flgPlot', false, 'verbose', true, 'flgStnd', false);


frml = 'amp ~ compartment * genotype + (1 | sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblEvent, frml, ...
    'flgPlot', false, 'verbose', true, 'flgStnd', false);

v = 'amp';
g = categories(tblEvent.genotype);
c = ["Mito" "Cyto"];
M = nan(2,6);
for i = 1:2
    for j = 1:2
        x = tblEvent.(v)(tblEvent.compartment==c(i) & tblEvent.genotype==g{j});
        M(i,(j-1)*3+(1:3)) = [mean(x,'omitnan'), std(x,'omitnan'), numel(x)];
    end
end


hFig = figure;

tblLme = tblEvent;
idx = tblLme.compartment == 'Cyto' & tblLme.paired;
% idx = tblLme.compartment == 'Cyto';
tblLme = tblLme(idx, :);

var = 'flux';
frml = ['fluxOther ~ ', var, ' * genotype + (1 | sbjID)'];
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLme, frml, 'dist', 'log-normal', ...
    'flgPlot', false, 'verbose', true, 'flgStnd', false);

hAx = nexttile; 
pdRes = lme_lsmeans(lmeMdl, {var, 'genotype'}, 'transParams', lmeInfo.transParams, ...
    'hAx', hAx);





tblLme = tblCell;
idx = tblLme.compartment == 'Cyto';
tblLme = tblLme(idx, :);

frml = 'pairAmp ~ amp * genotype + (1 | sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLme, frml, ...
    'flgPlot', false, 'verbose', true, 'flgStnd', false);

hAx = nexttile; 
pdRes = lme_lsmeans(lmeMdl, {'amp', 'genotype'}, 'transParams', lmeInfo.transParams, ...
    'hAx', hAx);



var = 'amp';
tblLme = tblCell(:, {'genotype', 'sbjID', 'compartment', var});
tblWide = unstack(tblLme, var, "compartment");
frml = 'Mito ~ Cyto * genotype';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblWide, frml, 'dist', 'log-normal', ...
    'flgPlot', false, 'verbose', true, 'flgStnd', false);


tblLme = tblCell;
frml = 'amp ~ compartment * genotype + (1 | sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLme, frml, ...
    'flgPlot', false, 'verbose', true, 'flgStnd', false);

hAx = nexttile;
pdRes = lme_lsmeans(lmeMdl, {'amp', 'genotype'}, 'transParams', lmeInfo.transParams, ...
    'hAx', hAx);


%% ========================================================================
%  BIOPHYSICS vs COMPENSATION  (within-between decomposition)
%  ========================================================================
% A single regression of pairAmp on amp mixes two distinct relationships
% that live on different scales:
%
%   WITHIN-cell  (per event)   MCU transfers Ca from cyto to mito. Bigger
%                              cyto event -> bigger mito partner. Slope
%                              positive. Blunted in KO. Biophysics.
%
%   BETWEEN-cell (cell mean)   Cells with effective MCU pull Ca into
%                              mito; cells with poor MCU let cyto pile
%                              up while mito stays small. Cross-cell
%                              slope negative. Altered in KO. Compensation.
%
% An event-level LME with (1|sbjID) reports the within slope only.
% A cell-level LME on tblCell reports the between slope only, with the
% aggregation issues below. The Mundlak decomposition fits BOTH in one
% event-level LME by splitting log(amp) into:
%
%       la_b = cell mean of log(amp)         between-cell predictor
%       la_w = log(amp) - la_b               within-cell deviation
%
% la_w and la_b are orthogonal by construction; genotype interactions
% test MCU dependence at each scale independently. This is the principled
% disentanglement and the headline result of this section.
%
% Variable choice: amp / pairAmp here. To repeat for flux, swap amp ->
% flux and pairAmp -> pairFlux throughout.


% ---- Build event-level table with within/between predictors -----------
% Paired cyto events only. amp and pairAmp are different events on
% different traces (mutual NN by time, not amplitude), so this regression
% has no algebraic ratio confound. log scale because both are skewed.

tblWB       = tblEvent(tblEvent.compartment == 'Cyto' & tblEvent.paired, :);
tblWB.la    = log(tblWB.amp);
tblWB.lp    = log(tblWB.pairAmp);
[G, ~]      = findgroups(tblWB.sbjID);
cellMeanLA  = splitapply(@mean, tblWB.la, G);
tblWB.la_b  = cellMeanLA(G);
tblWB.la_w  = tblWB.la - tblWB.la_b;


% ---- Fit Mundlak model -------------------------------------------------
% Predictors and response are pre-logged. dist='normal' skips response
% transform; flgStnd=false keeps slopes in log-log units (interpretable
% as power-law exponents).

frmlWB = 'lp ~ la_w * genotype + la_b * genotype + (1 | sbjID)';
[mdlWB, statsWB, infoWB] = lme_analyse(tblWB, frmlWB, ...
    'dist', 'normal', 'flgPlot', false, 'verbose', true, ...
    'flgStnd', false);


% ---- Print scale-by-scale slopes with genotype contrasts ---------------

c          = mdlWB.Coefficients;
nm         = string(c.Name);
sl_w_ctrl  = c.Estimate(nm == "la_w");
sl_b_ctrl  = c.Estimate(nm == "la_b");
sl_w_d_i   = find(contains(nm, "la_w") & contains(nm, "genotype"), 1);
sl_b_d_i   = find(contains(nm, "la_b") & contains(nm, "genotype"), 1);
sl_w_d     = c.Estimate(sl_w_d_i);  p_w_d = c.pValue(sl_w_d_i);
sl_b_d     = c.Estimate(sl_b_d_i);  p_b_d = c.pValue(sl_b_d_i);

fprintf('\n=== Within-between decomposition (log pairAmp ~ log amp) ===\n');
fprintf('  WITHIN-cell  (biophysical transfer)\n');
fprintf('     Ctrl slope   %+0.3f\n', sl_w_ctrl);
fprintf('     KO   slope   %+0.3f   (delta %+0.3f, p = %0.3g)\n', ...
    sl_w_ctrl + sl_w_d, sl_w_d, p_w_d);
fprintf('  BETWEEN-cell (compensation)\n');
fprintf('     Ctrl slope   %+0.3f\n', sl_b_ctrl);
fprintf('     KO   slope   %+0.3f   (delta %+0.3f, p = %0.3g)\n', ...
    sl_b_ctrl + sl_b_d, sl_b_d, p_b_d);
fprintf('\nReading the deltas:\n');
fprintf('  WITHIN delta  <0  -> KO has flatter per-event transfer\n');
fprintf('  BETWEEN delta <0  -> KO has stronger negative cross-cell slope\n\n');
disp(c);


% ---- Visualize the two slopes side by side -----------------------------
% Both panels stay on log-log scale; the fitted relationships are linear
% lines there. Slopes match the printed coefficients above.
%   Left:  vary la_w with la_b fixed at its mean (within-cell view).
%   Right: vary la_b with la_w fixed at 0 = cell's own mean.

hFig = figure('Color', 'w', 'Position', [200 200 1000 420]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

hAx = nexttile;
lme_lsmeans(mdlWB, {'la_w', 'genotype'}, ...
    'transParams', infoWB.transParams, 'hAx', hAx);
xlabel('log(amp) deviation from cell mean'); ylabel('log(pairAmp)');
title('Within-cell  (biophysical transfer)');

hAx = nexttile;
lme_lsmeans(mdlWB, {'la_b', 'genotype'}, ...
    'transParams', infoWB.transParams, 'hAx', hAx);
xlabel('log(amp) cell mean'); ylabel('log(pairAmp)');
title('Between-cell  (compensation)');


%% ========================================================================
%  ROBUSTNESS: pair-fraction confound + aggregation alignment
%  ========================================================================
% The within-cell slope is a clean per-event regression and needs no
% extra check. The between-cell slope is more fragile and has one
% specific aggregation risk worth confirming:
%
%   tblCell.amp     = mean over ALL cyto events.
%   tblCell.pairAmp = mean over PAIRED cyto events' mito partners.
%
% If cells with high mean amp pair a smaller fraction of their cyto
% events, then amp and pairAmp are summarized from drifting subsets and
% the cell-level slope is part biology, part bookkeeping. Two diagnostics:


% ---- Diagnostic 1: pair fraction vs cell-mean amp ----------------------
% Pearson r between pairing fraction and log cell-mean amp across cyto
% cells. Rule of thumb: |r| > 0.3 means aggregation contributes
% materially and you should rely on Diagnostic 2.

cellsC   = tblCell(tblCell.compartment == 'Cyto', :);
nCyto    = arrayfun(@(s) sum(tblEvent.sbjID == s & ...
    tblEvent.compartment == 'Cyto'), cellsC.sbjID);
nPaired  = arrayfun(@(s) sum(tblEvent.sbjID == s & ...
    tblEvent.compartment == 'Cyto' & tblEvent.paired), cellsC.sbjID);
pairFrac = nPaired ./ max(nCyto, 1);

[rPF, pPF] = corr(log(cellsC.amp), pairFrac, 'rows', 'complete');
fprintf('\n=== Pair-fraction confound check ===\n');
fprintf('  r(pairFrac, log cellMean amp) = %+0.3f   p = %0.3g\n', rPF, pPF);
if abs(rPF) > 0.3
    fprintf('  -> aggregation contributes. Rely on diagnostic 2 below.\n');
else
    fprintf('  -> aggregation negligible. Cell-level slope is trustworthy.\n');
end


% ---- Diagnostic 2: refit cell-level with aligned aggregation -----------
% Replace amp (mean over all cyto events) with ampPaired (mean over
% paired cyto events only). Now amp and pairAmp are summarized on the
% SAME event subset. If the negative slope and KO contrast survive,
% the cell-level compensation reading is robust to aggregation.

ampPaired = nan(height(cellsC), 1);
for k = 1:height(cellsC)
    sub = tblEvent(tblEvent.sbjID == cellsC.sbjID(k) & ...
                   tblEvent.compartment == 'Cyto' & tblEvent.paired, :);
    if ~isempty(sub), ampPaired(k) = mean(sub.amp); end
end
tblAlt           = cellsC;
tblAlt.ampPaired = ampPaired;
tblAlt           = tblAlt(~isnan(tblAlt.ampPaired) & ~isnan(tblAlt.pairAmp), :);

[mdlOrig, ~, infoOrig] = lme_analyse(cellsC, ...
    'pairAmp ~ amp * genotype + (1 | sbjID)', ...
    'flgPlot', false, 'verbose', true, 'flgStnd', false);

[mdlAlt, ~, infoAlt] = lme_analyse(tblAlt, ...
    'pairAmp ~ ampPaired * genotype + (1 | sbjID)', ...
    'flgPlot', false, 'verbose', true, 'flgStnd', false);

hFig = figure('Color', 'w', 'Position', [250 250 1000 420]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

hAx = nexttile;
lme_lsmeans(mdlOrig, {'amp', 'genotype'}, ...
    'transParams', infoOrig.transParams, 'hAx', hAx);
title('Cell level   amp = all cyto events   (original)');

hAx = nexttile;
lme_lsmeans(mdlAlt, {'ampPaired', 'genotype'}, ...
    'transParams', infoAlt.transParams, 'hAx', hAx);
title('Cell level   amp = paired cyto only   (aligned)');
