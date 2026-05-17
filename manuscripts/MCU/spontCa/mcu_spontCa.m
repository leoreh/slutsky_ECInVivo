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



%% ========================================================================
%  FILTERING (SELECTIVE)
%  ========================================================================
% Parameter screen for analysis-time refinement. Mandatory filtering
% (NaN events, cells with zero events in either compartment) already
% happened in spontCa_detectWrapper. Pass 1 metrics (amp/dur/flux/snr/
% fluxOther) are stable under filtering; only Pass 2 (pairing + cell
% aggregates) is refreshed afterwards via mode='cellOnly'.

[tblCell, tblEvent] = spontCa_filter(tblCell, tblEvent, ...
    'minAmp',    [0; 0], ...
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

% Legacy interactive triage (kept as a 1-liner reference):
% tblGUI_bar(tblCell, 'yVar', 'fluxRate', 'xVar', 'compartment', 'grpVar', 'genotype');


% Fraction of unpaired events is roughly the same
tblGrp = tblEvent(tblEvent.genotype == 'Control' & tblEvent.compartment == 'Cyto', :);
sum(isnan(tblGrp.pairIdx)) / height(tblGrp)
tblGrp = tblEvent(tblEvent.genotype == 'MCU-KO' & tblEvent.compartment == 'Cyto', :);
sum(isnan(tblGrp.pairIdx)) / height(tblGrp)


tblCmp = tblEvent(tblEvent.compartment == 'Mito', :);
unpairIdx = find(isnan(tblCmp.pairIdx))
tblCmp(unpairIdx, :)

tblCmp = tblEvent(tblEvent.compartment == 'Cyto', :);
iei = diff(tblCmp.start)
min(abs(iei))

tblCmp = tblEvent(tblEvent.compartment == 'Cyto', :);
tblCmp = tblCmp.genotype == 'MCU-KO';
sum(tblCmp.pairLag < 0.67 & tblCmp.pairLag > -0.38) / height(tblCmp) * 100


%% ========================================================================
%  STATS
%  ========================================================================
% Four predefined tables: T1 cell-level genotype contrasts, T2 event-level
% genotype contrasts, T3 cell-level relationships (main + interaction),
% T4 event-level relationships.

% stats = spontCa_stats(tblCell, tblEvent);


frml = 'amp ~ compartment * genotype + (1 | sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblCell, frml, ...
    'flgPlot', false, 'verbose', true, 'flgStnd', false);

% Pair-coupling LMEs at the event level. Two complementary views, NOT
% reciprocal regressions:
%   Mito subset: per mito event, trigger-cyto flux ~ mito flux.
%       pairFlux on a mito row is the max-amp cyto trigger's flux (amp*dt).
%   Cyto subset: per cyto event, paired-mito flux ~ cyto flux.
%       pairFlux on a cyto row is the paired mito's event-level flux (trapz).
% Different samples (cyto N > mito N due to cyto bursts claiming the
% same mito), different aggregation rules, different units of `flux`
% across compartments. Slopes are not expected to be 1/each-other.
% Using paired metrics means there's no need to pre-filter triggered
% events - rows without a partner already have NaN pairFlux and drop out.
tblLme = tblEvent;
tblCmp = tblLme(tblLme.compartment == 'Mito', :);
frml = 'pairFlux ~ flux * genotype + (1 | sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblCmp, frml, ...
    'flgPlot', false, 'verbose', true, 'flgStnd', false);

tblLme = tblCell;
tblCmp = tblLme(tblLme.compartment == 'Cyto', :);
frml = 'load ~ genotype + (1 | sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblCmp, frml, ...
    'flgPlot', false, 'verbose', true, 'flgStnd', false);


tblCmp = tblLme(tblLme.compartment == 'Cyto', :);
frml = 'pairFlux ~ flux * genotype + (1 | sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblCmp, frml, ...
    'flgPlot', false, 'verbose', true, 'flgStnd', false);


hFig = figure;
hAx = nexttile; pdRes = lme_lsmeans(lmeMdl, {'flux', 'genotype'}, 'transParams', lmeInfo.transParams, ...
    'hAx', hAx);
set(gca, 'YScale', 'log')


tblCmp = tblEvent(tblEvent.compartment == 'Cyto', :);
frml = 'amp ~ genotype * paired + (1 | sbjID)'
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblCmp, frml, ...
    'flgPlot', false, 'verbose', true, 'flgStnd', false);

frml = 'amp ~ genotype + (1 | sbjID)'
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblCmp, frml, ...
    'flgPlot', false, 'verbose', true, 'flgStnd', false);

% fluxOther variant: same scaling against the detection-free partner-
% trace integral over the pairing window (cyto row: mito trace in
% [s - winPair(1), s + winPair(2)]; mito row: cyto trace in
% [m - winPair(2), m + winPair(1)], causally flipped). Independent of
% partner-event detection; sanity check on the pairFlux-based slopes
% above.
tblCmp = tblLme(tblLme.compartment == 'Mito', :);
frml = 'fluxOther ~ flux * genotype + (1 | sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblCmp, frml, ...
    'flgPlot', false, 'verbose', true, 'flgStnd', false);

tblCmp = tblLme(tblLme.compartment == 'Cyto', :);
frml = 'fluxOther ~ flux * genotype + (1 | sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblCmp, frml, ...
    'flgPlot', false, 'verbose', true, 'flgStnd', false);





