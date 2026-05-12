%% mcu_spontCa.m  Spontaneous Ca2+ imaging pipeline (cyto + mito).
%
% PURPOSE
%   Read NF's SpontCa.xlsx into a long-format unit-level table, run
%   detection per row directly, optionally curate events via the manCur
%   GUI, then finalize (curation overlay + aggregates + ETA maps +
%   coupling). Reproduces Fig 1E,F + S1B-E and the transfer-function
%   preview.
%
% PIPELINE FILES (manuscripts/MCU/spontCa)
%   spontCa_load     Excel -> long-format table; returns fs separately
%   spontCa_detect   single-trace event detection (dF/F input)
%   spontCa_manCur   interactive per-cell event-curation GUI
%   spontCa_finalize curation overlay + per-row aggregates + ETA + coupling
%   spontCa_gui      per-cell QC viewer (post-curation)
%
% See also: MCU_TBLMEA, MEA_WRAPPER, TBLGUI_BAR, LME_ANALYSE


%% ========================================================================
%  LOAD
%  ========================================================================

[tbl, fs] = spontCa_load();


%% ========================================================================
%  DETECT (per-compartment params, tune here)
%  ========================================================================
% Detection runs inline per row so spontCa_detect is called directly with
% no wrapper hiding it; in the debugger the row's trace and chosen params
% are visible in scope. Outputs per-row event columns (start/stop/amp/
% dur/int) used as a pre-fill for spontCa_manCur and aggregated by
% spontCa_finalize.
%
% Detection is derivative-based with a local-baseline amplitude gate.
% Each event is a positive derivative crossing whose peak rises above
% the rolling 20th-percentile baseline by at least minAmp.
%   minAmp  - peak amplitude ABOVE LOCAL BASELINE (dF/F).
%   minIEI  - peak-to-peak distance for greedy max-suppression (s).
%   kNoise  - rise-threshold multiplier on per-cell derivative noise.
%   minDur  - minimum decay length, stop - peak (s).
paramsCyto = {'minAmp', 0.05, 'minIEI', 1.0, 'kNoise', 3.5, 'minDur', 0.4};
paramsMito = {'minAmp', 0.03, 'minIEI', 1.0, 'kNoise', 3.5, 'minDur', 0.4};

tbl = spontCa_detectAll(tbl, fs, ...
    'paramsCyto', paramsCyto, 'paramsMito', paramsMito);


%% ========================================================================
%  MANUAL CURATION (interactive)
%  ========================================================================
% Per-cell event-editing GUI. Opens with the most recent saved version
% for each cell (man/ if present, else auto/). Saves to spontCa/man/.

spontCa_manCur(tbl, fs);


%% ========================================================================
%  FINALIZE (overlay curation + aggregates + ETA + coupling)
%  ========================================================================

tbl = spontCa_finalize(tbl, fs, 'thrLag', 3);


%% ========================================================================
%  QC (per-cell viewer)
%  ========================================================================

spontCa_gui(tbl, fs);


%% ========================================================================
%  FIG 1E,F + S1B-E
%  ========================================================================
% Per-compartment LME + bar plot over genotype. With one observation per
% sbjID the random intercept is degenerate and the LME reduces to LM.
% meanAmp is computed inline.

metrics  = {'rate', 'meanAmp', 'flux', 'fluxInt'};
statsAll = struct();

for c = {'Cyto', 'Mito'}
    cmp = c{1};
    sub = tbl(tbl.compartment == cmp, ...
        {'genotype', 'sbjID', 'unitID', 'nEvents', 'rate', ...
         'amp', 'flux', 'fluxInt'});
    sub.meanAmp = cellfun(@(a) mean(a(~isnan(a))), sub.amp);
    sub.amp     = [];

    for m = metrics
        varRsp = m{1};
        frml = sprintf('%s ~ genotype + (1|sbjID)', varRsp);
        try
            [~, lmeStats, ~] = lme_analyse(sub, frml, ...
                'flgPlot', false, 'verbose', false);
            statsAll.(cmp).(varRsp) = lmeStats;
        catch ME
            warning('lme_analyse failed for %s/%s: %s', ...
                cmp, varRsp, ME.message);
        end
        hF = tblGUI_bar(sub, 'yVar', varRsp, 'xVar', 'genotype');
        set(hF, 'Name', sprintf('%s %s', cmp, varRsp));
    end
end


%% ========================================================================
%  VALIDATION : fraction of cyto-independent mito events
%  ========================================================================
% Per-cell fracIndep across genotypes tests the "mito is predominantly
% cyto-driven" model.

subM = tbl(tbl.compartment == 'Mito', ...
    {'genotype', 'sbjID', 'fracIndep'});
tblGUI_bar(subM, 'yVar', 'fracIndep', 'xVar', 'genotype');


%% ========================================================================
%  PREVIEW : transfer function T = flux_mito / flux_cyto per cell
%  ========================================================================
% Pivot to wide on (sbjID, compartment) so the ratio is one row per cell.
% See Atoms/MCU/MCU compensation model.md.

iC = find(tbl.compartment == 'Cyto');
iM = find(tbl.compartment == 'Mito');
assert(isequal(tbl.sbjID(iC), tbl.sbjID(iM)), ...
    'Cyto/Mito rows must be aligned per cell');

wide = table(tbl.sbjID(iC),     tbl.genotype(iC), ...
             tbl.flux(iC),      tbl.flux(iM), ...
             tbl.fluxInt(iC),   tbl.fluxInt(iM), ...
    'VariableNames', {'sbjID', 'genotype', ...
                      'flux_cyto', 'flux_mito', ...
                      'fluxInt_cyto', 'fluxInt_mito'});
wide.T_flux    = wide.flux_mito    ./ wide.flux_cyto;
wide.T_fluxInt = wide.fluxInt_mito ./ wide.fluxInt_cyto;

tblGUI_scatHist(wide, 'xVar', 'flux_cyto',    'yVar', 'T_flux',    'grpVar', 'genotype');
tblGUI_scatHist(wide, 'xVar', 'fluxInt_cyto', 'yVar', 'T_fluxInt', 'grpVar', 'genotype');
