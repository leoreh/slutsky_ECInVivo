%% mcu_spontCa.m  Spontaneous Ca2+ imaging pipeline (cyto + mito).
%
% PURPOSE
%   Read NF's SpontCa.xlsx into a long-format unit-level table, detect
%   events with a cyto-triggered scheme (one mito event per cyto event;
%   see Atoms/MCU/MCU compensation model.md), reproduce Fig 1E,F + S1B-E,
%   and stage the transfer-function preview.
%
% PIPELINE FILES (manuscripts/MCU/spontCa)
%   spontCa_load     Excel -> long-format table; returns fs separately
%   spontCa_detect   single-trace event detection (dF/F input)
%   spontCa_events   table-level orchestrator (cytoTrigger | independent)
%   spontCa_gui      per-cell QC viewer
%
% See also: MCU_TBLMEA, MEA_WRAPPER, TBLGUI_BAR, LME_ANALYSE


%% ========================================================================
%  LOAD & DETECT (primary: cyto-triggered mito)
%  ========================================================================

[tbl, fs] = spontCa_load();
tbl = spontCa_events(tbl, fs);


%% ========================================================================
%  QC (per-cell viewer)
%  ========================================================================

spontCa_gui(tbl, fs);


%% ========================================================================
%  FIG 1E,F + S1B-E
%  ========================================================================
% Per-compartment LME + bar plot over genotype. With one observation per
% sbjID the random intercept is degenerate and the LME reduces to LM.
% meanAmp / meanDur are computed inline.

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
%  VALIDATION : independent mito detection
%  ========================================================================
% Mito events detected without using cyto as trigger. Orphan mito events
% (those without a nearby cyto event) test the modelling assumption that
% mito is predominantly cyto-driven.

tblIndep = spontCa_events(tbl, fs, 'mode', 'independent');

ids = unique(tbl.sbjID);
cmp = table();
cmp.sbjID         = ids;
cmp.nCyto         = nan(length(ids), 1);
cmp.nMito_trig    = nan(length(ids), 1);
cmp.nMito_indep   = nan(length(ids), 1);
for i = 1:length(ids)
    cmp.nCyto(i)       = tbl.nEvents(tbl.sbjID == ids(i) & tbl.compartment == 'Cyto');
    cmp.nMito_trig(i)  = tbl.nEvents(tbl.sbjID == ids(i) & tbl.compartment == 'Mito');
    cmp.nMito_indep(i) = tblIndep.nEvents(tblIndep.sbjID == ids(i) & tblIndep.compartment == 'Mito');
end
fprintf('\nMito events per cell (median):\n');
fprintf('  cyto-triggered : %d\n', median(cmp.nMito_trig));
fprintf('  independent    : %d\n', median(cmp.nMito_indep));
fprintf('  orphans (indep - trig, may be negative if cyto fires without mito): %d\n', ...
    median(cmp.nMito_indep - cmp.nMito_trig));


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
