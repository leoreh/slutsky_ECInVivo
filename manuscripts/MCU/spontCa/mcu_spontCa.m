%% mcu_spontCa.m  Spontaneous Ca2+ imaging pipeline (cyto + mito).
%
% PURPOSE
%   Read NF's SpontCa.xlsx into a long-format unit-level table, run
%   detection independently per compartment, couple each mito event to a
%   preceding cyto event post-hoc, reproduce Fig 1E,F + S1B-E, and stage
%   the transfer-function preview.
%
% PIPELINE FILES (manuscripts/MCU/spontCa)
%   spontCa_load     Excel -> long-format table; returns fs separately
%   spontCa_detect   single-trace event detection (dF/F input)
%   spontCa_events   per-row independent detection + ETA maps
%   spontCa_couple   mito-to-preceding-cyto coupling (cytoIndependent flag)
%   spontCa_gui      per-cell QC viewer
%
% See also: MCU_TBLMEA, MEA_WRAPPER, TBLGUI_BAR, LME_ANALYSE


%% ========================================================================
%  LOAD
%  ========================================================================

[tbl, fs] = spontCa_load();


%% ========================================================================
%  DETECT (per-compartment params, tune independently)
%  ========================================================================
% Per-compartment detection: spontCa_detect operates on a single trace, so
% to debug a specific compartment pull a single row out of tbl and call it
% directly. Mito starts identical to cyto so the diff is visible from a
% single set of changes.

paramsCyto = {'kThr', 3, 'minAmp', 0.02, 'minDur', 0.4, 'minIEI', 0.4};
paramsMito = {'kThr', 3, 'minAmp', 0.02, 'minDur', 0.4, 'minIEI', 0.4};

tbl = spontCa_events(tbl, fs, ...
    'paramsCyto', paramsCyto, 'paramsMito', paramsMito);


%% ========================================================================
%  COUPLE (mito -> preceding cyto)
%  ========================================================================

tbl = spontCa_couple(tbl, 'thrLag', 3);


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
% Independent detection lets mito events stand on their own. Per cell,
% spontCa_couple flags each mito event as cytoIndependent if no cyto event
% precedes it within thrLag. fracIndep is the per-cell fraction; comparing
% it across genotypes tests the "mito is predominantly cyto-driven" model.

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
