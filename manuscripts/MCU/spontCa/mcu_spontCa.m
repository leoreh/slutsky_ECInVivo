%% mcu_spontCa.m  Spontaneous Ca2+ imaging pipeline (cyto + mito).
%
% PURPOSE
%   Read NF's SpontCa.xlsx into a long-format cell table and a long-format
%   events table, optionally curate via manCur, finalize (aggregates +
%   ETA maps + coupling). Reproduces Fig 1E,F + S1B-E and the
%   transfer-function preview.
%
% PIPELINE FILES (manuscripts/MCU/spontCa)
%   spontCa_load        Excel -> tblCell (traces + cell metadata)
%   spontCa_detect      single-trace event detection, returns a table
%   spontCa_writeEvents tblEvent -> per-cell <sbjID>.mat files
%   spontCa_readEvents  per-cell <sbjID>.mat files -> tblEvent
%   spontCa_manCur      interactive per-cell event-curation GUI
%   spontCa_finalize    aggregates on tblCell + coupling on tblEvent
%   spontCa_gui         per-cell QC viewer (post-finalize)
%
% See also: MCU_TBLMEA, MEA_WRAPPER, TBLGUI_BAR, LME_ANALYSE


%% ========================================================================
%  LOAD
%  ========================================================================

[tblCell, fs] = spontCa_load();


%% ========================================================================
%  DETECT
%  ========================================================================
% Per-row detection. Each call to spontCa_detect returns a table with
% rows = events; sbjID + compartment tags are added before vertcat into
% tblEvent. Result is written to <spontCa>/auto/<sbjID>.mat (one bare
% events table per cell) for the manCur Load button to pick up.
%
% Derivative-based detection with a local-baseline amplitude gate:
%   minAmp - peak amplitude above local baseline (dF/F)
%   minIEI - peak-to-peak distance for greedy max-suppression (s)
%   kNoise - rise-threshold multiplier on per-cell derivative noise
%   minDur - minimum decay length, stop - peak (s)

paramsCyto = {'minAmp', 0.05, 'minIEI', 1.0, 'kNoise', 3.5, 'minDur', 0.4};
paramsMito = {'minAmp', 0.03, 'minIEI', 1.0, 'kNoise', 3.5, 'minDur', 0.4};

n = height(tblCell);
chunks = cell(n, 1);
for iRow = 1:n
    if tblCell.compartment(iRow) == 'Cyto'
        rowEv = spontCa_detect(tblCell.trace(iRow, :), fs, paramsCyto{:});
    else
        rowEv = spontCa_detect(tblCell.trace(iRow, :), fs, paramsMito{:});
    end
    if height(rowEv) > 0
        rowEv.sbjID       = repmat(tblCell.sbjID(iRow),       height(rowEv), 1);
        rowEv.compartment = repmat(tblCell.compartment(iRow), height(rowEv), 1);
        chunks{iRow} = rowEv;
    end
end
tblEvent = vertcat(chunks{~cellfun(@isempty, chunks)});
tblEvent = tblEvent(:, ['sbjID', 'compartment', setdiff(...
    tblEvent.Properties.VariableNames, {'sbjID','compartment'}, 'stable')]);

autoDir = fullfile(fileparts(which('spontCa_detect')), 'auto');
spontCa_writeEvents(tblEvent, autoDir, 'backup', false);


%% ========================================================================
%  MANUAL CURATION (interactive)
%  ========================================================================
% Opens manCur on the in-memory tblEvent. Save writes per-cell bare
% events tables to man/<sbjID>.mat. After closing, re-read whatever's
% on disk in man/ (if any) and merge with the auto-detection rows for
% cells the user didn't curate.

spontCa_manCur(tblCell, tblEvent, fs);

manDir = fullfile(fileparts(which('spontCa_detect')), 'man');
tblEvent_man = spontCa_readEvents(manDir);
if height(tblEvent_man) > 0
    curatedCells = unique(tblEvent_man.sbjID);
    tblEvent = tblEvent(~ismember(tblEvent.sbjID, curatedCells), :);
    tblEvent = [tblEvent; tblEvent_man];
end


%% ========================================================================
%  FINALIZE (aggregates + ETA + coupling)
%  ========================================================================

[tblCell, tblEvent] = spontCa_finalize(tblCell, tblEvent, fs, 'thrLag', 3);


%% ========================================================================
%  QC (per-cell viewer)
%  ========================================================================

spontCa_gui(tblCell, tblEvent, fs);


%% ========================================================================
%  FIG 1E,F + S1B-E
%  ========================================================================
% Per-compartment LME + bar plot over genotype. With one observation per
% sbjID the random intercept is degenerate and the LME reduces to LM.

metrics  = {'rate', 'meanAmp', 'flux', 'fluxInt'};
statsAll = struct();

for c = {'Cyto', 'Mito'}
    cmp = c{1};
    sub = tblCell(tblCell.compartment == cmp, ...
        {'genotype', 'sbjID', 'unitID', 'nEvents', 'rate', ...
         'meanAmp', 'flux', 'fluxInt'});

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

subM = tblCell(tblCell.compartment == 'Mito', ...
    {'genotype', 'sbjID', 'fracIndep'});
tblGUI_bar(subM, 'yVar', 'fracIndep', 'xVar', 'genotype');


%% ========================================================================
%  PREVIEW : transfer function T = flux_mito / flux_cyto per cell
%  ========================================================================
% Pivot to wide on (sbjID, compartment) so the ratio is one row per cell.

iC = find(tblCell.compartment == 'Cyto');
iM = find(tblCell.compartment == 'Mito');
assert(isequal(tblCell.sbjID(iC), tblCell.sbjID(iM)), ...
    'Cyto/Mito rows must be aligned per cell');

wide = table(tblCell.sbjID(iC),     tblCell.genotype(iC), ...
             tblCell.flux(iC),      tblCell.flux(iM), ...
             tblCell.fluxInt(iC),   tblCell.fluxInt(iM), ...
    'VariableNames', {'sbjID', 'genotype', ...
                      'flux_cyto', 'flux_mito', ...
                      'fluxInt_cyto', 'fluxInt_mito'});
wide.T_flux    = wide.flux_mito    ./ wide.flux_cyto;
wide.T_fluxInt = wide.fluxInt_mito ./ wide.fluxInt_cyto;

tblGUI_scatHist(wide, 'xVar', 'flux_cyto',    'yVar', 'T_flux',    'grpVar', 'genotype');
tblGUI_scatHist(wide, 'xVar', 'fluxInt_cyto', 'yVar', 'T_fluxInt', 'grpVar', 'genotype');
