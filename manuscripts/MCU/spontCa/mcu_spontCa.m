%% mcu_spontCa.m  Spontaneous Ca2+ imaging pipeline (cyto + mito).
%
% PURPOSE
%   Read NF's SpontCa.xlsx into a long-format cell table and a long-format
%   events table, optionally curate via manCur, finalize (aggregates +
%   ETA maps + coupling). Reproduces Fig 1E,F + S1B-E and the
%   transfer-function preview.
%
% PIPELINE FILES (manuscripts/MCU/spontCa)
%   spontCa_loadXls     Excel -> tblCell (traces + cell metadata).
%                       Picks sheet 'f' (raw F) or 'dff' via flgRaw;
%                       optionally drops experimenter-excluded cells.
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

[tblCell, fs] = spontCa_loadXls('flgRaw', true, 'flgExclude', false);


%% ========================================================================
%  DF/F  (rolling 20th-percentile baseline)
%  ========================================================================
% F0(t) = 20th percentile of F over a ~30 s window centered at t.
% dF/F(t) = (F(t) - F0(t)) / F0(t). Deterministic recipe; no per-cell tuning.
% prctile skips NaN by default so excluded-cell rows pass through as NaN.

bslWin = 30;                    % baseline window (s)
bslQuant = 20;                  % percentile for baseline
winSamps = round(bslWin * fs);
halfWin = floor(winSamps / 2);
nSamps = size(tblCell.trace, 2);
nRows = height(tblCell);

for iRow = 1:nRows
    f  = tblCell.trace(iRow, :);
    f0 = zeros(1, nSamps);
    for iSamp = 1:nSamps
        i0 = max(1, iSamp - halfWin);
        i1 = min(nSamps, iSamp + halfWin);
        f0(iSamp) = prctile(f(i0:i1), bslQuant);
    end
    tblCell.trace(iRow, :) = (f - f0) ./ max(f0, eps);
end


%% ========================================================================
%  AUTOMATIC DETECT
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

% Params tuned by spontCa_tune against 4 curated cells (Ctrl_01/02/03/05).
paramsCyto = {'minAmp', 0.05, 'minIEI', 1.0, 'kNoise', 3.5, 'minDur', 0.4};
paramsMito = {'minAmp', 0.06, 'minIEI', 0.4, 'kNoise', 3.5, 'minDur', 0.2};

chunks = cell(nRows, 1);
for iRow = 1:nRows
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
%  MANUAL CURATION
%  ========================================================================
% Opens manCur on the in-memory tblEvent. Save writes per-cell bare
% events tables to man/<sbjID>.mat. After closing, re-read whatever's
% on disk in man/ (if any) and merge with the auto-detection rows for
% cells the user didn't curate.

spontCa_manCur(tblCell, tblEvent, fs);

% COMMENTS:
% Control_72, Control_73, and Control_76 appear exactly the same cell. Kept
% only 72.


%% ========================================================================
%  ORGANIZE
%  ========================================================================

manDir = fullfile(fileparts(which('spontCa_detect')), 'man');
tblEvent_man = spontCa_readEvents(manDir);
if height(tblEvent_man) > 0
    curatedCells = unique(tblEvent_man.sbjID);
    tblEvent = tblEvent(~ismember(tblEvent.sbjID, curatedCells), :);
    canonVars = {'sbjID', 'compartment', 'start', 'stop', 'amp', 'dur', 'int'};
    tblEvent     = tblEvent(:,     canonVars);
    tblEvent_man = tblEvent_man(:, canonVars);
    tblEvent = [tblEvent; tblEvent_man];
end

% Attach genotype to each event row (sbjID -> genotype lookup).
[~, idx] = ismember(tblEvent.sbjID, tblCell.sbjID);
tblEvent.genotype = tblCell.genotype(idx);

% Drop cells with zero events in either compartment.
rmvId  = unique(tblCell.sbjID(tblCell.nEvents == 0));
TblLme = tblCell(~ismember(tblCell.sbjID, rmvId), :);
tblEvent = tblEvent(~ismember(tblEvent.sbjID, rmvId), :);

% Remove cells excluded by Neta
% tblLme = tblCell(tblCell.excluded, :);


%% ========================================================================
%  FIG 1E,F + S1B-E
%  ========================================================================
% Per-compartment LME + bar plot over genotype. With one observation per
% sbjID the random intercept is degenerate and the LME reduces to LM.

tblGUI_bar(tblLme, 'yVar', 'meanAmp', 'xVar', 'compartment', 'grpVar', 'genotype');
tblGUI_bar(tblEvent, 'yVar', 'meanAmp', 'xVar', 'compartment', 'grpVar', 'genotype');

% LME (over cells)
frml = 'meanAmp ~ genotype * compartment + (1 | sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLme, frml, 'flgPlot', false, 'verbose', true);

% LME (over events)
frml = 'amp ~ genotype * compartment + (1 | sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblEvent, frml, 'flgPlot', false, 'verbose', true);



%% ========================================================================
%  FINALIZE (aggregates + ETA + coupling)
%  ========================================================================

[tblCell, tblEvent] = spontCa_finalize(tblCell, tblEvent, fs, 'thrLag', 3);



% Maybe best transfer is the integral of cytoCa for each mitoCa event. 

%% ========================================================================
%  QC (per-cell viewer)
%  ========================================================================

spontCa_gui(tblCell, tblEvent, fs);





%% ========================================================================
%  VALIDATION : fraction of cyto-independent mito events
%  ========================================================================
% Obsolete after manual curation

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
