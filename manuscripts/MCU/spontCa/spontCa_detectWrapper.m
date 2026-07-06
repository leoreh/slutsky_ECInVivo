%% spontCa_detectWrapper.m  Detection-pipeline wrapper.
%
% PURPOSE
%   Take NetaF's SpontCa.xlsx through:
%       loadXls -> (autodetect) -> (manual curation) -> assemble ->
%       mandatory filter -> save.
%   Output is one cached .mat at cache/spontCa_tbl.mat holding
%   tblCell, tblEvent, fs. mcu_spontCa loads from there without ever
%   touching the .xlsx.
%
%   Run as a script: sections are designed to be executed
%   independently. AUTO DETECT and MANUAL CURATION are commented out by
%   default - they are re-run only when changing detection params or
%   adding manual curation.
%
% PIPELINE FILES
%   spontCa_loadXls     Excel -> tblCell (traces + cell metadata).
%   spontCa_detect      single-trace event detection.
%   spontCa_writeEvents tblEvent -> per-cell <sbjID>.mat (minimal schema).
%   spontCa_readEvents  per-cell <sbjID>.mat -> tblEvent.
%   spontCa_manCur      interactive per-cell event-curation GUI.
%
% OUTPUT (cache/spontCa_tbl.mat)
%   tblCell   - per-cell table with traces, genotype, compartment, fs.
%   tblEvent  - long-format events: {sbjID, genotype, compartment, start, stop}.
%   fs        - sampling rate (Hz).
%
% See also: MCU_SPONTCA, SPONTCA2_METRICS.


%% ========================================================================
%  LOAD
%  ========================================================================

[tblCell, fs] = spontCa_loadXls('flgRaw', true, 'flgExclude', false, 'flgCorrF0', false);


%% ========================================================================
%  AUTO DETECT
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
%
% Re-run only when changing detection params.

paramsCyto = {'minAmp', 0.05, 'minIEI', 1.0, 'kNoise', 3.5, 'minDur', 0.4};
paramsMito = {'minAmp', 0.06, 'minIEI', 0.4, 'kNoise', 3.5, 'minDur', 0.2};

nRows = height(tblCell);
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
% Opens manCur on the in-memory tblEvent (or on the assembled tblEvent
% from the ASSEMBLE section below). Save writes per-cell bare events
% tables to man/<sbjID>.mat. After closing, re-run ASSEMBLE to refresh
% the on-disk view.
%
% NOTES (curation log):
%   Control_72, 73, 76 appear to be the same cell. Kept only 72.
%   Mito events can occur even 7 s after a cyto event (Control _33). Only
%   happened once, though.

spontCa_manCur(tblCell, tblEvent, fs);



%% ========================================================================
%  ASSEMBLE
%  ========================================================================
% Build tblEvent from disk: man/<sbjID>.mat for curated cells, falling
% back to auto/<sbjID>.mat for uncurated ones. Stands alone - no need
% to re-run AUTO DETECT or MANUAL CURATION as long as those folders
% are populated.

spDir = fileparts(which('spontCa_detect'));
tblEvent_auto = spontCa_readEvents(fullfile(spDir, 'auto'));
tblEvent_man  = spontCa_readEvents(fullfile(spDir, 'man'));

curatedCells  = unique(tblEvent_man.sbjID);
tblEvent_auto = tblEvent_auto(~ismember(tblEvent_auto.sbjID, curatedCells), :);
tblEvent      = [tblEvent_auto; tblEvent_man];
tblEvent      = sortrows(tblEvent);

% Attach genotype to each event row (sbjID -> genotype lookup).
[~, idx] = ismember(tblEvent.sbjID, tblCell.sbjID);
tblEvent.genotype = tblCell.genotype(idx);
tblEvent = movevars(tblEvent, 'genotype', 'before', 1);


%% ========================================================================
%  MANDATORY FILTER
%  ========================================================================
% (1) Drop rows with NaN start or stop - corruption from manual edits
%     that should never reach analysis.
% (2) Drop cells with zero events in either compartment. Cascading:
%     a cell that has no mito events also loses its cyto events from
%     downstream analyses (and vice-versa).
%
% These are presence/absence checks only. Amplitude-based filtering
% (e.g., sub-threshold cyto removal) is selective and lives in
% mcu_spontCa.

badRows = isnan(tblEvent.start) | isnan(tblEvent.stop);
if any(badRows)
    fprintf('Dropped %d events with NaN start or stop\n', sum(badRows));
    tblEvent(badRows, :) = [];
end

nCper = arrayfun(@(s) sum(tblEvent.sbjID == s & ...
    tblEvent.compartment == 'Cyto'), tblCell.sbjID);
nMper = arrayfun(@(s) sum(tblEvent.sbjID == s & ...
    tblEvent.compartment == 'Mito'), tblCell.sbjID);
keepCell = nCper > 0 & nMper > 0;
fprintf('Dropped %d cells with zero events in cyto or mito\n', ...
    sum(~keepCell) / 2);
tblCell  = tblCell(keepCell, :);
tblEvent = tblEvent(ismember(tblEvent.sbjID, tblCell.sbjID), :);



%% ========================================================================
%  SAVE
%  ========================================================================
% Single cached .mat for mcu_spontCa to load. Holds tblCell (with
% traces), tblEvent (minimal schema + genotype), and fs.

cacheDir = fullfile(spDir, 'cache');
if ~exist(cacheDir, 'dir'), mkdir(cacheDir); end
cachePath = fullfile(cacheDir, 'spontCa_tbl.mat');
save(cachePath, 'tblCell', 'tblEvent', 'fs', '-v7.3');
fprintf('Saved %d cells, %d events to %s\n', ...
    height(tblCell), height(tblEvent), cachePath);


%% ========================================================================
%  F0 CHECK
%  ========================================================================
% Visual check whether baseline brightness F0 differs systematically
% between genotypes per compartment. dF/F = (F - F0)/F0, so a dimmer
% baseline in one genotype mechanically inflates dF/F amps in that
% group. If F0 looks well-separated by genotype within a compartment,
% the dF/F-based amp comparison is partially confounded by
% indicator-expression bias; report ΔF (no division) alongside dF/F or
% restrict claims to slope/ratio analyses. F0 is populated by
% spontCa_loadXls when flgRaw=true (median of raw F per cell).

guiTbl_bar(tblCell, 'yVar', 'F0', 'xVar', 'compartment', 'grpVar', 'genotype');