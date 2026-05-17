function [tblCell, tblEvent, info] = spontCa_filter(tblCell, tblEvent, varargin)
% SPONTCA_FILTER  Selective event/cell filter with summary report.
%
% [tblCell, tblEvent, info] = SPONTCA_FILTER(tblCell, tblEvent, ...) applies
% per-compartment amp / snr / lag filters and a per-cell minimum-events
% rule. Returns the filtered tables and an info struct with per-rule
% drop counts. Prints a per-rule summary by default.
%
% Does NOT call spontCa2_metrics. Caller is responsible for running
% Pass 2 (mode='cellOnly') afterwards to refresh pairing and cell-level
% aggregates.
%
% OPTIONAL (Name-Value):
%   'minAmp'    - [Cyto; Mito] min event amp.        Default [0; 0].
%   'minSNR'    - [Cyto; Mito] min event snr.        Default [0; 0].
%   'minEvents' - [Cyto; Mito] min events per cell.  Default [0; 0].
%   'flgPair'   - 'all' (default) | 'paired' | 'unpaired'. Filters on
%                 the tblEvent.paired logical (built upstream by
%                 spontCa2_metrics from winPair). 'paired' keeps
%                 paired==true, 'unpaired' keeps paired==false, 'all'
%                 keeps everything.
%   'verbose'   - print per-rule counts to stdout.   Default true.
%
% INFO struct fields:
%   .amp        - struct with .Cyto, .Mito (events dropped by amp rule).
%   .snr        - struct with .Cyto, .Mito (events dropped by snr rule).
%   .pair       - struct with .Cyto, .Mito (events dropped by flgPair).
%   .minEvents  - scalar (cell pairs dropped by minEvents rule).
%   .nEventsBefore, .nEventsAfter - struct with .Cyto, .Mito.
%   .nCellsBefore, .nCellsAfter   - scalar (cell pairs).
%
% See also: SPONTCA2_METRICS.


%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tblCell',  @istable);
addRequired(p, 'tblEvent', @istable);
addParameter(p, 'minAmp',    [0; 0],       @(x) isnumeric(x) && numel(x) == 2);
addParameter(p, 'minSNR',    [0; 0],       @(x) isnumeric(x) && numel(x) == 2);
addParameter(p, 'minEvents', [0; 0],       @(x) isnumeric(x) && numel(x) == 2);
addParameter(p, 'flgPair',   'all',        @(x) any(strcmpi(x, {'all', 'paired', 'unpaired'})));
addParameter(p, 'verbose',   true,         @islogical);
parse(p, tblCell, tblEvent, varargin{:});
P = p.Results;

minAmp    = P.minAmp(:);
minSNR    = P.minSNR(:);
minEvents = P.minEvents(:);


%% ========================================================================
%  BEFORE COUNTS
%  ========================================================================

isC = tblEvent.compartment == 'Cyto';
isM = tblEvent.compartment == 'Mito';

info.nEventsBefore = struct('Cyto', sum(isC), 'Mito', sum(isM));
info.nCellsBefore  = height(tblCell) / 2;   % cell pairs


%% ========================================================================
%  AMP RULE (per-compartment)
%  ========================================================================

keepAmp = (isC & tblEvent.amp >= minAmp(1)) | (isM & tblEvent.amp >= minAmp(2));
dropAmpCyto = sum(isC & ~keepAmp);
dropAmpMito = sum(isM & ~keepAmp);
tblEvent = tblEvent(keepAmp, :);
info.amp = struct('Cyto', dropAmpCyto, 'Mito', dropAmpMito);


%% ========================================================================
%  SNR RULE (per-compartment)
%  ========================================================================

isC = tblEvent.compartment == 'Cyto';
isM = tblEvent.compartment == 'Mito';

if ismember('snr', tblEvent.Properties.VariableNames)
    keepSnr = (isC & tblEvent.snr >= minSNR(1)) | (isM & tblEvent.snr >= minSNR(2));
else
    keepSnr = true(height(tblEvent), 1);
end
dropSnrCyto = sum(isC & ~keepSnr);
dropSnrMito = sum(isM & ~keepSnr);
tblEvent = tblEvent(keepSnr, :);
info.snr = struct('Cyto', dropSnrCyto, 'Mito', dropSnrMito);


%% ========================================================================
%  PAIRED RULE (uses tblEvent.paired set upstream by spontCa2_metrics)
%  ========================================================================

isC = tblEvent.compartment == 'Cyto';
isM = tblEvent.compartment == 'Mito';

switch lower(P.flgPair)
    case 'all'
        keepPair = true(height(tblEvent), 1);
    case 'paired'
        keepPair = tblEvent.paired;
    case 'unpaired'
        keepPair = ~tblEvent.paired;
end
dropPairCyto = sum(isC & ~keepPair);
dropPairMito = sum(isM & ~keepPair);
tblEvent = tblEvent(keepPair, :);
info.pair = struct('Cyto', dropPairCyto, 'Mito', dropPairMito);


%% ========================================================================
%  MIN-EVENTS RULE (per cell pair)
%  ========================================================================

cellSbj = unique(tblCell.sbjID);
nCper = arrayfun(@(s) sum(tblEvent.sbjID == s & ...
    tblEvent.compartment == 'Cyto'), cellSbj);
nMper = arrayfun(@(s) sum(tblEvent.sbjID == s & ...
    tblEvent.compartment == 'Mito'), cellSbj);
keepSbj = cellSbj(nCper >= minEvents(1) & nMper >= minEvents(2));

dropCells = numel(cellSbj) - numel(keepSbj);
tblCell  = tblCell(ismember(tblCell.sbjID, keepSbj), :);
tblEvent = tblEvent(ismember(tblEvent.sbjID, keepSbj), :);
info.minEvents = dropCells;


%% ========================================================================
%  AFTER COUNTS + REPORT
%  ========================================================================

isC = tblEvent.compartment == 'Cyto';
isM = tblEvent.compartment == 'Mito';
info.nEventsAfter = struct('Cyto', sum(isC), 'Mito', sum(isM));
info.nCellsAfter  = height(tblCell) / 2;

if P.verbose
    pctC = @(d, n) 100 * d / max(1, n);
    nC0 = info.nEventsBefore.Cyto;
    nM0 = info.nEventsBefore.Mito;
    fprintf('=== spontCa_filter ===\n');
    fprintf('  amp gate         Cyto: -%d / %d (%.1f%%)   Mito: -%d / %d (%.1f%%)\n', ...
        info.amp.Cyto, nC0, pctC(info.amp.Cyto, nC0), ...
        info.amp.Mito, nM0, pctC(info.amp.Mito, nM0));
    fprintf('  snr gate         Cyto: -%d / %d (%.1f%%)   Mito: -%d / %d (%.1f%%)\n', ...
        info.snr.Cyto, nC0, pctC(info.snr.Cyto, nC0), ...
        info.snr.Mito, nM0, pctC(info.snr.Mito, nM0));
    if ~strcmpi(P.flgPair, 'all')
        fprintf('  pair gate        Cyto: -%d / %d (%.1f%%)   Mito: -%d / %d (%.1f%%)   [mode=%s]\n', ...
            info.pair.Cyto, nC0, pctC(info.pair.Cyto, nC0), ...
            info.pair.Mito, nM0, pctC(info.pair.Mito, nM0), ...
            P.flgPair);
    end
    fprintf('  minEvents drop   %d / %d cell pairs\n', ...
        dropCells, info.nCellsBefore);
    fprintf('  Remaining: %d cyto events, %d mito events, %d cell pairs.\n', ...
        info.nEventsAfter.Cyto, info.nEventsAfter.Mito, info.nCellsAfter);
end

end     % SPONTCA_FILTER
