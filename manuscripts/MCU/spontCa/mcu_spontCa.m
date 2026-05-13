%% mcu_spontCa.m  Spontaneous Ca2+ imaging pipeline (cyto + mito).
%
% PURPOSE
%   Read NF's SpontCa.xlsx into a long-format cell table and a long-format
%   events table, optionally curate via manCur, finalize (aggregates +
%   ETA maps + coupling). Reproduces Fig 1E,F + S1B-E and the
%   transfer-function preview.
%
% PIPELINE FILES 
%   spontCa_loadXls     Excel -> tblCell (traces + cell metadata).
%                       Picks sheet 'f' (raw F) or 'dff' via flgRaw;
%                       optionally drops experimenter-excluded cells.
%   spontCa_detect      single-trace event detection, returns a table
%   spontCa_writeEvents tblEvent -> per-cell <sbjID>.mat files
%   spontCa_readEvents  per-cell <sbjID>.mat files -> tblEvent
%   spontCa_manCur      interactive per-cell event-curation GUI
%   spontCa_gui         per-cell QC viewer


%% ========================================================================
%  LOAD
%  ========================================================================

[tblCell, fs] = spontCa_loadXls('flgRaw', true, 'flgExclude', false);

% Params
dt = 1 / fs;
nRows = height(tblCell);
nSamps = size(tblCell.trace, 2);
recDur = nSamps / fs;


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

% % Params tuned by spontCa_tune against 4 curated cells (Ctrl_01/02/03/05).
% paramsCyto = {'minAmp', 0.05, 'minIEI', 1.0, 'kNoise', 3.5, 'minDur', 0.4};
% paramsMito = {'minAmp', 0.06, 'minIEI', 0.4, 'kNoise', 3.5, 'minDur', 0.2};
% 
% chunks = cell(nRows, 1);
% for iRow = 1:nRows
%     if tblCell.compartment(iRow) == 'Cyto'
%         rowEv = spontCa_detect(tblCell.trace(iRow, :), fs, paramsCyto{:});
%     else
%         rowEv = spontCa_detect(tblCell.trace(iRow, :), fs, paramsMito{:});
%     end
%     if height(rowEv) > 0
%         rowEv.sbjID       = repmat(tblCell.sbjID(iRow),       height(rowEv), 1);
%         rowEv.compartment = repmat(tblCell.compartment(iRow), height(rowEv), 1);
%         chunks{iRow} = rowEv;
%     end
% end
% tblEvent = vertcat(chunks{~cellfun(@isempty, chunks)});
% tblEvent = tblEvent(:, ['sbjID', 'compartment', setdiff(...
%     tblEvent.Properties.VariableNames, {'sbjID','compartment'}, 'stable')]);
% 
% autoDir = fullfile(fileparts(which('spontCa_detect')), 'auto');
% spontCa_writeEvents(tblEvent, autoDir, 'backup', false);


%% ========================================================================
%  MANUAL CURATION
%  ========================================================================
% Opens manCur on the in-memory tblEvent. Save writes per-cell bare
% events tables to man/<sbjID>.mat. After closing, re-read whatever's
% on disk in man/ (if any) and merge with the auto-detection rows for
% cells the user didn't curate.

% COMMENTS:
% Control_72, 73, and 76 appear exactly the same cell. Kept only 72.

% spontCa_manCur(tblCell, tblEvent, fs);


%% ========================================================================
%  LOAD & ORGANIZE EVENTS
%  ========================================================================
% Build tblEvent from disk: man/<sbjID>.mat for curated cells, falling
% back to auto/<sbjID>.mat for uncurated ones. This block stands alone -
% no need to re-run AUTOMATIC DETECT or MANUAL CURATION as long as those
% folders are populated.

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

% Sanity check - events with zero or nan amplitude
badEvents = find(tblEvent.amp < eps | isnan(tblEvent.amp));
if ~isempty(badEvents)
    tblEvent(badEvents, :)
end

% Sanity check - events with unreasonable high amplitude
badEvents = find(tblEvent.amp > 5);
if ~isempty(badEvents)
    tblEvent(badEvents, :)
end

% Rename int -> flux for naming consistency with cell-level.
tblEvent = renamevars(tblEvent, 'int', 'flux');

% Cyto has no real decay at fs=3: the event lives in one sample. Define
% the event "integral" as amp * dt - units dF/F*s, dimensionally
% consistent with the mito event integral. Stop = start (point), dur = dt
% (one frame). 
isCyto = tblEvent.compartment == 'Cyto';
tblEvent.stop(isCyto) = tblEvent.start(isCyto);
tblEvent.dur(isCyto)  = dt;
tblEvent.flux(isCyto)  = tblEvent.amp(isCyto) * dt;


%% ========================================================================
%  FILTER & EXCLUDE
%  ========================================================================
% (1) Drop sub-threshold cyto events. minAmpCyto is set by CYTO THRESHOLD
%     section below.
% (2) Drop cells that end up with zero events in either compartment.
%     This cascades: a cell whose only cyto events were sub-threshold
%     also loses its mito events from downstream analyses.

minAmpCyto = 0;
drop = tblEvent.compartment == 'Cyto' & tblEvent.amp < minAmpCyto;
fprintf('Dropped %d cyto events with amp < %g\n', sum(drop), minAmpCyto);
tblEvent(drop, :) = [];

nCper = arrayfun(@(s) sum(tblEvent.sbjID == s & tblEvent.compartment == 'Cyto'), tblCell.sbjID);
nMper = arrayfun(@(s) sum(tblEvent.sbjID == s & tblEvent.compartment == 'Mito'), tblCell.sbjID);
keepCell = nCper > 0 & nMper > 0;
fprintf('Dropped %d cells with zero kept events in cyto or mito\n', sum(~keepCell) / 2);
tblCell  = tblCell(keepCell, :);
tblEvent = tblEvent(ismember(tblEvent.sbjID, tblCell.sbjID), :);
nRows    = height(tblCell);

%% ========================================================================
%  PAIR & CROSS-FLUX
%  ========================================================================
% Couple cyto and mito events and compute the cross-compartment flux
% integral for each event.
%
%   cyto row: pairIdx  -> first mito event whose start falls within
%                         [start, start + win_c2m], else NaN.
%             pairFlux -> integral of mito dF/F over the same window.
%   mito row: pairIdx  -> closest preceding cyto within maxLag of mito
%                         start (the "trigger"), else NaN.
%             pairFlux -> integral of cyto dF/F over [trigger.start,
%                         this.stop].
%
% Per-event transfer:
%   cyto row: T = crossInt / amp   (mito response per cyto attempt; s)
%   mito row: T = crossInt / int   (cyto input per mito response; -)

win_c2m = 10;       % cyto -> mito response window (s)
maxLag  = 5;        % max cyto -> mito lag to be called a trigger (s)

nEv = height(tblEvent);
tblEvent.pairIdx  = nan(nEv, 1);
tblEvent.pairLag  = nan(nEv, 1);
tblEvent.pairFlux = nan(nEv, 1);
tblEvent.tf        = nan(nEv, 1);

cells = unique(tblEvent.sbjID);
for iCell = 1:numel(cells)
    sid = cells(iCell);
    iC = find(tblCell.sbjID == sid & tblCell.compartment == 'Cyto', 1);
    iM = find(tblCell.sbjID == sid & tblCell.compartment == 'Mito', 1);

    cytoTrace = tblCell.trace(iC, :);
    mitoTrace = tblCell.trace(iM, :);

    rowsC = find(tblEvent.sbjID == sid & tblEvent.compartment == 'Cyto');
    rowsM = find(tblEvent.sbjID == sid & tblEvent.compartment == 'Mito');
    startsC = tblEvent.start(rowsC);
    startsM = tblEvent.start(rowsM);
    stopsM  = tblEvent.stop(rowsM);

    % Cyto -> mito response integrated over a fixed window after each cyto.
    for iEvent = 1:numel(rowsC)
        s = startsC(iEvent);
        i0 = max(1, round(s * fs) + 1);
        i1 = min(nSamps, round((s + win_c2m) * fs) + 1);
        if i1 >= i0
            tblEvent.pairFlux(rowsC(iEvent)) = sum(mitoTrace(i0:i1)) * dt;
        end
        j = find(startsM >= s & startsM <= s + win_c2m, 1, 'first');
        if ~isempty(j)
            tblEvent.pairIdx(rowsC(iEvent)) = rowsM(j);
            tblEvent.pairLag(rowsC(iEvent)) = startsM(j) - s;
        end
    end

    % Mito -> trigger cyto, cyto input integrated from trigger to mito stop.
    for iEvent = 1:numel(rowsM)
        ms = startsM(iEvent);
        j = find(startsC <= ms & startsC >= ms - maxLag, 1, 'last');
        if isempty(j), continue; end
        trigStart = startsC(j);
        tblEvent.pairIdx(rowsM(iEvent)) = rowsC(j);
        tblEvent.pairLag(rowsM(iEvent)) = ms - trigStart;
        i0 = max(1, round(trigStart * fs) + 1);
        i1 = min(nSamps, round(stopsM(iEvent) * fs) + 1);
        if i1 >= i0
            tblEvent.pairFlux(rowsM(iEvent)) = sum(cytoTrace(i0:i1)) * dt;
        end
    end
end

% Per-event transfer ratio. Uniform across compartments because cyto.flux
% (= amp*dt) and mito.flux (= event integral) share units of dF/F*s.
%   cyto row: T = mito response integral / cyto integral
%   mito row: T = cyto input integral    / mito integral
tblEvent.tf = tblEvent.pairFlux ./ tblEvent.flux;

% Clip negative pairFlux / tf to zero. These are events whose cross-window
% integrated the opposite compartment around its noise floor and came out
% slightly below zero. They carry no biological information.
tblEvent.pairFlux(tblEvent.pairFlux < 0) = 0;
tblEvent.tf(tblEvent.tf < 0) = 0;


%% ========================================================================
%  PER-CELL SUMMARY
%  ========================================================================

% tblEvent.flux has units of dF/F*s (event integral). tblCell.fluxRate
% is the sum of event flux over recording duration, units dF/F per s.

tblCell.nEvents  = zeros(nRows, 1);
tblCell.rate     = zeros(nRows, 1);
tblCell.amp      = nan(nRows, 1);
tblCell.dur      = nan(nRows, 1);
tblCell.fluxRate = zeros(nRows, 1);
tblCell.tf       = nan(nRows, 1);
for iRow = 1:nRows
    mask = tblEvent.sbjID == tblCell.sbjID(iRow) & ...
           tblEvent.compartment == tblCell.compartment(iRow);
    tblCell.nEvents(iRow)  = sum(mask);
    tblCell.rate(iRow)     = tblCell.nEvents(iRow) / recDur;
    tblCell.amp(iRow)      = mean(tblEvent.amp(mask));
    tblCell.dur(iRow)      = mean(tblEvent.dur(mask));
    tblCell.fluxRate(iRow) = sum(tblEvent.flux(mask)) / recDur;
    tblCell.tf(iRow)       = mean(tblEvent.tf(mask), 'omitnan');
end





%% ========================================================================
%  INSPECT RESULTS
%  ========================================================================

mode = 'cell';

if strcmp(mode, 'event')
    tbl = tblEvent;
    vars = {'amp', 'flux', 'pairFlux', 'tf'};
elseif strcmp(mode, 'cell')
    tbl = tblCell;
    vars = {'amp', 'fluxRate', 'rate', 'tf'};
end

yVar = vars{3};
xVar = vars{2};

% Bar
tblGUI_bar(tbl, 'yVar', yVar, 'xVar', 'compartment', 'grpVar', 'genotype');

% Across compartments
tblGUI_scatHist(tbl, 'yVar', yVar, 'xVar', xVar, 'grpVar', 'compartment');

% Per compartment
cmp = 'Mito';
tblGUI_scatHist(tbl(tbl.compartment == cmp, :), 'yVar', yVar, 'xVar', xVar, 'grpVar', 'genotype');

cmp = 'Cyto';
tblGUI_scatHist(tbl(tbl.compartment == cmp, :), 'yVar', yVar, 'xVar', xVar, 'grpVar', 'genotype');

% Cyto versus Mito. For 'cell' mode each cell has one cyto + one mito
% row; pair them by position. For 'event' mode pair via pairIdx (cyto
% events that have a paired mito event).
var = 'amp';
if strcmp(mode, 'cell')
    isC = tbl.compartment == 'Cyto';
    isM = tbl.compartment == 'Mito';
    tblSub = tbl(isC, {'sbjID', 'genotype'});
    tblSub.cyto = tbl.(var)(isC);
    tblSub.mito = tbl.(var)(isM);
else
    isCyto = tbl.compartment == 'Cyto' & ~isnan(tbl.pairIdx);
    tblSub = tbl(isCyto, {'sbjID', 'genotype'});
    tblSub.cyto = tbl.(var)(isCyto);
    tblSub.mito = tbl.(var)(tbl.pairIdx(isCyto));
end
tblGUI_scatHist(tblSub, 'xVar', 'cyto', 'yVar', 'mito', 'grpVar', 'genotype');


%  LME
frml = [vars{1}, ' ~ genotype * compartment + (1 | sbjID)'];
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tbl, frml, 'flgPlot', false, 'verbose', true);







%% ========================================================================
%  CYTO THRESHOLD (decision)
%  ========================================================================
% Compare two distributions of cyto amps in Control cells:
%   - all cyto events (every detected cyto, regardless of triggering)
%   - triggering cytos (those paired with a mito event via pairIdx)
% The trigger distribution is the empirical "real cyto" distribution;
% below its lower tail is detector noise that inflates cyto rate and
% drags T toward zero. Pick minAmpCyto from its lower tail and set it
% at the top of LOAD & ORGANIZE, then re-run from there.

mask    = tblEvent.compartment == 'Cyto' & tblEvent.genotype == 'Control';
ampAll  = tblEvent.amp(mask);
ampTrig = tblEvent.amp(mask & ~isnan(tblEvent.pairIdx));
thrSugg = prctile(ampTrig, 5);

figure('Name', 'Cyto threshold (Control)', 'Color', 'w');
hold on
histogram(ampAll,  'Normalization', 'pdf', 'FaceAlpha', 0.4, ...
    'DisplayName', sprintf('all cytos (n=%d)', numel(ampAll)));
histogram(ampTrig, 'Normalization', 'pdf', 'FaceAlpha', 0.4, ...
    'DisplayName', sprintf('triggers (n=%d)', numel(ampTrig)));
xline(thrSugg, 'r--', sprintf('  5%% trig = %.3f', thrSugg), ...
    'LabelOrientation', 'horizontal');
xlabel('cyto amp (dF/F)'); ylabel('pdf');
legend('Location', 'best');
title(sprintf('Suggested minAmpCyto = %.3f', thrSugg));
hold off


%% ========================================================================
%  QC (per-cell viewer)
%  ========================================================================

spontCa_gui(tblCell, tblEvent, fs);




