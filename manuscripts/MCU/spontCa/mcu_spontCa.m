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
% Control_72, 73, and 76 appear exactly the same cell. Kept only 72.


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

% Remove stop/dur/int from cyto events. Cyto decays aren't biologically
% real at fs=3.
isCyto = tblEvent.compartment == 'Cyto';
tblEvent.stop(isCyto) = nan;
tblEvent.dur(isCyto)  = nan;
tblEvent.int(isCyto)  = nan;


%% ========================================================================
%  PAIR & CROSS-FLUX
%  ========================================================================
% Couple cyto and mito events and compute the cross-compartment flux
% integral for each event.
%
%   cyto row: pairIdx  -> first mito event whose start falls within
%                         [start, start + win_c2m], else NaN.
%             crossInt -> integral of mito dF/F over the same window.
%   mito row: pairIdx  -> closest preceding cyto within maxLag of mito
%                         start (the "trigger"), else NaN.
%             crossInt -> integral of cyto dF/F over [trigger.start,
%                         this.stop].
%
% Per-event transfer:
%   cyto row: T = crossInt / amp   (mito response per cyto attempt; s)
%   mito row: T = crossInt / int   (cyto input per mito response; -)

win_c2m = 5;        % cyto->mito response window (s)
maxLag  = 3;        % max cyto->mito lag to be called a trigger (s)

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
    for k = 1:numel(rowsC)
        s = startsC(k);
        i0 = max(1, round(s * fs) + 1);
        i1 = min(nSamps, round((s + win_c2m) * fs) + 1);
        if i1 >= i0
            tblEvent.pairFlux(rowsC(k)) = sum(mitoTrace(i0:i1)) * dt;
        end
        j = find(startsM >= s & startsM <= s + win_c2m, 1, 'first');
        if ~isempty(j)
            tblEvent.pairIdx(rowsC(k)) = rowsM(j);
            tblEvent.pairLag(rowsC(k)) = startsM(j) - s;
        end
    end

    % Mito -> trigger cyto, cyto input integrated from trigger to mito stop.
    for k = 1:numel(rowsM)
        ms = startsM(k);
        j = find(startsC <= ms & startsC >= ms - maxLag, 1, 'last');
        if isempty(j), continue; end
        trigStart = startsC(j);
        tblEvent.pairIdx(rowsM(k)) = rowsC(j);
        tblEvent.pairLag(rowsM(k)) = ms - trigStart;
        i0 = max(1, round(trigStart * fs) + 1);
        i1 = min(nSamps, round(stopsM(k) * fs) + 1);
        if i1 >= i0
            tblEvent.pairFlux(rowsM(k)) = sum(cytoTrace(i0:i1)) * dt;
        end
    end
end

% Per-event transfer ratio. cyto: response / attempt. mito: input / response.
isC = tblEvent.compartment == 'Cyto';
isM = tblEvent.compartment == 'Mito';
tblEvent.tf(isC) = tblEvent.pairFlux(isC) ./ tblEvent.amp(isC);
tblEvent.tf(isM) = tblEvent.pairFlux(isM) ./ tblEvent.int(isM);


%% ========================================================================
%  PER-CELL SUMMARY
%  ========================================================================

tblCell.nEvents = zeros(nRows, 1);
tblCell.rate    = zeros(nRows, 1);
tblCell.amp     = nan(nRows, 1);
tblCell.dur     = nan(nRows, 1);
tblCell.flux    = zeros(nRows, 1);
tblCell.tf   = nan(nRows, 1);
for iRow = 1:nRows
    mask = tblEvent.sbjID == tblCell.sbjID(iRow) & ...
           tblEvent.compartment == tblCell.compartment(iRow);
    tblCell.nEvents(iRow) = sum(mask);
    tblCell.rate(iRow)    = tblCell.nEvents(iRow) / recDur;
    tblCell.amp(iRow)     = mean(tblEvent.amp(mask));
    tblCell.dur(iRow)     = mean(tblEvent.dur(mask));
    tblCell.flux(iRow)    = sum(tblEvent.amp(mask)) / recDur;
    tblCell.tf(iRow)   = mean(tblEvent.tf(mask), 'omitnan');
end

% Drop cells with zero events in either compartment.
rmvId    = unique(tblCell.sbjID(tblCell.nEvents == 0));
tblLme   = tblCell(~ismember(tblCell.sbjID, rmvId), :);
tblEvent = tblEvent(~ismember(tblEvent.sbjID, rmvId), :);




%% ========================================================================
%  FIG 1E,F + S1B-E
%  ========================================================================
% Per-compartment LME + bar plot over genotype. With one observation per
% sbjID the random intercept is degenerate and the LME reduces to LM.

tblGUI_bar(tblLme, 'yVar', 'amp', 'xVar', 'compartment', 'grpVar', 'genotype');
tblGUI_bar(tblEvent, 'yVar', 'amp', 'xVar', 'compartment', 'grpVar', 'genotype');

% LME (over cells)
frml = 'meanAmp ~ genotype * compartment + (1 | sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLme, frml, 'flgPlot', false, 'verbose', true);

% LME (over events)
frml = 'amp ~ genotype * compartment + (1 | sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblEvent, frml, 'flgPlot', false, 'verbose', true);

% NOTE: currently, rate is not different between genotypes, but it could be
% biased by small events. Rate is problematic because unlike spikes, Ca
% transients are not binary. Need to find a solution. 

tblGUI_scatHist(tblEvent, 'yVar', 'amp', 'xVar', 'dur', 'grpVar', 'compartment');
tblGUI_scatHist(tblCell, 'yVar', 'amp', 'xVar', 'dur', 'grpVar', 'compartment');


%% ========================================================================
%  QC (per-cell viewer)
%  ========================================================================

spontCa_gui(tblCell, tblEvent, fs);




