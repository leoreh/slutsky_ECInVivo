%% run_characterize.m  Diagnostic characterization of spontCa traces + detector.
%
% PURPOSE
%   Read-only audit of the spontaneous-Ca2+ pipeline. Loads the long-format
%   table via spontCa_load, runs spontCa_detect per row with the canonical
%   per-compartment params, and writes per-cell and per-event statistics
%   plus a human-readable summary into _diag/. The goal is to surface
%   cells that the detector likely mishandles (plateau-heavy traces,
%   suspiciously few events, bimodal amplitude distributions) so curation
%   can be prioritised.
%
% INPUTS
%   None. Reads the canonical SpontCa.xlsx via spontCa_load.
%
% OUTPUTS (written to manuscripts/MCU/spontCa/_diag/)
%   cellStats.csv   - one row per (sbjID, compartment).
%   eventStats.csv  - one row per detected event.
%   summary.md      - human-readable report.
%
% USAGE
%   Run from anywhere; the script resolves its own output directory
%   relative to its own location.
%       >> run('D:/Code/slutsky_ECInVivo/manuscripts/MCU/spontCa/_diag/run_characterize.m')
%   The pipeline files (spontCa_load, spontCa_detect) must be on the path
%   or accessible from the script's parent directory.
%
% NOTES
%   - Does NOT modify any pipeline file.
%   - Does NOT depend on spontCa_finalize; detection is inlined.
%   - Per-cell event columns (start/stop/amp/dur/int) are populated in a
%     local copy of the table only; the caller's workspace is not touched.

%% ========================================================================
%  PATHS + OUTPUT DIR
%  ========================================================================

thisFile  = mfilename('fullpath');
if isempty(thisFile)
    % Fallback when run via "Run Section" with an unsaved buffer.
    diagDir = fullfile('D:', 'Code', 'slutsky_ECInVivo', ...
        'manuscripts', 'MCU', 'spontCa', '_diag');
else
    diagDir = fileparts(thisFile);
end
parentDir = fileparts(diagDir);

% Make sure the pipeline directory is on path so spontCa_load / _detect
% resolve regardless of caller cwd.
if exist(parentDir, 'dir') && ~contains(lower(path), lower(parentDir))
    addpath(parentDir);
end

if ~exist(diagDir, 'dir')
    mkdir(diagDir);
end

fprintf('[run_characterize] Output dir: %s\n', diagDir);


%% ========================================================================
%  LOAD
%  ========================================================================

[tbl, fs] = spontCa_load();
nRows = height(tbl);
nT    = size(tbl.trace, 2);
tVec  = (0:nT-1) / fs;
recDur = (nT - 1) / fs;
fprintf('[run_characterize] Loaded %d rows, fs=%.3f Hz, nT=%d (%.1f s)\n', ...
    nRows, fs, nT, recDur);


%% ========================================================================
%  DETECTOR PARAMS (canonical, from mcu_spontCa.m)
%  ========================================================================

paramsCyto = {'minAmp', 0.05, 'minIEI', 1.0, 'kNoise', 3.5, 'minDur', 0.4};
paramsMito = {'minAmp', 0.03, 'minIEI', 1.0, 'kNoise', 3.5, 'minDur', 0.4};


%% ========================================================================
%  ALLOCATE PER-ROW EVENT CELLS
%  ========================================================================

tbl.start = cell(nRows, 1);
tbl.stop  = cell(nRows, 1);
tbl.amp   = cell(nRows, 1);
tbl.dur   = cell(nRows, 1);
tbl.int   = cell(nRows, 1);


%% ========================================================================
%  PER-CELL DIAGNOSTIC ACCUMULATORS
%  ========================================================================

sbjStr      = strings(nRows, 1);
genStr      = strings(nRows, 1);
cmpStr      = strings(nRows, 1);
traceMed    = nan(nRows, 1);
traceMad    = nan(nRows, 1);
plateauFrac = nan(nRows, 1);
dNoise      = nan(nRows, 1);
nEventsVec  = zeros(nRows, 1);
meanAmpVec  = nan(nRows, 1);
meanDurVec  = nan(nRows, 1);
meanIEIVec  = nan(nRows, 1);
fluxVec     = nan(nRows, 1);
fluxIntVec  = nan(nRows, 1);
plateauHeavy = false(nRows, 1);
fewEvents    = false(nRows, 1);
bimodalAmp   = false(nRows, 1);

% Per-event accumulator (build as a cell-of-rows then vertcat at the end).
evRowsAll = cell(nRows, 1);


%% ========================================================================
%  MAIN LOOP : DETECT + CHARACTERIZE
%  ========================================================================

minPlatSec = 5;     % plateau run must be >= 5 s
minPlatSmp = max(1, round(minPlatSec * fs));

for iR = 1:nRows
    trace = tbl.trace(iR, :);
    sbj   = char(tbl.sbjID(iR));
    gen   = char(tbl.genotype(iR));
    cmp   = char(tbl.compartment(iR));
    sbjStr(iR) = sbj;
    genStr(iR) = gen;
    cmpStr(iR) = cmp;

    % Defensive: all-NaN trace.
    if all(isnan(trace))
        warning('run_characterize:nanTrace', ...
            'Row %d (%s/%s) trace is all NaN, skipping.', iR, sbj, cmp);
        continue;
    end

    % Trace-level robust stats.
    traceVal = trace(~isnan(trace));
    if isempty(traceVal)
        continue;
    end
    medT = median(traceVal);
    madT = 1.4826 * median(abs(traceVal - medT));
    traceMed(iR) = medT;
    traceMad(iR) = madT;

    % Plateau fraction: runs of (trace > med + 3*MAD) >= 5 s.
    thrPlat = medT + 3 * madT;
    isHi = trace > thrPlat;
    isHi(isnan(trace)) = false;
    runDur = 0;
    runs   = [];
    inRun  = false;
    runStart = 0;
    for k = 1:numel(isHi)
        if isHi(k)
            if ~inRun
                inRun = true;
                runStart = k;
            end
        else
            if inRun
                runs(end+1, :) = [runStart, k-1]; %#ok<AGROW>
                inRun = false;
            end
        end
    end
    if inRun
        runs(end+1, :) = [runStart, numel(isHi)]; %#ok<AGROW>
    end
    if ~isempty(runs)
        lens = runs(:, 2) - runs(:, 1) + 1;
        runDur = sum(lens(lens >= minPlatSmp)) / fs;
    end
    if recDur > 0
        plateauFrac(iR) = runDur / recDur;
    end

    % Derivative noise (per-cell): MAD of diff(movmean(trace, 3)).
    sm = movmean(trace, 3, 'omitnan');
    dsm = diff(sm);
    dsm = dsm(~isnan(dsm));
    if ~isempty(dsm)
        dNoise(iR) = 1.4826 * median(abs(dsm - median(dsm)));
    end

    % Run detector with compartment-appropriate params.
    if strcmpi(cmp, 'Cyto')
        params = paramsCyto;
    else
        params = paramsMito;
    end
    try
        ev = spontCa_detect(trace, fs, params{:});
    catch ME
        warning('run_characterize:detectFail', ...
            'spontCa_detect failed on row %d (%s/%s): %s', ...
            iR, sbj, cmp, ME.message);
        continue;
    end

    nEv = numel(ev.start);
    tbl.start{iR} = ev.start;
    tbl.stop{iR}  = ev.stop;
    tbl.amp{iR}   = ev.amp;
    tbl.dur{iR}   = ev.dur;
    tbl.int{iR}   = ev.int;
    nEventsVec(iR) = nEv;

    if nEv == 0
        % Nothing more to compute; flags below still meaningful.
    else
        amps = ev.amp(:);
        durs = ev.dur(:);
        starts = ev.start(:);
        ints   = ev.int(:);

        meanAmpVec(iR) = mean(amps, 'omitnan');
        meanDurVec(iR) = mean(durs, 'omitnan');
        if nEv >= 2
            meanIEIVec(iR) = mean(diff(starts), 'omitnan');
        end

        % Flux / fluxInt definitions: per-cell totals normalised by
        % recording duration. This matches the spirit of "rate-weighted"
        % summaries used downstream without requiring spontCa_finalize.
        if recDur > 0
            fluxVec(iR)    = sum(amps, 'omitnan') / recDur;
            fluxIntVec(iR) = sum(ints, 'omitnan') / recDur;
        end

        % Bimodality flag: skip if < 6 events.
        if nEv >= 6
            ampsSorted = sort(amps);
            gaps = diff(ampsSorted);
            ampRange = ampsSorted(end) - ampsSorted(1);
            if ampRange > 0
                [maxGap, gIdx] = max(gaps);
                gapAmp = 0.5 * (ampsSorted(gIdx) + ampsSorted(gIdx + 1));
                if (maxGap / ampRange) > 0.3 && gapAmp > median(amps)
                    bimodalAmp(iR) = true;
                end
            end
        end

        % Per-event table rows.
        ieiCol = nan(nEv, 1);
        if nEv >= 2
            ieiCol(2:end) = diff(starts);
        end
        evRowsAll{iR} = table( ...
            repmat(string(sbj), nEv, 1), ...
            repmat(string(gen), nEv, 1), ...
            repmat(string(cmp), nEv, 1), ...
            (1:nEv)', ...
            starts, ...
            amps, ...
            durs, ...
            ints, ...
            ieiCol, ...
            'VariableNames', {'sbjID', 'genotype', 'compartment', ...
                              'eventIdx', 'peakTime', 'amp', 'dur', ...
                              'intg', 'iei'});
    end

    % Flags (computed regardless of nEv).
    plateauHeavy(iR) = (plateauFrac(iR) > 0.5);
    fewEvents(iR)    = (nEv < 5);
end


%% ========================================================================
%  ASSEMBLE PER-CELL TABLE
%  ========================================================================

cellStats = table( ...
    sbjStr, genStr, cmpStr, ...
    traceMed, traceMad, plateauFrac, dNoise, ...
    nEventsVec, meanAmpVec, meanDurVec, meanIEIVec, ...
    fluxVec, fluxIntVec, ...
    plateauHeavy, fewEvents, bimodalAmp, ...
    'VariableNames', {'sbjID', 'genotype', 'compartment', ...
                      'traceMed', 'traceMad', 'plateauFrac', 'dNoise', ...
                      'nEvents', 'meanAmp', 'meanDur', 'meanIEI', ...
                      'flux', 'fluxInt', ...
                      'plateauHeavy', 'fewEvents', 'bimodalAmp'});

writetable(cellStats, fullfile(diagDir, 'cellStats.csv'));
fprintf('[run_characterize] Wrote %s\n', fullfile(diagDir, 'cellStats.csv'));


%% ========================================================================
%  ASSEMBLE PER-EVENT TABLE
%  ========================================================================

evRowsAll = evRowsAll(~cellfun(@isempty, evRowsAll));
if isempty(evRowsAll)
    eventStats = table( ...
        strings(0,1), strings(0,1), strings(0,1), ...
        zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
        'VariableNames', {'sbjID', 'genotype', 'compartment', ...
                          'eventIdx', 'peakTime', 'amp', 'dur', ...
                          'intg', 'iei'});
else
    eventStats = vertcat(evRowsAll{:});
end
writetable(eventStats, fullfile(diagDir, 'eventStats.csv'));
fprintf('[run_characterize] Wrote %s (%d events)\n', ...
    fullfile(diagDir, 'eventStats.csv'), height(eventStats));


%% ========================================================================
%  PER GENOTYPE x COMPARTMENT SUMMARIES
%  ========================================================================

genoLevels = unique(genStr(genStr ~= ""), 'stable');
cmpLevels  = unique(cmpStr(cmpStr ~= ""), 'stable');

% Helper: 5-number summary -> struct.
fiveNum = @(x) struct( ...
    'n',    sum(~isnan(x)), ...
    'med',  median(x, 'omitnan'), ...
    'q1',   prctile(x, 25), ...
    'q3',   prctile(x, 75), ...
    'mn',   min(x, [], 'omitnan'), ...
    'mx',   max(x, [], 'omitnan'));

distSumm = struct();
nCellsGC  = zeros(numel(genoLevels), numel(cmpLevels));
nEventsGC = zeros(numel(genoLevels), numel(cmpLevels));
meanEvPerCell = nan(numel(genoLevels), numel(cmpLevels));
nEvtsIQR  = cell(numel(genoLevels), numel(cmpLevels));

for ig = 1:numel(genoLevels)
    gName = char(genoLevels(ig));
    for ic = 1:numel(cmpLevels)
        cName = char(cmpLevels(ic));
        fldG  = matlab.lang.makeValidName(gName);
        fldC  = matlab.lang.makeValidName(cName);

        rowSel = (genStr == gName) & (cmpStr == cName);
        evSel  = (eventStats.genotype == gName) & ...
                 (eventStats.compartment == cName);

        nCellsGC(ig, ic)      = sum(rowSel);
        nEventsGC(ig, ic)     = sum(evSel);
        if nCellsGC(ig, ic) > 0
            meanEvPerCell(ig, ic) = mean(nEventsVec(rowSel), 'omitnan');
        end

        distSumm.(fldG).(fldC).amp   = fiveNum(eventStats.amp(evSel));
        distSumm.(fldG).(fldC).dur   = fiveNum(eventStats.dur(evSel));
        distSumm.(fldG).(fldC).iei   = fiveNum(eventStats.iei(evSel));
        distSumm.(fldG).(fldC).plat  = fiveNum(plateauFrac(rowSel));

        % nEvents IQR for the "expected per-cell event count range" section.
        nEvtsIQR{ig, ic} = fiveNum(nEventsVec(rowSel));
    end
end


%% ========================================================================
%  WRITE summary.md
%  ========================================================================

mdPath = fullfile(diagDir, 'summary.md');
fid = fopen(mdPath, 'w');
if fid < 0
    error('run_characterize:openFail', ...
        'Cannot open %s for writing.', mdPath);
end
cleanupObj = onCleanup(@() fclose(fid));

fprintf(fid, '# spontCa diagnostic summary\n\n');
fprintf(fid, 'Generated by `run_characterize.m`.\n\n');
fprintf(fid, 'Input: SpontCa.xlsx via `spontCa_load`. fs=%.3f Hz, nT=%d, ', ...
    fs, nT);
fprintf(fid, 'recording duration=%.1f s.\n\n', recDur);

% Detector setup.
fprintf(fid, '## Detector setup\n\n');
fprintf(fid, 'Per-row detection via `spontCa_detect` with compartment-specific params.\n\n');
fprintf(fid, 'Cyto: ');
for k = 1:2:numel(paramsCyto)
    fprintf(fid, '`%s`=%g ', paramsCyto{k}, paramsCyto{k+1});
end
fprintf(fid, '\n\n');
fprintf(fid, 'Mito: ');
for k = 1:2:numel(paramsMito)
    fprintf(fid, '`%s`=%g ', paramsMito{k}, paramsMito{k+1});
end
fprintf(fid, '\n\n');

% Counts per genotype x compartment.
fprintf(fid, '## Counts per genotype x compartment\n\n');
fprintf(fid, '| genotype | compartment | n cells | total events | mean events / cell |\n');
fprintf(fid, '|---|---|---|---|---|\n');
for ig = 1:numel(genoLevels)
    for ic = 1:numel(cmpLevels)
        fprintf(fid, '| %s | %s | %d | %d | %.2f |\n', ...
            char(genoLevels(ig)), char(cmpLevels(ic)), ...
            nCellsGC(ig, ic), nEventsGC(ig, ic), meanEvPerCell(ig, ic));
    end
end
fprintf(fid, '\n');

% Flagged cells.
fprintf(fid, '## Flagged cells\n\n');
isFlagged = plateauHeavy | fewEvents | bimodalAmp;
flagIdx = find(isFlagged);
if isempty(flagIdx)
    fprintf(fid, 'No cells were flagged.\n\n');
else
    fprintf(fid, '| sbjID | genotype | compartment | flags | nEvents | plateauFrac |\n');
    fprintf(fid, '|---|---|---|---|---|---|\n');
    for k = 1:numel(flagIdx)
        iR = flagIdx(k);
        flagList = strings(0, 1);
        if plateauHeavy(iR), flagList(end+1) = "plateauHeavy"; end %#ok<AGROW>
        if fewEvents(iR),    flagList(end+1) = "fewEvents"; end %#ok<AGROW>
        if bimodalAmp(iR),   flagList(end+1) = "bimodalAmp"; end %#ok<AGROW>
        fprintf(fid, '| %s | %s | %s | %s | %d | %.2f |\n', ...
            char(sbjStr(iR)), char(genStr(iR)), char(cmpStr(iR)), ...
            char(strjoin(flagList, ', ')), nEventsVec(iR), plateauFrac(iR));
    end
    fprintf(fid, '\n');
end

% Distribution summaries.
fprintf(fid, '## Distribution summaries (per genotype x compartment)\n\n');
fprintf(fid, 'Reported as median [Q1, Q3] (min, max); n is # observations.\n\n');
metricLabels = {'amp', 'dur', 'iei', 'plat'};
metricTitles = {'Amplitude (dF/F)', 'Duration (s)', 'IEI (s)', 'Plateau fraction'};
for ig = 1:numel(genoLevels)
    for ic = 1:numel(cmpLevels)
        gName = char(genoLevels(ig));
        cName = char(cmpLevels(ic));
        fldG  = matlab.lang.makeValidName(gName);
        fldC  = matlab.lang.makeValidName(cName);
        fprintf(fid, '### %s / %s\n\n', gName, cName);
        fprintf(fid, '| metric | n | median | Q1 | Q3 | min | max |\n');
        fprintf(fid, '|---|---|---|---|---|---|---|\n');
        for m = 1:numel(metricLabels)
            s = distSumm.(fldG).(fldC).(metricLabels{m});
            fprintf(fid, '| %s | %d | %.4g | %.4g | %.4g | %.4g | %.4g |\n', ...
                metricTitles{m}, s.n, s.med, s.q1, s.q3, s.mn, s.mx);
        end
        fprintf(fid, '\n');
    end
end

% Cells to curate first.
fprintf(fid, '## Cells to curate first\n\n');
fprintf(fid, 'These cells carry the most weight in changing downstream results.\n');
fprintf(fid, 'Two priority criteria:\n\n');
fprintf(fid, '- `plateauHeavy` AND `fewEvents` simultaneously: the detector is likely missing real activity hidden in plateaus.\n');
fprintf(fid, '- `bimodalAmp` AND high event count (nEvents >= median over flagged set): amplitude clustering suggests mixed populations the detector lumps together.\n\n');

priority1 = find(plateauHeavy & fewEvents);
medN = median(nEventsVec(bimodalAmp), 'omitnan');
if isnan(medN), medN = 0; end
priority2 = find(bimodalAmp & (nEventsVec >= medN) & ~ismember((1:nRows)', priority1));

if isempty(priority1) && isempty(priority2)
    fprintf(fid, 'No cells meet the priority criteria.\n\n');
else
    fprintf(fid, '| priority | sbjID | genotype | compartment | nEvents | plateauFrac | flags |\n');
    fprintf(fid, '|---|---|---|---|---|---|---|\n');
    for k = 1:numel(priority1)
        iR = priority1(k);
        flagList = strings(0, 1);
        if plateauHeavy(iR), flagList(end+1) = "plateauHeavy"; end %#ok<AGROW>
        if fewEvents(iR),    flagList(end+1) = "fewEvents"; end %#ok<AGROW>
        if bimodalAmp(iR),   flagList(end+1) = "bimodalAmp"; end %#ok<AGROW>
        fprintf(fid, '| 1 | %s | %s | %s | %d | %.2f | %s |\n', ...
            char(sbjStr(iR)), char(genStr(iR)), char(cmpStr(iR)), ...
            nEventsVec(iR), plateauFrac(iR), ...
            char(strjoin(flagList, ', ')));
    end
    for k = 1:numel(priority2)
        iR = priority2(k);
        flagList = strings(0, 1);
        if plateauHeavy(iR), flagList(end+1) = "plateauHeavy"; end %#ok<AGROW>
        if fewEvents(iR),    flagList(end+1) = "fewEvents"; end %#ok<AGROW>
        if bimodalAmp(iR),   flagList(end+1) = "bimodalAmp"; end %#ok<AGROW>
        fprintf(fid, '| 2 | %s | %s | %s | %d | %.2f | %s |\n', ...
            char(sbjStr(iR)), char(genStr(iR)), char(cmpStr(iR)), ...
            nEventsVec(iR), plateauFrac(iR), ...
            char(strjoin(flagList, ', ')));
    end
    fprintf(fid, '\n');
end

% Expected per-cell event count range.
fprintf(fid, '## Expected per-cell event count range\n\n');
fprintf(fid, 'IQR of nEvents within each genotype x compartment cell.\n');
fprintf(fid, 'Cells that fall outside the IQR by a wide margin are worth a manual look.\n\n');
fprintf(fid, '| genotype | compartment | n cells | median | Q1 | Q3 | min | max |\n');
fprintf(fid, '|---|---|---|---|---|---|---|---|\n');
for ig = 1:numel(genoLevels)
    for ic = 1:numel(cmpLevels)
        s = nEvtsIQR{ig, ic};
        fprintf(fid, '| %s | %s | %d | %.1f | %.1f | %.1f | %.1f | %.1f |\n', ...
            char(genoLevels(ig)), char(cmpLevels(ic)), ...
            s.n, s.med, s.q1, s.q3, s.mn, s.mx);
    end
end
fprintf(fid, '\n');

% Done.
clear cleanupObj;
fprintf('[run_characterize] Wrote %s\n', mdPath);
fprintf('[run_characterize] Done.\n');

% EOF
