function [outArr, report] = spontCa_extractRaw(varargin)
% SPONTCA_EXTRACTRAW Build sponCa_raw.xlsx from NetaF raw fluorescence files.
%
%   [OUTARR, REPORT] = SPONTCA_EXTRACTRAW(...) reads per-cell raw
%   fluorescence from the two NetaF source workbooks, attaches the
%   experimenter's manual exclusion flag from CellsExcluded.xlsx, and
%   writes one paired cyto+mito raw-trace workbook.
%
%   Cell IDs match the experimenter's labeling (and SpontCa.xlsx): the k-th
%   cell in a raw sheet (skipping the gap column that separates two
%   recording blocks) gets ID = k for file 1, ID = N_file1 + k for file 2.
%   So `Control 65` is the last cell of file 1 ctrl; `Control 66` is the
%   first cell of file 2 ctrl.
%
%   The script also reproduces Boaz's dF/F (msbackadj -> clip negatives ->
%   percentile-normalize, with a lag search to recover frame cropping) and
%   correlates it against the *dff sheet column. Correlations confirm
%   raw <-> dff column correspondence on a per-cell basis.
%
%   FORMULA (from Ca_EventDetector.m, lines 28-64):
%       dataAdjust    = msbackadj(t, raw, 'WindowSize', 0.1*nT)
%       dataAdjust(dataAdjust < 0) = 0
%       F0            = prctile(raw, p)        % p=20 cyto, p=5 mito
%       dff           = dataAdjust / F0
%
%   OUTPUT XLSX LAYOUT (rows 1-indexed):
%       row 1   : genotype     (Control / MCU-KO), repeated over each pair
%       row 2   : cell ID      (numeric, above Cyto col; empty above Mito)
%       row 3   : excluded     (logical, the experimenter's flag from
%                               CellsExcluded.xlsx; repeated over the pair)
%       row 4   : compartment  (Cyto / Mito)
%       row 5+  : data         col 1 = time (s); col 2+ = raw F values
%
%   OPTIONAL (Name-Value):
%       'srcFiles'      - (cell)  Paths to raw_*.xlsx files
%                         {default: canonical NetaF files}
%       'outFile'       - (char)  Path for sponCa_raw.xlsx
%                         {default: NetaF folder \ sponCa_raw.xlsx}
%       'cellsExclPath' - (char)  Path to CellsExcluded.xlsx (used to set
%                         per-cell excluded flag). If '', the flag is set
%                         to false for all cells.   {<NetaF>\CellsExcluded.xlsx}
%       'fs'            - (num)   Sampling rate (Hz)                {3}
%       'pCyto'         - (num)   Percentile for cyto F0            {20}
%       'pMito'         - (num)   Percentile for mito F0            {5}
%       'corrTol'       - (num)   Min repro-vs-sheet correlation    {0.90}
%       'msWinFrac'     - (num)   msbackadj WindowSize / nT         {0.1}
%       'verbose'       - (log)   Print progress                    {true}
%       'diagDir'       - (char)  Folder for diagnostic PNGs + report
%                         {default: <thisFolder>\_diag\spontCa_extractRaw}
%
%   OUTPUT:
%       outArr - cell array written to the xlsx
%       report - struct with .perFile, .corrSummary, .failCells, .excluded
%                .ctrl_total / .ko_total / .ctrl_kept / .ko_kept, etc.
%
%   See also: SPONTCA_LOAD, CA_EVENTDETECTOR

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addParameter(p, 'srcFiles',      {},   @(x) iscell(x) || isstring(x));
addParameter(p, 'outFile',       '',   @(x) ischar(x) || isstring(x));
addParameter(p, 'cellsExclPath', '',   @(x) ischar(x) || isstring(x));
addParameter(p, 'fs',            3,    @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'pCyto',         20,   @isnumeric);
addParameter(p, 'pMito',         5,    @isnumeric);
addParameter(p, 'corrTol',       0.90, @isnumeric);
addParameter(p, 'msWinFrac',     0.1,  @isnumeric);
addParameter(p, 'verbose',       true, @islogical);
addParameter(p, 'diagDir',       '',   @(x) ischar(x) || isstring(x));
parse(p, varargin{:});
P = p.Results;

baseDataDir = fullfile('D:\OneDrive - Tel-Aviv University', ...
    'PhD\Slutsky\Manuscripts\MCU\Results\Data_NetaF');

if isempty(P.srcFiles)
    srcFiles = { ...
        autoSrcFile(fullfile(baseDataDir, '20250505+06')), ...
        autoSrcFile(fullfile(baseDataDir, '20250907+09')) };
else
    srcFiles = cellstr(P.srcFiles);
end

if isempty(P.outFile)
    outFile = fullfile(baseDataDir, 'sponCa_raw.xlsx');
else
    outFile = char(P.outFile);
end

if isempty(P.cellsExclPath)
    cellsExclPath = fullfile(baseDataDir, 'CellsExcluded.xlsx');
else
    cellsExclPath = char(P.cellsExclPath);
end

if isempty(P.diagDir)
    diagDir = fullfile(fileparts(mfilename('fullpath')), ...
        '_diag', 'spontCa_extractRaw');
else
    diagDir = char(P.diagDir);
end
if ~isfolder(diagDir), mkdir(diagDir); end

for iF = 1:numel(srcFiles)
    assert(isfile(srcFiles{iF}), 'spontCa_extractRaw:fileNotFound', ...
        'Source file not found: %s', srcFiles{iF});
end

if P.verbose
    fprintf('[spontCa_extractRaw] %d source file(s)\n', numel(srcFiles));
end

hasMsBackAdj = exist('msbackadj', 'file') == 2;
if ~hasMsBackAdj && P.verbose
    fprintf('  msbackadj unavailable - falling back to movmedian baseline\n');
end


%% ========================================================================
%  PARSE EXCLUSION FILE
%  ========================================================================

if ~isempty(cellsExclPath) && isfile(cellsExclPath)
    exclTbl = parseCellsExcluded(cellsExclPath);
    if P.verbose
        fprintf('  parsed %s: %d entries\n', cellsExclPath, height(exclTbl));
    end
else
    exclTbl = table('Size', [0, 4], ...
        'VariableTypes', {'string', 'string', 'double', 'logical'}, ...
        'VariableNames', {'experiment', 'genotype', 'position', 'excluded'});
    if P.verbose
        fprintf('  no CellsExcluded.xlsx - all cells marked excluded=false\n');
    end
end


%% ========================================================================
%  SHEET CONFIG
%  ========================================================================

genoCfg = struct( ...
    'gen',     {'Control',         'MCU-KO'}, ...
    'shCytoR', {'MCU ctrl - cyto', 'MCU KO - cyto'}, ...
    'shCytoD', {'MCU ctrl - cyto dff', 'MCU KO - cyto dff'}, ...
    'shMitoR', {'MCU ctrl - mito', 'MCU KO - mito'}, ...
    'shMitoD', {'MCU ctrl - mito dff', 'MCU KO - mito dff'}, ...
    'shCytoX', {'cyto Ca influx - ctrl', 'cyto Ca influx - KO'}, ...
    'shMitoX', {'mito Ca influx - ctrl', 'mito ca influx - KO'});

% Per-file experiment label (used as key into CellsExcluded)
expLabel = {'I', 'II'};
assert(numel(srcFiles) <= numel(expLabel), 'Need to extend expLabel');


%% ========================================================================
%  COLLECT CELLS (all with valid raw in both compartments)
%  ========================================================================

cells = struct('genotype', {}, 'srcFile', {}, 'experiment', {}, ...
               'positionInExp', {}, 'origCol', {}, 'excluded', {}, ...
               'cyto', {}, 'mito', {}, ...
               'corrCyto', {}, 'corrMito', {}, ...
               'cytoDffRepro', {}, 'mitoDffRepro', {}, ...
               'cytoDffSheet', {}, 'mitoDffSheet', {});

perFile = struct('file', {}, 'experiment', {}, 'genotype', {}, ...
                 'nTotal', {}, 'nDffOk', {}, 'nWithEvts', {}, ...
                 'nExpExcluded', {}, 'nExpKept', {});

allFails = struct('file', {}, 'genotype', {}, 'cellID', {}, ...
                  'compartment', {}, 'corr', {});

for iFile = 1:numel(srcFiles)
    fpath = srcFiles{iFile};
    fname = srcFileTag(fpath);
    expS  = expLabel{iFile};
    if P.verbose, fprintf('[%s = Exp %s]\n', fname, expS); end

    for iG = 1:numel(genoCfg)
        gen = genoCfg(iG).gen;
        sh  = genoCfg(iG);

        rawC = readmatrix(fpath, 'Sheet', sh.shCytoR, 'NumHeaderLines', 1);
        dffC = readmatrix(fpath, 'Sheet', sh.shCytoD, 'NumHeaderLines', 1);
        rawM = readmatrix(fpath, 'Sheet', sh.shMitoR, 'NumHeaderLines', 1);
        dffM = readmatrix(fpath, 'Sheet', sh.shMitoD, 'NumHeaderLines', 1);

        nT = min([size(rawC,1), size(dffC,1), size(rawM,1), size(dffM,1)]);
        rawC = rawC(1:nT, :); dffC = dffC(1:nT, :);
        rawM = rawM(1:nT, :); dffM = dffM(1:nT, :);

        assert(size(rawC,2) == size(rawM,2), ...
            'spontCa_extractRaw:colMismatch', ...
            'Cyto/Mito column count mismatch in %s/%s', fname, gen);
        nColsRead = size(rawC, 2);

        % Real cells: BOTH cyto and mito raw columns are non-empty.
        % Gap columns (separators between recording blocks) fail this test.
        isCell = false(1, nColsRead);
        for k = 1:nColsRead
            isCell(k) = ~isColEmpty(rawC(:, k)) && ~isColEmpty(rawM(:, k));
        end
        cellCols = find(isCell);

        % Influx-based event presence (informational; used for report)
        hasCytoEvt = readEventMask(fpath, sh.shCytoX, cellCols);
        hasMitoEvt = readEventMask(fpath, sh.shMitoX, cellCols);

        % Per-cell exclusion flag from CellsExcluded.xlsx (matched on
        % experiment x genotype x positionInExp)
        exclFlags = lookupExclusion(exclTbl, expS, gen, numel(cellCols));

        nDffOk = 0;
        for posExp = 1:numel(cellCols)
            k = cellCols(posExp);
            cytoRaw = rawC(:, k);
            mitoRaw = rawM(:, k);
            cytoDff = dffC(:, k);
            mitoDff = dffM(:, k);

            dffOk = ~isColEmpty(cytoDff) && ~isColEmpty(mitoDff);
            if dffOk, nDffOk = nDffOk + 1; end

            cytoRepr = boazDff(cytoRaw, P.pCyto, P.msWinFrac, hasMsBackAdj);
            mitoRepr = boazDff(mitoRaw, P.pMito, P.msWinFrac, hasMsBackAdj);

            maxLag = round(0.50 * nT);
            if dffOk
                cCyto = corrBestLag(cytoRepr, cytoDff, maxLag);
                cMito = corrBestLag(mitoRepr, mitoDff, maxLag);
            else
                cCyto = NaN; cMito = NaN;
            end

            entry = struct();
            entry.genotype      = gen;
            entry.srcFile       = fname;
            entry.experiment    = expS;
            entry.positionInExp = posExp;
            entry.origCol       = k;
            entry.excluded      = exclFlags(posExp);
            entry.cyto          = cytoRaw;
            entry.mito          = mitoRaw;
            entry.corrCyto      = cCyto;
            entry.corrMito      = cMito;
            entry.cytoDffRepro  = cytoRepr;
            entry.mitoDffRepro  = mitoRepr;
            entry.cytoDffSheet  = cytoDff;
            entry.mitoDffSheet  = mitoDff;
            cells(end+1) = entry; %#ok<AGROW>
        end

        pf = struct();
        pf.file         = fname;
        pf.experiment   = expS;
        pf.genotype     = gen;
        pf.nTotal       = numel(cellCols);
        pf.nDffOk       = nDffOk;
        pf.nWithEvts    = sum(hasCytoEvt & hasMitoEvt);
        pf.nExpExcluded = sum(exclFlags);
        pf.nExpKept     = sum(~exclFlags);
        perFile(end+1) = pf; %#ok<AGROW>

        if P.verbose
            fprintf('  %-8s: %d cells | dff-ok=%d | with-events=%d | exp-kept=%d / exp-excl=%d\n', ...
                gen, pf.nTotal, pf.nDffOk, pf.nWithEvts, pf.nExpKept, pf.nExpExcluded);
        end
    end
end


%% ========================================================================
%  ORDER CELLS AND ASSIGN GLOBAL IDs
%  ========================================================================
%
%   Order: Control first, MCU-KO second; within each genotype, by
%   (file index, position-in-experiment). Global ID = sequential rank within
%   each genotype across files - so the k-th ctrl cell in file 1 has ID k,
%   and file-2 ctrl cells continue from N_file1_ctrl + 1.

genOrder  = arrayfun(@(c) double(strcmp(c.genotype, 'MCU-KO')), cells);
expOrder  = arrayfun(@(c) double(strcmp(c.experiment, 'II')), cells);
posOrder  = [cells.positionInExp];
[~, sortIdx] = sortrows([genOrder(:), expOrder(:), posOrder(:)]);
cells = cells(sortIdx);

idCtrl = 0; idKo = 0;
for k = 1:numel(cells)
    if strcmp(cells(k).genotype, 'Control')
        idCtrl = idCtrl + 1;
        cells(k).cellID = idCtrl;
    else
        idKo = idKo + 1;
        cells(k).cellID = idKo;
    end
end


%% ========================================================================
%  COLLECT CORR FAILURES
%  ========================================================================

for k = 1:numel(cells)
    if cells(k).corrCyto < P.corrTol
        allFails(end+1) = struct('file', cells(k).srcFile, ...
            'genotype', cells(k).genotype, 'cellID', cells(k).cellID, ...
            'compartment', 'Cyto', 'corr', cells(k).corrCyto); %#ok<AGROW>
    end
    if cells(k).corrMito < P.corrTol
        allFails(end+1) = struct('file', cells(k).srcFile, ...
            'genotype', cells(k).genotype, 'cellID', cells(k).cellID, ...
            'compartment', 'Mito', 'corr', cells(k).corrMito); %#ok<AGROW>
    end
end


%% ========================================================================
%  BUILD OUTPUT CELL ARRAY
%  ========================================================================
%
%   Row 1 = genotype, Row 2 = ID, Row 3 = excluded, Row 4 = Cyto/Mito,
%   Row 5+ = time + data (col 1 = time).

nCellsOut = numel(cells);
nT        = min(arrayfun(@(c) numel(c.cyto), cells));
nCols     = 1 + 2 * nCellsOut;
outArr    = cell(4 + nT, nCols);

t = (0:nT-1)' / P.fs;
outArr{4, 1} = 'Time';
for iT = 1:nT
    outArr{4 + iT, 1} = t(iT);
end

for k = 1:nCellsOut
    ic = 2*k;
    im = 2*k + 1;

    outArr{1, ic} = cells(k).genotype;
    outArr{1, im} = cells(k).genotype;
    outArr{2, ic} = cells(k).cellID;
    % outArr{2, im} stays empty (matches SpontCa.xlsx)
    outArr{3, ic} = cells(k).excluded;
    outArr{3, im} = cells(k).excluded;
    outArr{4, ic} = 'Cyto';
    outArr{4, im} = 'Mito';

    cytoTrim = cells(k).cyto(1:nT);
    mitoTrim = cells(k).mito(1:nT);
    for iT = 1:nT
        outArr{4 + iT, ic} = cytoTrim(iT);
        outArr{4 + iT, im} = mitoTrim(iT);
    end
end


%% ========================================================================
%  WRITE XLSX
%  ========================================================================

if isfile(outFile), delete(outFile); end
writecell(outArr, outFile, 'Sheet', 'Sheet1');

if P.verbose
    fprintf('  wrote %s (%dx%d)\n', outFile, size(outArr,1), size(outArr,2));
end


%% ========================================================================
%  REPORT + DIAGNOSTIC PNGs
%  ========================================================================

corrSummary = buildCorrSummary(cells);
failTbl     = struct2tableSafe(allFails);

ctrlTot = sum(strcmp({cells.genotype}, 'Control'));
koTot   = sum(strcmp({cells.genotype}, 'MCU-KO'));
ctrlKept = sum(strcmp({cells.genotype}, 'Control') & ~[cells.excluded]);
koKept   = sum(strcmp({cells.genotype}, 'MCU-KO')  & ~[cells.excluded]);

reportPath = fullfile(diagDir, 'report.txt');
writeReport(reportPath, srcFiles, cells, perFile, corrSummary, failTbl, ...
    ctrlTot, koTot, ctrlKept, koKept, P, outFile);

savePngs(cells, diagDir, P.corrTol);

report = struct();
report.cells         = cells;
report.ctrl_total    = ctrlTot;
report.ko_total      = koTot;
report.ctrl_kept     = ctrlKept;
report.ko_kept       = koKept;
report.perFile       = perFile;
report.corrSummary   = corrSummary;
report.failCells     = failTbl;
report.outFile       = outFile;
report.diagDir       = diagDir;
report.reportPath    = reportPath;

if P.verbose
    fprintf('\nTotals:\n');
    fprintf('  Ctrl: %d cells, %d kept by experimenter (target 53)\n', ctrlTot, ctrlKept);
    fprintf('  KO  : %d cells, %d kept by experimenter (target 48)\n', koTot, koKept);
    fprintf('  report  : %s\n', reportPath);
    fprintf('  diagDir : %s\n', diagDir);
end

end     % main function


%% ========================================================================
%  HELPERS
%  ========================================================================

function tf = isColEmpty(col)
%ISCOLEMPTY  True if >= 95% of column is NaN.
tf = mean(isnan(col)) >= 0.95;
end


function hasEvt = readEventMask(fpath, sheetName, cellCols)
%READEVENTMASK  Per-cell flag: amplitude > 0 in the Ca-influx sheet.
infl = readmatrix(fpath, 'Sheet', sheetName, 'NumHeaderLines', 1);
nRows = size(infl, 1);
amp = nan(nRows, 1);
if size(infl, 2) >= 1, amp = infl(:, 1); end
nCells = numel(cellCols);
if nRows < nCells
    amp(nCells, 1) = NaN;
elseif nRows > nCells
    amp = amp(1:nCells);
end
hasEvt = (~isnan(amp) & amp > 0).';
end


function flags = lookupExclusion(exclTbl, expS, gen, nExpected)
%LOOKUPEXCLUSION  Vector of excluded flags for each (expS,gen) cell.
%   Returns false() if no entry exists (no exclusion file).
flags = false(1, nExpected);
if isempty(exclTbl) || height(exclTbl) == 0, return; end
mask = strcmp(string(exclTbl.experiment), expS) ...
     & strcmp(string(exclTbl.genotype),   gen);
sub = exclTbl(mask, :);
if isempty(sub), return; end
sub = sortrows(sub, 'position');
nFound = height(sub);
if nFound ~= nExpected
    warning('spontCa_extractRaw:exclMismatch', ...
        'CellsExcluded has %d entries for Exp %s / %s but raw has %d cells', ...
        nFound, expS, gen, nExpected);
end
for r = 1:min(nFound, nExpected)
    flags(sub.position(r)) = sub.excluded(r);
end
end


function dff = boazDff(rawCol, pPct, winFrac, hasMsBackAdj)
%BOAZDFF  Reproduce Ca_EventDetector.m dF/F (lines 28-64).
nT = numel(rawCol);
if any(isnan(rawCol))
    rawCol = fillmissing(rawCol, 'linear', 'EndValues', 'nearest');
end
winSize = max(round(winFrac * nT), 5);
if hasMsBackAdj
    try
        b = msbackadj((1:nT)', rawCol, 'WindowSize', winSize, 'ShowPlot', 0);
    catch
        b = rawCol - movmedian(rawCol, winSize);
    end
else
    b = rawCol - movmedian(rawCol, winSize);
end
b(b < 0) = 0;
F0 = prctile(rawCol, pPct);
if F0 == 0 || isnan(F0), F0 = eps; end
dff = b / F0;
end


function r = corrNaN(x, y)
%CORRNAN  Pearson correlation ignoring NaNs in either input.
x = x(:); y = y(:);
nL = min(numel(x), numel(y));
x = x(1:nL); y = y(1:nL);
m = ~(isnan(x) | isnan(y));
if sum(m) < 3, r = NaN; return; end
xs = x(m) - mean(x(m));
ys = y(m) - mean(y(m));
d = sqrt(sum(xs.^2) * sum(ys.^2));
if d == 0, r = NaN; else, r = sum(xs .* ys) / d; end
end


function r = corrBestLag(x, y, maxLag)
%CORRBESTLAG  Best Pearson correlation across integer lags in [-maxLag, maxLag].
x = x(:); y = y(:);
r = corrNaN(x, y);
nx = numel(x);
maxLag = min(maxLag, nx - 10);
if maxLag <= 0, return; end
step = max(1, round(maxLag / 60));
for lag = step:step:maxLag
    rPos = corrNaN(x(1+lag:end), y(1:end-lag));
    rNeg = corrNaN(x(1:end-lag), y(1+lag:end));
    r = max([r, rPos, rNeg]);
end
end


function fp = autoSrcFile(folder)
%AUTOSRCFILE  Pick the raw_*.xlsx in folder (excluding Excel lock files).
d = dir(fullfile(folder, 'raw_*.xlsx'));
d = d(~startsWith({d.name}, '~$'));
assert(~isempty(d), 'spontCa_extractRaw:noSrcFile', ...
    'No raw_*.xlsx in %s', folder);
fp = fullfile(folder, d(1).name);
end


function tag = srcFileTag(fpath)
%SRCFILETAG  Compact identifier: <parentFolder>/<base>.
[parent, base] = fileparts(fpath);
[~, parentName] = fileparts(parent);
if isempty(parentName), tag = base; else, tag = [parentName '/' base]; end
end


function exclTbl = parseCellsExcluded(path)
%PARSECELLSEXCLUDED  Parse CellsExcluded.xlsx into a tidy table.
%   Layout: col A has 'Exp I' / 'Exp II' labels; cols B/C hold values for
%   Ctrl / KO. A trailing asterisk in a value means experimenter-excluded.
raw = readcell(path);
nR = size(raw, 1);
expIStart = 0; expIIStart = 0;
for r = 1:nR
    s = cellToText(raw{r, 1});
    if contains(s, 'Exp') && contains(s, 'II')
        expIIStart = r;
    elseif contains(s, 'Exp') && expIStart == 0
        expIStart = r;
    end
end
assert(expIStart > 0 && expIIStart > 0, 'parseCellsExcluded:noMarkers', ...
    'Could not find "Exp I" and "Exp II" markers in %s', path);

genCols  = [2, 3];
genNames = {'Control', 'MCU-KO'};
expRanges = {expIStart:expIIStart-1, expIIStart:nR};
expNames  = {'I', 'II'};

experiment = strings(0, 1);
genotype   = strings(0, 1);
position   = zeros(0, 1);
excluded   = false(0, 1);

for iE = 1:2
    rng = expRanges{iE};
    for iG = 1:2
        pos = 0;
        for r = rng
            if r > size(raw, 2)*0 + size(raw, 1)
                break;
            end
            flag = parseExclFlag(raw{r, genCols(iG)});
            if isnan(flag), continue; end
            pos = pos + 1;
            experiment(end+1, 1) = expNames{iE}; %#ok<AGROW>
            genotype(end+1,   1) = genNames{iG}; %#ok<AGROW>
            position(end+1,   1) = pos;          %#ok<AGROW>
            excluded(end+1,   1) = logical(flag); %#ok<AGROW>
        end
    end
end

exclTbl = table(experiment, genotype, position, excluded);
end


function flag = parseExclFlag(v)
%PARSEEXCLFLAG  NaN if empty; false if numeric; true if string ending in '*'.
if isnumeric(v)
    if isempty(v) || (isscalar(v) && isnan(v))
        flag = NaN;
    else
        flag = false;
    end
elseif ischar(v) || isstring(v)
    s = string(v);
    if strlength(s) == 0
        flag = NaN;
    else
        flag = endsWith(s, "*");
    end
else
    % missing or unsupported type
    flag = NaN;
end
end


function s = cellToText(v)
%CELLTOTEXT  Convert a cell value to char for text searching (empty if not text).
if ischar(v)
    s = v;
elseif isstring(v)
    s = char(v);
else
    s = '';
end
end


function S = buildCorrSummary(cells)
%BUILDCORRSUMMARY  Per genotype x compartment: min / median / max corr
%   (over cells where dff exists, i.e. corr is not NaN).
genL  = {'Control', 'MCU-KO'};
comp  = {'Cyto', 'Mito'};
rows  = {};
for ig = 1:numel(genL)
    cMask = strcmp({cells.genotype}, genL{ig});
    if ~any(cMask), continue; end
    cytoCorrs = [cells(cMask).corrCyto];
    mitoCorrs = [cells(cMask).corrMito];
    cytoOk = cytoCorrs(~isnan(cytoCorrs));
    mitoOk = mitoCorrs(~isnan(mitoCorrs));
    rows(end+1, :) = {genL{ig}, comp{1}, ...
        safeMin(cytoOk), safeMed(cytoOk), safeMax(cytoOk), numel(cytoOk)}; %#ok<AGROW>
    rows(end+1, :) = {genL{ig}, comp{2}, ...
        safeMin(mitoOk), safeMed(mitoOk), safeMax(mitoOk), numel(mitoOk)}; %#ok<AGROW>
end
S = cell2table(rows, 'VariableNames', ...
    {'genotype', 'compartment', 'minCorr', 'medCorr', 'maxCorr', 'n'});
end


function v = safeMin(x), if isempty(x), v = NaN; else, v = min(x); end, end
function v = safeMed(x), if isempty(x), v = NaN; else, v = median(x); end, end
function v = safeMax(x), if isempty(x), v = NaN; else, v = max(x); end, end


function T = struct2tableSafe(S)
if isempty(S), T = table(); else, T = struct2table(S); end
end


function writeReport(reportPath, srcFiles, cells, perFile, corrSummary, ...
                     failTbl, ctrlTot, koTot, ctrlKept, koKept, P, outFile)
fid = fopen(reportPath, 'w');
cu = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'spontCa_extractRaw report\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS')); %#ok<TNOW1,DATST>

fprintf(fid, 'Source files:\n');
for i = 1:numel(srcFiles), fprintf(fid, '  %s\n', srcFiles{i}); end
fprintf(fid, '\n');

fprintf(fid, 'Params: fs=%g, pCyto=%g, pMito=%g, msWinFrac=%g, corrTol=%g\n\n', ...
    P.fs, P.pCyto, P.pMito, P.msWinFrac, P.corrTol);

fprintf(fid, '--- Per-file counts ---\n');
fprintf(fid, '%-26s %-3s %-8s %6s %6s %6s %6s %6s\n', ...
    'file', 'exp', 'genotype', 'nTot', 'dffOk', 'wEvts', 'kept', 'excl');
for i = 1:numel(perFile)
    pf = perFile(i);
    fprintf(fid, '%-26s %-3s %-8s %6d %6d %6d %6d %6d\n', ...
        pf.file, pf.experiment, pf.genotype, pf.nTotal, pf.nDffOk, ...
        pf.nWithEvts, pf.nExpKept, pf.nExpExcluded);
end
fprintf(fid, '  nTot   = cells with valid raw in both compartments\n');
fprintf(fid, '  dffOk  = cells with non-empty dff in both compartments\n');
fprintf(fid, '  wEvts  = cells with detected events in both compartments\n');
fprintf(fid, '  kept   = NOT excluded by experimenter (per CellsExcluded.xlsx)\n');
fprintf(fid, '  excl   = excluded by experimenter (asterisk in CellsExcluded)\n\n');

fprintf(fid, '--- Totals ---\n');
fprintf(fid, '  Ctrl: %d cells, %d kept (target 53), %d excluded\n', ...
    ctrlTot, ctrlKept, ctrlTot - ctrlKept);
fprintf(fid, '  KO  : %d cells, %d kept (target 48), %d excluded\n\n', ...
    koTot, koKept, koTot - koKept);

fprintf(fid, '--- Cross-check: dff-empty vs experimenter-excluded ---\n');
fprintf(fid, '  (counts of cells in each combination of the two flags)\n');
for ig = 1:2
    genL = {'Control', 'MCU-KO'};
    gen = genL{ig};
    mask = strcmp({cells.genotype}, gen);
    sub = cells(mask);
    dffEmpty = arrayfun(@(c) isnan(c.corrCyto) || isnan(c.corrMito), sub);
    excl     = [sub.excluded];
    n11 = sum( dffEmpty &  excl);
    n10 = sum( dffEmpty & ~excl);
    n01 = sum(~dffEmpty &  excl);
    n00 = sum(~dffEmpty & ~excl);
    fprintf(fid, '  %-8s   dff-empty/excl=%d   dff-empty/kept=%d   dff-ok/excl=%d   dff-ok/kept=%d\n', ...
        gen, n11, n10, n01, n00);
end
fprintf(fid, '\n');

fprintf(fid, '--- Correlation summary (Boaz-dff repro vs sheet) ---\n');
if ~isempty(corrSummary)
    for i = 1:height(corrSummary)
        fprintf(fid, '  %-8s %-5s  n=%3d  min=%.3f  med=%.3f  max=%.3f\n', ...
            corrSummary.genotype{i}, corrSummary.compartment{i}, ...
            corrSummary.n(i), corrSummary.minCorr(i), ...
            corrSummary.medCorr(i), corrSummary.maxCorr(i));
    end
end
fprintf(fid, '\n');

fprintf(fid, '--- Sub-tolerance cells (corr < %.2f) ---\n', P.corrTol);
if isempty(failTbl) || height(failTbl) == 0
    fprintf(fid, '  none\n');
else
    for i = 1:height(failTbl)
        fprintf(fid, '  %-8s %s%d  (%s)  corr=%.3f\n', ...
            failTbl.genotype{i}, ...
            shortGen(failTbl.genotype{i}), ...
            failTbl.cellID(i), failTbl.compartment{i}, failTbl.corr(i));
    end
end
fprintf(fid, '\n');

fprintf(fid, 'Output: %s\n', outFile);
end


function s = shortGen(g)
if strcmpi(g, 'Control'), s = 'Ctrl_'; else, s = 'KO_'; end
end


function savePngs(cells, diagDir, corrTol)
%SAVEPNGS  Save 3 worst-correlation overlays per genotype x compartment
%   (only over cells where dff exists).
genL = {'Control', 'MCU-KO'};
comp = {'Cyto', 'Mito'};
for ig = 1:numel(genL)
    mask = strcmp({cells.genotype}, genL{ig});
    if ~any(mask), continue; end
    idx = find(mask);
    for ic = 1:numel(comp)
        if strcmp(comp{ic}, 'Cyto')
            corrs = [cells(idx).corrCyto];
        else
            corrs = [cells(idx).corrMito];
        end
        ok = ~isnan(corrs);
        if ~any(ok), continue; end
        idxOk = idx(ok);
        corrsOk = corrs(ok);
        [sortedC, sIdx] = sort(corrsOk);
        nPlot = min(3, numel(sIdx));
        for j = 1:nPlot
            cIdx = idxOk(sIdx(j));
            c = cells(cIdx);
            f = figure('Visible', 'off', 'Position', [100, 100, 1200, 700]);
            tl = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
            title(tl, sprintf('%s %s ID=%d (%s, pos %d) corr=%.3f%s', ...
                c.genotype, comp{ic}, c.cellID, c.srcFile, c.positionInExp, sortedC(j), ...
                tern(sortedC(j) < corrTol, '  [BELOW TOL]', '')), ...
                'Interpreter', 'none');

            if strcmp(comp{ic}, 'Cyto')
                raw  = c.cyto;
                repr = c.cytoDffRepro;
                sht  = c.cytoDffSheet;
            else
                raw  = c.mito;
                repr = c.mitoDffRepro;
                sht  = c.mitoDffSheet;
            end

            nexttile; plot(raw, 'k'); ylabel('raw F'); grid on;
            nexttile; plot(repr, 'b'); ylabel('dF/F (repro)'); grid on;
            nexttile; plot(sht, 'r'); ylabel('dF/F (sheet)'); grid on; xlabel('sample');

            fname = sprintf('worst_%s_%s_%02d_ID%d.png', ...
                genL{ig}, comp{ic}, j, c.cellID);
            fname = strrep(fname, '-', '_');
            exportgraphics(f, fullfile(diagDir, fname), 'Resolution', 120);
            close(f);
        end
    end
end
end


function s = tern(cond, a, b)
if cond, s = a; else, s = b; end
end
