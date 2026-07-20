function [tbl, fs] = spontCa_loadXls(varargin)
% SPONTCA_LOADXLS Reads spontCa.xlsx into a long-format unit-level table.
%
%   [tbl, fs] = SPONTCA_LOADXLS(...) returns one row per (cell x compartment).
%   The xlsx has two sheets: 'f' (raw fluorescence) and 'dff' (Boaz dF/F).
%   Picks which sheet via flgRaw and optionally drops experimenter-excluded
%   cells via flgExclude. Cell numbering is preserved across exclusions:
%   if cell #4 is excluded, cell #5 still has unitID=5.
%
%   OUTPUT SCHEMA:
%       genotype     (categorical) 'Control' / 'MCU-KO'
%       sbjID        (categorical) 'Ctrl_05', 'KO_27', ...
%       unitID       (categorical) cell number from row 2 of the xlsx
%       compartment  (categorical) 'Cyto' / 'Mito'
%       excluded     (logical)     experimenter's flag from row 3
%       trace        (1 x nT)      dF/F (see flgRaw for source)
%
%   OPTIONAL (Name-Value):
%       'xlsPath'    - (char) Path to spontCa.xlsx
%                      {default: NetaF folder \ spontCa.xlsx}
%       'flgRaw'     - (log)  true: read sheet 'f' (raw F) and compute
%                             dF/F here via a rolling 20th-percentile
%                             baseline over a 30 s window;
%                             false: read sheet 'dff' as-is (Boaz dF/F).
%                             Either way the returned trace is dF/F.  {false}
%       'flgExclude' - (log)  true: drop rows where excluded=true;
%                             unitIDs are preserved (gaps remain)    {false}
%       'flgCorrF0'  - (log)  true: apply per-cell F0 trace correction
%                             after the dF/F conversion. Each cell's
%                             trace is scaled by F0_ref/F0_cell, with
%                             F0_ref = median F0 in its compartment, so
%                             baseline-brightness differences between
%                             cells (and genotypes) don't bias dF/F
%                             amplitude. Requires flgRaw=true.       {false}
%       'verbose'    - (log)  Print progress                         {true}
%
%   FILE FORMAT (rows 1-indexed, both sheets):
%       row 1   : genotype, repeated over each Cyto/Mito pair
%       row 2   : numeric cell ID, above Cyto col only
%       row 3   : excluded flag (TRUE/FALSE), repeated over the pair
%       row 4   : 'Time' (col 1) / 'Cyto' / 'Mito'
%       row 5+  : data. col 1 = time (s); col 2+ = signal columns
%
%   See also: SPONTCA_DETECT, SPONTCA_GUI

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addParameter(p, 'xlsPath',    '',    @(x) ischar(x) || isstring(x));
addParameter(p, 'flgRaw',     false, @islogical);
addParameter(p, 'flgExclude', false, @islogical);
addParameter(p, 'flgCorrF0',  false, @islogical);
addParameter(p, 'verbose',    true,  @islogical);
parse(p, varargin{:});

xlsPath    = char(p.Results.xlsPath);
flgRaw     = p.Results.flgRaw;
flgExclude = p.Results.flgExclude;
flgCorrF0  = p.Results.flgCorrF0;
verbose    = p.Results.verbose;

if isempty(xlsPath)
    xlsPath = fullfile('D:\OneDrive - Tel-Aviv University', ...
        'PhD\Slutsky\Manuscripts\MCU\Results\Data_NetaF', 'spontCa.xlsx');
end
assert(isfile(xlsPath), 'spontCa_loadXls:fileNotFound', ...
    'Excel file not found: %s', xlsPath);

if flgRaw, sheetName = 'f'; else, sheetName = 'dff'; end

if verbose
    fprintf('[spontCa_loadXls] %s  sheet=%s  flgExclude=%d  flgCorrF0=%d\n', ...
        xlsPath, sheetName, flgExclude, flgCorrF0);
end


%% ========================================================================
%  READ EXCEL
%  ========================================================================

raw = readcell(xlsPath, 'Sheet', sheetName);
[nRowsRaw, nColsRaw] = size(raw);

nSig    = nColsRaw - 1;
hdrGen  = strings(1, nSig);
hdrCid  = nan(1, nSig);
hdrExcl = false(1, nSig);
hdrSig  = strings(1, nSig);
for iC = 1:nSig
    val = raw{1, iC + 1};
    if ~ismissing(val), hdrGen(iC) = string(val); end
    val = raw{2, iC + 1};
    if isnumeric(val) && ~isnan(val), hdrCid(iC) = val; end
    val = raw{3, iC + 1};
    hdrExcl(iC) = coerceLogical(val);
    val = raw{4, iC + 1};
    if ~ismissing(val), hdrSig(iC) = string(val); end
end

dataMat = nan(nRowsRaw - 4, nColsRaw);
for iR = 5:nRowsRaw
    for iC = 1:nColsRaw
        val = raw{iR, iC};
        if isnumeric(val) && ~isempty(val) && ~ismissing(val)
            dataMat(iR - 4, iC) = val;
        end
    end
end
dataMat = dataMat(~isnan(dataMat(:, 1)), :);    % drop ragged tail

t   = dataMat(:, 1)';
sig = dataMat(:, 2:end);
fs  = 1 / median(diff(t), 'omitnan');
nT  = length(t);


%% ========================================================================
%  BUILD LONG-FORMAT TABLE
%  ========================================================================

% Pair each Cyto column with the Mito column immediately to its right.
cytoIdx = find(strcmpi(hdrSig, 'Cyto'));
nCells  = length(cytoIdx);
n       = 2 * nCells;

genotype    = strings(n, 1);
sbjID       = strings(n, 1);
compartment = strings(n, 1);
unitID      = nan(n, 1);
excluded    = false(n, 1);
traces      = nan(n, nT);

iRow = 0;
for iC = 1:nCells
    cC   = cytoIdx(iC);
    cM   = cC + 1;
    gen  = char(hdrGen(cC));
    cnum = hdrCid(cC);
    excl = hdrExcl(cC);

    if strcmpi(gen, 'Control')
        prefix = 'Ctrl';
    elseif strcmpi(gen, 'MCU-KO')
        prefix = 'KO';
    else
        prefix = regexprep(gen, '\W', '');
    end
    sName = sprintf('%s_%02d', prefix, cnum);

    iRow = iRow + 1;
    genotype(iRow)    = gen;
    sbjID(iRow)       = sName;
    compartment(iRow) = 'Cyto';
    unitID(iRow)      = cnum;
    excluded(iRow)    = excl;
    traces(iRow, :)   = sig(:, cC)';

    iRow = iRow + 1;
    genotype(iRow)    = gen;
    sbjID(iRow)       = sName;
    compartment(iRow) = 'Mito';
    unitID(iRow)      = cnum;
    excluded(iRow)    = excl;
    if cM <= nSig && strcmpi(hdrSig(cM), 'Mito') && hdrGen(cM) == hdrGen(cC)
        traces(iRow, :) = sig(:, cM)';
    else
        warning('spontCa_loadXls:noMito', ...
            'No paired Mito column for cell %d (%s).', cnum, gen);
    end
end

% cfg.lbl.grp also lists the in vivo viral cohort, which imaging does not
% have. removecats drops that empty level - fitlme rejects a predictor
% holding a level with no observations
cfg = mcu_cfg;
tbl = table(...
    removecats(categorical(genotype, cfg.lbl.grp)), ...
    categorical(sbjID), ...
    categorical(unitID), ...
    categorical(compartment, {'Cyto', 'Mito'}), ...
    excluded, ...
    traces, ...
    'VariableNames', ...
    {'genotype', 'sbjID', 'unitID', 'compartment', 'excluded', 'trace'});

% If we loaded raw F, convert to dF/F here via a rolling 20th-percentile
% baseline (deterministic recipe; no per-cell tuning). When flgRaw is
% false the trace is already dF/F from the 'dff' sheet.
if flgRaw
    % Per-cell median of raw F before the dF/F conversion. Used by the
    % F0 diagnostic in detectWrapper and by the per-cell F0 trace
    % correction below.
    tbl.F0 = median(tbl.trace, 2, 'omitnan');

    bslWin   = 30;                  % baseline window (s)
    bslQuant = 20;                  % percentile for baseline
    winSamps = round(bslWin * fs);
    halfWin  = floor(winSamps / 2);
    for iR = 1:height(tbl)
        f  = tbl.trace(iR, :);
        f0 = zeros(1, nT);
        for iS = 1:nT
            i0 = max(1, iS - halfWin);
            i1 = min(nT, iS + halfWin);
            f0(iS) = prctile(f(i0:i1), bslQuant);
        end
        tbl.trace(iR, :) = (f - f0) ./ max(f0, eps);
    end

    % Per-cell F0 trace correction (opt-in via flgCorrF0). Cells with
    % above-median F0 in their compartment get their dF/F values scaled
    % down (and vice-versa) so baseline brightness no longer biases
    % dF/F amplitude across cells. Reference is the per-compartment
    % median F0 across both genotypes (neutral - doesn't privilege one
    % group).
    %
    % Caveat: this assumes the F0 difference between genotypes is
    % technical (indicator concentration / imaging conditions). If KO
    % has elevated resting cyto Ca2+ that drives both F0 and dF/F amp
    % (biological), the correction removes part of the real signal.
    % Trade-off deliberately favors conservative claims.
    if flgCorrF0
        for c = {'Cyto', 'Mito'}
            rows  = tbl.compartment == c{1};
            F0ref = median(tbl.F0(rows), 'omitnan');
            scale = F0ref ./ tbl.F0(rows);
            tbl.trace(rows, :) = tbl.trace(rows, :) .* scale;
        end
    end
end

% Clip dF/F at zero. Sub-baseline excursions are noise; keeping them
% biases the flux integral (mito decays integrated below 0) and clutters
% manCur display (mito return-to-baseline is hard to read on a fluctuating
% floor). NaNs are preserved.
tbl.trace = max(tbl.trace, 0);

nKept = nCells;
if flgExclude
    keep = ~tbl.excluded;
    tbl  = tbl(keep, :);
    nKept = sum(keep) / 2;
end

if verbose
    fprintf('  %d rows (%d cells x 2 comp), fs=%.2f Hz, nT=%d', ...
        height(tbl), nKept, fs, nT);
    if flgExclude
        fprintf('  [%d cells dropped]', nCells - nKept);
    end
    fprintf('\n');
end

end     % main


function tf = coerceLogical(v)
%COERCELOGICAL  Map a cell value to logical false unless it parses as true.
if islogical(v)
    tf = v;
elseif isnumeric(v) && isscalar(v) && ~isnan(v)
    tf = v ~= 0;
elseif ischar(v) || isstring(v)
    s = upper(strtrim(char(v)));
    tf = strcmp(s, 'TRUE') || strcmp(s, '1');
else
    tf = false;
end
end
