function [tbl, fs] = spontCa_load(varargin)
% SPONTCA_LOAD Reads SpontCa.xlsx into a long-format unit-level table.
%
%   [tbl, fs] = SPONTCA_LOAD(...) returns one row per (cell x compartment).
%   The table is the canonical input to SPONTCA_EVENTS and SPONTCA_GUI.
%   Sampling rate is returned alongside; the time vector is derivable as
%   (0:nT-1)/fs and is intentionally not stored in the table.
%
%   OUTPUT SCHEMA:
%       genotype     (categorical) 'Control' / 'MCU-KO'
%       sbjID        (categorical) per-cell identifier ('Ctrl_01', 'KO_06')
%       unitID       (categorical) raw Excel cell number
%       compartment  (categorical) 'Cyto' / 'Mito'
%       trace        (1 x nT)      dF/F signal
%
%   OPTIONAL (Name-Value):
%       'xlsPath' - (char) Path to Excel file. Defaults to canonical
%                   location under D:\OneDrive ...\Data_NetaF\SpontCa.xlsx.
%       'verbose' - (log)  Print progress {true}.
%
%   FILE FORMAT (rows 1-indexed):
%       row 1   : genotype label, repeated over each Cyto/Mito pair
%       row 2   : cell ID number above each Cyto col; missing over Mito
%       row 3   : signal type ('Cyto' / 'Mito')
%       row 4+  : data. col 1 is time (s); col 2+ are signal columns
%
%   See also: SPONTCA_DETECT, SPONTCA_EVENTS, SPONTCA_GUI

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addParameter(p, 'xlsPath', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'verbose', true, @islogical);
parse(p, varargin{:});

xlsPath = char(p.Results.xlsPath);
verbose = p.Results.verbose;

if isempty(xlsPath)
    xlsPath = fullfile('D:\OneDrive - Tel-Aviv University', ...
        'PhD\Slutsky\Manuscripts\MCU\Results\Data_NetaF', 'SpontCa.xlsx');
end

assert(isfile(xlsPath), 'spontCa_load:fileNotFound', ...
    'Excel file not found: %s', xlsPath);

if verbose
    fprintf('[spontCa_load] Reading %s\n', xlsPath);
end


%% ========================================================================
%  READ EXCEL
%  ========================================================================

raw = readcell(xlsPath);
[nRowsRaw, nColsRaw] = size(raw);

% Headers (col 1 is time; signal cols start at 2)
nSig = nColsRaw - 1;
hdrGen = strings(1, nSig);
hdrCid = nan(1, nSig);
hdrSig = strings(1, nSig);
for iC = 1:nSig
    val = raw{1, iC + 1};
    if ~ismissing(val), hdrGen(iC) = string(val); end
    val = raw{2, iC + 1};
    if isnumeric(val) && ~isnan(val), hdrCid(iC) = val; end
    val = raw{3, iC + 1};
    if ~ismissing(val), hdrSig(iC) = string(val); end
end

% Numeric data block, padding non-numeric entries with NaN.
dataMat = nan(nRowsRaw - 3, nColsRaw);
for iR = 4:nRowsRaw
    for iC = 1:nColsRaw
        val = raw{iR, iC};
        if isnumeric(val) && ~isempty(val) && ~ismissing(val)
            dataMat(iR - 3, iC) = val;
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
traces      = nan(n, nT);

iRow = 0;
for iC = 1:nCells
    cC   = cytoIdx(iC);
    cM   = cC + 1;
    gen  = char(hdrGen(cC));
    cnum = hdrCid(cC);

    % Short sbjID: 'Ctrl_XX' for Control, 'KO_YY' for MCU-KO
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
    traces(iRow, :)   = sig(:, cC)';

    iRow = iRow + 1;
    genotype(iRow)    = gen;
    sbjID(iRow)       = sName;
    compartment(iRow) = 'Mito';
    unitID(iRow)      = cnum;
    if cM <= nSig && strcmpi(hdrSig(cM), 'Mito') && hdrGen(cM) == hdrGen(cC)
        traces(iRow, :) = sig(:, cM)';
    else
        warning('spontCa_load:noMito', ...
            'No paired Mito column for cell %d (%s).', cnum, gen);
    end
end

cfg = mcu_cfg;
tbl = table(...
    categorical(genotype, cfg.lbl.grp), ...
    categorical(sbjID), ...
    categorical(unitID), ...
    categorical(compartment, {'Cyto', 'Mito'}), ...
    traces, ...
    'VariableNames', ...
    {'genotype', 'sbjID', 'unitID', 'compartment', 'trace'});

if verbose
    fprintf('  %d rows (%d cells x 2 compartments), fs=%.2f Hz, nT=%d\n', ...
        height(tbl), nCells, fs, nT);
end

end     % EOF
