function out = spontCa_tune(varargin)
% SPONTCA_TUNE  Sweep spontCa_detect parameters against curated cells
% in man/ and report per-cell + median F1.
%
% Tunes Cyto and Mito independently. Sweep is a two-pass coordinate
% descent: first sweep (minAmp, kNoise) at default (minIEI, minDur),
% then sweep (minIEI, minDur) at the winning (minAmp, kNoise). Total
% ~70 detect calls per cell per compartment.
%
% USAGE
%   out = spontCa_tune()                       % defaults
%   out = spontCa_tune('cells', {'Ctrl_01','Ctrl_03'})
%   out = spontCa_tune('tolSec', 0.5)          % stricter matching
%
% OPTIONAL (Name-Value):
%   'cells'      - cellstr of sbjIDs to score against. Default: all
%                  files in man/.
%   'tolSec'     - peak-time match tolerance (s). Default 0.7.
%   'minAmpCyto' - sweep values. Default 0.02:0.01:0.10.
%   'minAmpMito' - sweep values. Default 0.01:0.005:0.06.
%   'kNoise'     - sweep values. Default 2.5:0.5:5.0.
%   'minIEI'     - sweep values. Default 0.4:0.2:2.0.
%   'minDur'     - sweep values. Default 0.2:0.2:0.8.
%
% RETURNS
%   out.Cyto, out.Mito - each with:
%     .baseline      - struct {params, perCellF1, medF1, meanF1}
%     .best          - same shape, at the winning params
%     .grid1, .grid2 - tables with all swept combinations + per-cell F1
%
% Prints a tune-report summary. See spontCa_compare for the matching
% logic used here.

%% ARGUMENTS

p = inputParser;
addParameter(p, 'cells', {}, @(x) iscell(x) || ischar(x) || isstring(x));
addParameter(p, 'tolSec', 0.7, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'minAmpCyto', 0.02:0.01:0.10, @isnumeric);
addParameter(p, 'minAmpMito', 0.01:0.005:0.06, @isnumeric);
addParameter(p, 'kNoise', 2.5:0.5:5.0, @isnumeric);
addParameter(p, 'minIEI', 0.4:0.2:2.0, @isnumeric);
addParameter(p, 'minDur', 0.2:0.2:0.8, @isnumeric);
parse(p, varargin{:});
P = p.Results;
if ischar(P.cells) || isstring(P.cells), P.cells = cellstr(P.cells); end

thisDir = fileparts(mfilename('fullpath'));
manDir  = fullfile(thisDir, 'man');
if ~exist(manDir, 'dir')
    error('No man/ folder at %s', manDir);
end

% Discover curated cells.
manFiles = dir(fullfile(manDir, '*.mat'));
manCells = cellfun(@(s) erase(s, '.mat'), {manFiles.name}, ...
    'UniformOutput', false);
if isempty(P.cells)
    cellsToScore = manCells;
else
    cellsToScore = intersect(manCells, P.cells, 'stable');
end
if isempty(cellsToScore)
    error('No matching cells in man/');
end


%% LOAD TRACES + GOLD STARTS

[tblCell, fs] = spontCa_load();
nC = numel(cellsToScore);
traces  = struct('Cyto', cell(nC, 1), 'Mito', cell(nC, 1));
goldStarts = struct('Cyto', cell(nC, 1), 'Mito', cell(nC, 1));

for iC = 1:nC
    sid = cellsToScore{iC};
    iCyto = find(tblCell.sbjID == sid & tblCell.compartment == 'Cyto');
    iMito = find(tblCell.sbjID == sid & tblCell.compartment == 'Mito');
    if isempty(iCyto) || isempty(iMito)
        warning('Skipped %s: missing trace row', sid);
        continue;
    end
    traces(iC).Cyto = tblCell.trace(iCyto, :);
    traces(iC).Mito = tblCell.trace(iMito, :);

    fpath = fullfile(manDir, [sid '.mat']);
    L = load(fpath);
    if isfield(L, 'events')
        events = L.events;
    elseif isfield(L, 'cur') && isfield(L.cur, 'events')
        events = L.cur.events;
    else
        warning('Skipped %s: no events', sid);
        continue;
    end
    goldStarts(iC).Cyto = events.start(events.compartment == 'Cyto');
    goldStarts(iC).Mito = events.start(events.compartment == 'Mito');
end

fprintf('[spontCa_tune] %d curated cells: %s\n', nC, ...
    strjoin(cellsToScore, ', '));


%% BASELINE (current defaults)

defaultCyto = struct('minAmp', 0.05, 'minIEI', 1.0, ...
                     'kNoise', 3.5, 'minDur', 0.4);
defaultMito = struct('minAmp', 0.03, 'minIEI', 1.0, ...
                     'kNoise', 3.5, 'minDur', 0.4);

out.Cyto.baseline = scoreParams(defaultCyto, ...
    {traces.Cyto}, {goldStarts.Cyto}, fs, P.tolSec);
out.Mito.baseline = scoreParams(defaultMito, ...
    {traces.Mito}, {goldStarts.Mito}, fs, P.tolSec);


%% PASS 1: sweep (minAmp, kNoise) at default (minIEI, minDur)

fprintf('[spontCa_tune] Pass 1: minAmp x kNoise ...\n');
out.Cyto.grid1 = sweepGrid('Cyto', P.minAmpCyto, P.kNoise, ...
    defaultCyto.minIEI, defaultCyto.minDur, ...
    {traces.Cyto}, {goldStarts.Cyto}, fs, P.tolSec, cellsToScore);
out.Mito.grid1 = sweepGrid('Mito', P.minAmpMito, P.kNoise, ...
    defaultMito.minIEI, defaultMito.minDur, ...
    {traces.Mito}, {goldStarts.Mito}, fs, P.tolSec, cellsToScore);

bestC1 = pickBest(out.Cyto.grid1);
bestM1 = pickBest(out.Mito.grid1);


%% PASS 2: sweep (minIEI, minDur) at winning (minAmp, kNoise)

fprintf('[spontCa_tune] Pass 2: minIEI x minDur ...\n');
out.Cyto.grid2 = sweepGrid2('Cyto', bestC1.minAmp, bestC1.kNoise, ...
    P.minIEI, P.minDur, ...
    {traces.Cyto}, {goldStarts.Cyto}, fs, P.tolSec, cellsToScore);
out.Mito.grid2 = sweepGrid2('Mito', bestM1.minAmp, bestM1.kNoise, ...
    P.minIEI, P.minDur, ...
    {traces.Mito}, {goldStarts.Mito}, fs, P.tolSec, cellsToScore);

bestC = pickBest(out.Cyto.grid2);
bestM = pickBest(out.Mito.grid2);

out.Cyto.best = struct( ...
    'params', struct('minAmp', bestC.minAmp, 'kNoise', bestC.kNoise, ...
                     'minIEI', bestC.minIEI, 'minDur', bestC.minDur), ...
    'perCellF1', bestC.perCellF1, ...
    'medF1', bestC.medF1, 'meanF1', bestC.meanF1);
out.Mito.best = struct( ...
    'params', struct('minAmp', bestM.minAmp, 'kNoise', bestM.kNoise, ...
                     'minIEI', bestM.minIEI, 'minDur', bestM.minDur), ...
    'perCellF1', bestM.perCellF1, ...
    'medF1', bestM.medF1, 'meanF1', bestM.meanF1);


%% REPORT

fprintf('\n[spontCa_tune] ===== TUNE REPORT (tol=%.2fs) =====\n', P.tolSec);
reportLine('Cyto', out.Cyto, cellsToScore);
reportLine('Mito', out.Mito, cellsToScore);

end


%% =====================================================================

function grid = sweepGrid(cmp, minAmpVec, kNoiseVec, minIEI, minDur, ...
                          tr, gold, fs, tol, sbjIDs)
% Sweep over (minAmp x kNoise) at fixed (minIEI, minDur).
nA = numel(minAmpVec); nK = numel(kNoiseVec);
rows = {};
fprintf('  %s: %d combos x %d cells\n', cmp, nA * nK, numel(tr));
for iA = 1:nA
    for iK = 1:nK
        params = struct('minAmp', minAmpVec(iA), ...
                        'kNoise', kNoiseVec(iK), ...
                        'minIEI', minIEI, 'minDur', minDur);
        s = scoreParams(params, tr, gold, fs, tol);
        rows{end+1} = mergeStructs(params, s); %#ok<AGROW>
    end
end
grid = struct2table(vertcat(rows{:}));
end


function grid = sweepGrid2(cmp, minAmp, kNoise, minIEIvec, minDurVec, ...
                           tr, gold, fs, tol, sbjIDs)
% Sweep (minIEI x minDur) at fixed (minAmp, kNoise).
nI = numel(minIEIvec); nD = numel(minDurVec);
rows = {};
fprintf('  %s: %d combos x %d cells\n', cmp, nI * nD, numel(tr));
for iI = 1:nI
    for iD = 1:nD
        params = struct('minAmp', minAmp, 'kNoise', kNoise, ...
                        'minIEI', minIEIvec(iI), 'minDur', minDurVec(iD));
        s = scoreParams(params, tr, gold, fs, tol);
        rows{end+1} = mergeStructs(params, s); %#ok<AGROW>
    end
end
grid = struct2table(vertcat(rows{:}));
end


function s = scoreParams(params, tr, gold, fs, tol)
% Run spontCa_detect with params on each trace; score F1 against gold.
nC = numel(tr);
f1 = nan(nC, 1);
for iC = 1:nC
    if isempty(tr{iC}) || isempty(gold{iC})
        continue;
    end
    evTbl = spontCa_detect(tr{iC}, fs, ...
        'minAmp', params.minAmp, 'minIEI', params.minIEI, ...
        'kNoise', params.kNoise, 'minDur', params.minDur);
    f1(iC) = scoreF1(gold{iC}, evTbl.start, tol);
end
s = struct( ...
    'perCellF1', f1, ...
    'medF1',  median(f1, 'omitnan'), ...
    'meanF1', mean(f1, 'omitnan'));
end


function F1 = scoreF1(gold, test, tol)
gold = sort(gold(:));
test = sort(test(:));
matched = false(numel(test), 1);
matchedGold = false(numel(gold), 1);
for k = 1:numel(gold)
    d = abs(test - gold(k));
    d(matched) = Inf;
    [dmin, j] = min(d);
    if ~isempty(j) && dmin <= tol
        matched(j) = true;
        matchedGold(k) = true;
    end
end
TP = sum(matched);
FP = numel(test) - TP;
FN = numel(gold) - TP;
P = TP / max(TP + FP, 1);
R = TP / max(TP + FN, 1);
if P + R > 0
    F1 = 2 * P * R / (P + R);
else
    F1 = 0;
end
end


function best = pickBest(grid)
% Pick the row maximizing medF1.
[~, idx] = max(grid.medF1);
best = table2struct(grid(idx, :));
end


function s = mergeStructs(a, b)
s = a;
fb = fieldnames(b);
for k = 1:numel(fb), s.(fb{k}) = b.(fb{k}); end
end


function reportLine(cmp, info, sbjIDs)
b = info.baseline; w = info.best;
fprintf('  %s : baseline med F1 = %.2f | best med F1 = %.2f\n', ...
    cmp, b.medF1, w.medF1);
fprintf('         best params: minAmp=%.3f  kNoise=%.1f  minIEI=%.1f  minDur=%.1f\n', ...
    w.params.minAmp, w.params.kNoise, w.params.minIEI, w.params.minDur);
fprintf('         per-cell baseline -> best:\n');
for k = 1:numel(sbjIDs)
    fprintf('           %s : %.2f -> %.2f\n', ...
        sbjIDs{k}, b.perCellF1(k), w.perCellF1(k));
end
end
