function out = spontCa_compare(sbjID, varargin)
% SPONTCA_COMPARE  Quantitative match between two curation files for
% one cell. Generic over the source being tested: pass any .mat in the
% events-table format as 'testPath' (LLM output, autodetector output,
% an older curation, etc.) and compare against the gold standard.
%
% USAGE
%   out = spontCa_compare('Ctrl_01')                          % defaults
%   out = spontCa_compare('Ctrl_01', 'testPath', PATH, ...)
%   out = spontCa_compare('Ctrl_01', 'goldPath', G, 'testPath', T)
%
% By default loads:
%   GOLD: <spontCa>/man/<sbjID>.mat
%   TEST: <spontCa>/llm/<sbjID>.mat
%
% Reports per compartment: gold count, test count, matched, precision,
% recall, F1. Match tolerance default 0.7 s (~2 samples at fs=3).
%
% OPTIONAL (Name-Value):
%   'goldPath' - explicit gold .mat path.
%   'testPath' - explicit test .mat path.
%   'tolSec'   - peak-time match tolerance in seconds. Default 0.7.
%   'verbose'  - print summary to stdout. Default true.
%
% RETURNS
%   out struct with Cyto, Mito fields; each {nGold, nTest, TP, FP, FN,
%   precision, recall, F1, unmatchedGold, unmatchedTest}.

p = inputParser;
addRequired(p, 'sbjID', @(x) ischar(x) || isstring(x));
addParameter(p, 'goldPath', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'testPath', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'tolSec', 0.7, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'verbose', true, @islogical);
parse(p, sbjID, varargin{:});
P = p.Results;
sbjID = char(P.sbjID);

thisDir = fileparts(mfilename('fullpath'));

if isempty(P.goldPath)
    P.goldPath = fullfile(thisDir, 'man', [sbjID '.mat']);
end
if isempty(P.testPath)
    P.testPath = fullfile(thisDir, 'llm', [sbjID '.mat']);
end

[goldCyto, goldMito] = loadStartsByCompartment(P.goldPath);
[testCyto, testMito] = loadStartsByCompartment(P.testPath);

out.Cyto = scorePair(goldCyto, testCyto, P.tolSec);
out.Mito = scorePair(goldMito, testMito, P.tolSec);

if P.verbose
    reportLine(sbjID, 'Cyto', out.Cyto, P.goldPath, P.testPath);
    reportLine(sbjID, 'Mito', out.Mito, P.goldPath, P.testPath);
end

end


function [sCyto, sMito] = loadStartsByCompartment(fpath)
S = load(fpath);
if isfield(S, 'events')
    events = S.events;
elseif isfield(S, 'cur') && isfield(S.cur, 'events')
    events = S.cur.events;  % legacy v4 wrapping
else
    error('spontCa_compare:badFile', ...
        'File missing events variable: %s', fpath);
end
sCyto = events.start(events.compartment == 'Cyto');
sMito = events.start(events.compartment == 'Mito');
end


function r = scorePair(gold, test, tol)
gold = sort(gold);
test = sort(test);
nG = numel(gold); nT = numel(test);
matched = false(nT, 1);
matchedGold = false(nG, 1);
for k = 1:nG
    d = abs(test - gold(k));
    d(matched) = Inf;
    [dmin, j] = min(d);
    if ~isempty(j) && dmin <= tol
        matched(j) = true;
        matchedGold(k) = true;
    end
end
TP = sum(matched);
FP = nT - TP;
FN = nG - TP;
precision = TP / max(TP + FP, 1);
recall    = TP / max(TP + FN, 1);
if precision + recall > 0
    F1 = 2 * precision * recall / (precision + recall);
else
    F1 = 0;
end
r = struct('nGold', nG, 'nTest', nT, 'TP', TP, 'FP', FP, 'FN', FN, ...
    'precision', precision, 'recall', recall, 'F1', F1, ...
    'unmatchedGold', gold(~matchedGold), ...
    'unmatchedTest', test(~matched));
end


function reportLine(sbjID, cmp, r, gPath, tPath)
[~, gName] = fileparts(gPath);
[~, tName] = fileparts(tPath);
fprintf('[%s %s] gold=%d  test=%d  TP=%d  FP=%d  FN=%d  P=%.2f  R=%.2f  F1=%.2f\n', ...
    sbjID, cmp, r.nGold, r.nTest, r.TP, r.FP, r.FN, ...
    r.precision, r.recall, r.F1);
fprintf('   gold: %s\n   test: %s\n', gName, tName);
if ~isempty(r.unmatchedGold)
    fprintf('   missed by test (t,s): %s\n', ...
        strjoin(arrayfun(@(x) sprintf('%.1f', x), ...
        r.unmatchedGold, 'UniformOutput', false), ', '));
end
if ~isempty(r.unmatchedTest)
    fprintf('   extra in test (t,s):  %s\n', ...
        strjoin(arrayfun(@(x) sprintf('%.1f', x), ...
        r.unmatchedTest, 'UniformOutput', false), ', '));
end
end
