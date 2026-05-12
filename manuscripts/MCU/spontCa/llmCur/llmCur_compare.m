function out = llmCur_compare(sbjID, varargin)
% LLMCUR_COMPARE  Quantitative match between two curation files for one
% cell. Designed for pilot validation: LLM-curation vs user's manCur
% gold standard.
%
% USAGE
%   out = llmCur_compare('Ctrl_01')                           % defaults
%   out = llmCur_compare('Ctrl_01', 'goldPath', PATH, ...)
%
% By default loads:
%   GOLD: <spontCa>/man/<sbjID>.mat
%   LLM:  <spontCa>/llm/<sbjID>.mat
%
% Reports per compartment: gold count, LLM count, matched, precision,
% recall, F1. Match tolerance default 0.7 s (~2 samples at fs=3).
%
% OPTIONAL (Name-Value):
%   'goldPath' - explicit gold .mat path.
%   'llmPath'  - explicit LLM .mat path.
%   'tolSec'   - peak-time match tolerance in seconds. Default 0.7.
%   'verbose'  - print summary to stdout. Default true.
%
% RETURNS
%   out struct with Cyto, Mito fields; each {nGold, nLLM, TP, FP, FN,
%   precision, recall, F1, unmatchedGold, unmatchedLLM}.

p = inputParser;
addRequired(p, 'sbjID', @(x) ischar(x) || isstring(x));
addParameter(p, 'goldPath', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'llmPath',  '', @(x) ischar(x) || isstring(x));
addParameter(p, 'tolSec', 0.7, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'verbose', true, @islogical);
parse(p, sbjID, varargin{:});
P = p.Results;
sbjID = char(P.sbjID);

thisDir    = fileparts(mfilename('fullpath'));
spontCaDir = fileparts(thisDir);

if isempty(P.goldPath)
    P.goldPath = fullfile(spontCaDir, 'man', [sbjID '.mat']);
end
if isempty(P.llmPath)
    P.llmPath = fullfile(spontCaDir, 'llm', [sbjID '.mat']);
end

[goldCyto, goldMito] = loadStartsByCompartment(P.goldPath);
[llmCyto,  llmMito]  = loadStartsByCompartment(P.llmPath);

out.Cyto = scorePair(goldCyto, llmCyto, P.tolSec);
out.Mito = scorePair(goldMito, llmMito, P.tolSec);

if P.verbose
    reportLine(sbjID, 'Cyto', out.Cyto, P.goldPath, P.llmPath);
    reportLine(sbjID, 'Mito', out.Mito, P.goldPath, P.llmPath);
end

end


function [sCyto, sMito] = loadStartsByCompartment(fpath)
S = load(fpath, 'cur');
events = S.cur.events;
sCyto = events.start(events.compartment == 'Cyto');
sMito = events.start(events.compartment == 'Mito');
end


function r = scorePair(gold, llm, tol)
gold = sort(gold);
llm  = sort(llm);
nG = numel(gold); nL = numel(llm);
matched = false(nL, 1);
matchedGold = false(nG, 1);
for k = 1:nG
    d = abs(llm - gold(k));
    d(matched) = Inf;
    [dmin, j] = min(d);
    if ~isempty(j) && dmin <= tol
        matched(j) = true;
        matchedGold(k) = true;
    end
end
TP = sum(matched);
FP = nL - TP;
FN = nG - TP;
precision = TP / max(TP + FP, 1);
recall    = TP / max(TP + FN, 1);
if precision + recall > 0
    F1 = 2 * precision * recall / (precision + recall);
else
    F1 = 0;
end
r = struct('nGold', nG, 'nLLM', nL, 'TP', TP, 'FP', FP, 'FN', FN, ...
    'precision', precision, 'recall', recall, 'F1', F1, ...
    'unmatchedGold', gold(~matchedGold), ...
    'unmatchedLLM',  llm(~matched));
end


function reportLine(sbjID, cmp, r, gPath, lPath)
[~, gName] = fileparts(gPath);
[~, lName] = fileparts(lPath);
fprintf('[%s %s] gold=%d  llm=%d  TP=%d  FP=%d  FN=%d  P=%.2f  R=%.2f  F1=%.2f\n', ...
    sbjID, cmp, r.nGold, r.nLLM, r.TP, r.FP, r.FN, ...
    r.precision, r.recall, r.F1);
fprintf('   gold: %s\n   llm:  %s\n', gName, lName);
if ~isempty(r.unmatchedGold)
    fprintf('   missed by LLM (t,s): %s\n', ...
        strjoin(arrayfun(@(x) sprintf('%.1f', x), ...
        r.unmatchedGold, 'UniformOutput', false), ', '));
end
if ~isempty(r.unmatchedLLM)
    fprintf('   extra in LLM (t,s):  %s\n', ...
        strjoin(arrayfun(@(x) sprintf('%.1f', x), ...
        r.unmatchedLLM, 'UniformOutput', false), ', '));
end
end
