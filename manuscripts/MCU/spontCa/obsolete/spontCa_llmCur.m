%% spontCa_llmCur.m  Workflow for LLM-based curation.
%
% Three-stage pipeline. MATLAB owns rendering and assembly; Claude
% (the parent agent) dispatches per-image marking subagents in between.
%
%   1. RENDER images (MATLAB):
%        spontCa_render();
%
%      Writes <spontCa>/llm/images/<sbjID>_w<NN>.png. Default 10
%      equal-width windows per cell with 20% overlap. Each image
%      encodes its absolute time range in the title - no sidecar
%      metadata file.
%
%   2. MARK events (Claude):
%        Open a Claude Code conversation and ask:
%        "Dispatch marking subagents on llm/images/. Save outputs
%         to llm/raw/."
%
%      Each subagent reads one PNG (reads the time range from the
%      title) and writes llm/raw/<sbjID>_w<NN>.json with cyto peak
%      times and mito (peak, stop) pairs in absolute seconds.
%
%   3. ASSEMBLE per-cell curation files (MATLAB):
%        spontCa_json2mat();
%
%      Reads raw JSONs, dedups cross-window duplicates, computes
%      amp/dur/int from the trace, and writes llm/<sbjID>.mat (one
%      file per cell, events table with both compartments).
%
% After step 3: run mcu_spontCa.m as usual; spontCa_finalize will
% overlay <spontCa>/man/<sbjID>.mat if present. To inspect an LLM
% curation in manCur, click Load and navigate to llm/. To promote
% LLM to gold for a cell, open it in manCur and click Save (writes
% to man/<sbjID>.mat).

spontCa_render();

llmDir = fullfile(fileparts(mfilename('fullpath')), 'llm');
rawDir = fullfile(llmDir, 'raw');
if ~exist(rawDir, 'dir') || isempty(dir(fullfile(rawDir, '*.json')))
    fprintf('\n');
    fprintf('STAGE 2 (Claude): ask it to dispatch marking subagents on\n');
    fprintf('%s, writing JSONs to %s.\n', ...
        fullfile(llmDir, 'images'), rawDir);
    fprintf('Then re-run spontCa_json2mat() here.\n');
    return;
end

spontCa_json2mat();
