%% llmCur_run.m  Workflow for LLM-based curation.
%
% Three-stage pipeline. MATLAB owns rendering and assembly; Claude (the
% parent agent) dispatches per-image marking subagents in between.
%
%   1. RENDER images (MATLAB):
%        llmCur_render();
%
%      Writes <spontCa>/llm/images/<sbjID>_w<NN>.{png,json}. Default
%      10 windows per cell, 20% overlap, equal width. Adjust via
%      name-value args.
%
%   2. MARK events (Claude):
%        Open a Claude Code conversation and ask:
%        "Dispatch llmCur marking subagents on llm/images/. Save
%         outputs to llm/raw/."
%
%      Each subagent reads one PNG + its sidecar JSON and writes
%      <spontCa>/llm/raw/<sbjID>_w<NN>.json with cyto peak times and
%      mito (peak, stop) pairs.
%
%   3. ASSEMBLE per-cell curation files (MATLAB):
%        llmCur_assemble();
%
%      Reads raw JSONs, dedups cross-window duplicates, computes
%      amp/dur/int from the trace, and writes
%      <spontCa>/llm/<sbjID>.mat (single file per cell, events table
%      with both compartments). Cyto stops are placeholders.
%
% After step 3: run mcu_spontCa.m as usual; spontCa_finalize will
% overlay <spontCa>/man/<sbjID>.mat if present. To inspect an LLM
% curation in manCur, hit Load and navigate to the llm/ folder. To
% accept LLM as gold, copy llm/<sbjID>.mat -> man/<sbjID>.mat (or open
% it in manCur and Save, which writes to man/).

llmCur_render();

llmDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'llm');
rawDir = fullfile(llmDir, 'raw');
if ~exist(rawDir, 'dir') || isempty(dir(fullfile(rawDir, '*.json')))
    fprintf('\n');
    fprintf('STAGE 2 (Claude): ask it to dispatch llmCur marking\n');
    fprintf('subagents on %s, writing JSONs to %s.\n', ...
        fullfile(llmDir, 'images'), rawDir);
    fprintf('Then re-run llmCur_assemble() here.\n');
    return;
end

llmCur_assemble();
