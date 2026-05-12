%% llmCur_run.m  Workflow for LLM-based curation.
%
% Three-stage pipeline. MATLAB owns rendering and assembly; Claude (as
% the parent agent) dispatches the per-image marking subagents in
% between.
%
%   1. RENDER images (MATLAB):
%        llmCur_render();
%
%      Writes <this_dir>/images/<sbjID>_w<NN>.{png,json}. Default
%      window 120 s, 20% overlap. Adjust via name-value args. PNGs are
%      gitignored.
%
%   2. MARK events (Claude):
%        Open a Claude Code conversation and ask:
%        "Dispatch llmCur marking subagents on llmCur/images/. Save
%         outputs to llmCur/raw/."
%
%      Each subagent reads one PNG + its sidecar JSON and writes
%      <this_dir>/raw/<sbjID>_w<NN>.json with cyto peak times and
%      mito (peak, stop) pairs.
%
%   3. ASSEMBLE per-cell curation files (MATLAB):
%        llmCur_assemble();
%
%      Reads raw JSONs, dedups cross-window duplicates, computes
%      amp/dur/int from the trace, and writes
%      spontCa_curated/<sbjID>_<cmp>.mat in the same format
%      spontCa_manCur produces. Cyto stops are placeholders.
%
% After step 3: run mcu_spontCa.m as usual. spontCa_finalize will
% overlay the curated files automatically. Use spontCa_manCur to
% review/polish per cell.

% Convenience entry: try to run all stages in sequence. Stages 1 and 3
% live in MATLAB and run here; stage 2 must be invoked from Claude
% Code separately - the script just prints a reminder.

llmCur_render();

rawDir = fullfile(fileparts(mfilename('fullpath')), 'raw');
if ~exist(rawDir, 'dir') || isempty(dir(fullfile(rawDir, '*.json')))
    fprintf('\n');
    fprintf('STAGE 2 (Claude): ask it to dispatch llmCur marking\n');
    fprintf('subagents on llmCur/images/, writing JSONs to\n');
    fprintf('llmCur/raw/. Then re-run llmCur_assemble() here.\n');
    return;
end

llmCur_assemble();
