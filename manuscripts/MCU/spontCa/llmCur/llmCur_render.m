function llmCur_render(varargin)
% LLMCUR_RENDER  Generate per-cell, per-window trace images for LLM-based
% event marking.
%
% For each cell, slides a fixed window (default 120 s with 20% overlap)
% across the recording and renders one PNG per window showing both
% compartments stacked (cyto top, mito bottom). No detector overlay -
% clean traces only so the LLM is not anchored by autodetection output.
%
% USAGE
%   llmCur_render()                       % default params, all cells
%   llmCur_render('cells', {'Ctrl_03'})   % subset
%   llmCur_render('windowSec', 180)       % wider windows
%   llmCur_render('overwrite', true)      % re-render existing PNGs
%
% OPTIONAL (Name-Value):
%   'cells'      - cellstr of sbjIDs to render. Default: all in tbl.
%   'windowSec'  - window width in seconds. Default 120.
%   'overlap'    - fractional overlap with previous window. Default 0.2.
%   'imageWH'    - [W H] in pixels for the figure. Default [1600, 900].
%   'outDir'     - root output dir; PNGs go in <outDir>/images/.
%                  Default: alongside this file.
%   'overwrite'  - re-render images that already exist. Default false.
%
% OUTPUTS (under outDir/images/)
%   <sbjID>_w<NN>.png   - rendered trace pair
%   <sbjID>_w<NN>.json  - sidecar metadata (sbjID, fs, t_start, t_end)
%
% See also: LLMCUR_RUN, LLMCUR_ASSEMBLE, SPONTCA_LOAD, MCU_CFG

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
p.addParameter('cells', {}, @(x) iscell(x) || ischar(x) || isstring(x));
p.addParameter('windowSec', 120, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('overlap', 0.2, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 1);
p.addParameter('imageWH', [1600, 900], ...
    @(x) isnumeric(x) && numel(x) == 2);
p.addParameter('outDir', '', @(x) ischar(x) || isstring(x));
p.addParameter('overwrite', false, @islogical);
parse(p, varargin{:});
P = p.Results;

% Normalize cells argument
if ischar(P.cells) || isstring(P.cells)
    P.cells = cellstr(P.cells);
end

if isempty(P.outDir)
    P.outDir = fileparts(mfilename('fullpath'));
end
P.outDir = char(P.outDir);
imgDir = fullfile(P.outDir, 'images');
if ~exist(imgDir, 'dir'), mkdir(imgDir); end


%% ========================================================================
%  LOAD
%  ========================================================================
% Add the parent spontCa directory to path so spontCa_load / mcu_cfg
% resolve regardless of caller cwd.

spontCaDir = fileparts(P.outDir);
if exist(spontCaDir, 'dir') && ~contains(lower(path), lower(spontCaDir))
    addpath(spontCaDir);
end

[tbl, fs] = spontCa_load();
nT     = size(tbl.trace, 2);
recDur = (nT - 1) / fs;
tVec   = (0:nT-1) / fs;


%% ========================================================================
%  CELL LIST
%  ========================================================================

allCells = unique(cellstr(string(tbl.sbjID)), 'stable');
if isempty(P.cells)
    cellsToRender = allCells;
else
    cellsToRender = intersect(allCells, P.cells, 'stable');
    missing = setdiff(P.cells, allCells);
    if ~isempty(missing)
        warning('Cells not in table: %s', strjoin(missing, ', '));
    end
end


%% ========================================================================
%  WINDOW GEOMETRY
%  ========================================================================

W    = P.windowSec;
step = W * (1 - P.overlap);
nW   = max(1, ceil((recDur - W) / step) + 1);
fprintf('[llmCur_render] %d windows per cell (W=%.0fs, step=%.0fs)\n', ...
    nW, W, step);


%% ========================================================================
%  COLORS
%  ========================================================================

cfg     = mcu_cfg;
clrCyto = cfg.clr.cmp(1, :);
clrMito = cfg.clr.cmp(2, :);


%% ========================================================================
%  RENDER
%  ========================================================================

for iCell = 1:numel(cellsToRender)
    sName = cellsToRender{iCell};
    iC = find(tbl.sbjID == sName & tbl.compartment == 'Cyto');
    iM = find(tbl.sbjID == sName & tbl.compartment == 'Mito');
    if isempty(iC) || isempty(iM)
        warning('Cell %s missing cyto or mito row; skipped', sName);
        continue;
    end
    cyTrace = tbl.trace(iC, :);
    miTrace = tbl.trace(iM, :);

    nRendered = 0;
    for iW = 1:nW
        t0 = (iW - 1) * step;
        t1 = min(recDur, t0 + W);
        if t0 >= recDur, break; end

        pngPath  = fullfile(imgDir, sprintf('%s_w%02d.png',  sName, iW));
        jsonPath = fullfile(imgDir, sprintf('%s_w%02d.json', sName, iW));
        if exist(pngPath, 'file') && ~P.overwrite
            continue;
        end

        mask  = tVec >= t0 & tVec <= t1;
        tWin  = tVec(mask);
        cyWin = cyTrace(mask);
        miWin = miTrace(mask);

        f = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
            'Position', [50, 50, P.imageWH(1), P.imageWH(2)], ...
            'InvertHardcopy', 'off');

        axCyto = axes('Parent', f, 'Position', [0.06, 0.55, 0.91, 0.40]);
        axMito = axes('Parent', f, 'Position', [0.06, 0.10, 0.91, 0.40]);

        plot(axCyto, tWin, cyWin, 'Color', clrCyto, 'LineWidth', 1.0);
        xlim(axCyto, [t0, t1]);
        ylabel(axCyto, 'cyto dF/F');
        set(axCyto, 'XTickLabel', []);
        title(axCyto, sprintf('%s | t = %.1f - %.1f s', ...
            sName, t0, t1), 'Interpreter', 'none', 'FontWeight', 'normal');

        plot(axMito, tWin, miWin, 'Color', clrMito, 'LineWidth', 1.0);
        xlim(axMito, [t0, t1]);
        ylabel(axMito, 'mito dF/F');
        xlabel(axMito, 'Time (s)');

        exportgraphics(f, pngPath, 'Resolution', 100, ...
            'BackgroundColor', 'white');
        close(f);

        meta = struct( ...
            'sbjID',      sName, ...
            'window_idx', iW, ...
            'fs',         fs, ...
            't_start',    t0, ...
            't_end',      t1, ...
            'png',        pngPath);
        fjson = fopen(jsonPath, 'w');
        fprintf(fjson, '%s\n', jsonencode(meta));
        fclose(fjson);

        nRendered = nRendered + 1;
    end
    fprintf('[llmCur_render] %s : %d windows rendered (of %d)\n', ...
        sName, nRendered, nW);
end

end     % EOF
