function spontCa_render(varargin)
% SPONTCA_RENDER  Generate per-cell, per-window trace images for
% LLM-based event marking.
%
% Splits each recording into nWindows equal-width overlapping windows
% and renders one PNG per window showing both compartments stacked
% (cyto top, mito bottom). No detector overlay - clean traces only so
% the LLM is not anchored by autodetection output. Each image's title
% encodes the absolute time range so the LLM can map pixel positions
% back to seconds without needing a sidecar metadata file.
%
% USAGE
%   spontCa_render()                       % default params, all cells
%   spontCa_render('cells', {'Ctrl_03'})   % subset
%   spontCa_render('nWindows', 12)         % more windows per cell
%   spontCa_render('overwrite', true)      % re-render existing PNGs
%
% OPTIONAL (Name-Value):
%   'cells'      - cellstr of sbjIDs to render. Default: all in tbl.
%   'nWindows'   - number of windows per cell. Default 10. All windows
%                  have the same width (recDur / (N - (N-1)*overlap));
%                  the last window ends exactly at recDur.
%   'overlap'    - fractional overlap between consecutive windows.
%                  Default 0.2.
%   'imageWH'    - [W H] in pixels for the figure. Default [1600, 900].
%   'outDir'     - directory to write images/<sbjID>_w<NN>.png into.
%                  Default: <spontCa>/llm/.
%   'overwrite'  - re-render images that already exist. Default false.
%
% See also: SPONTCA_LLMCUR, SPONTCA_JSON2MAT, SPONTCA_LOAD, MCU_CFG

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
p.addParameter('cells', {}, @(x) iscell(x) || ischar(x) || isstring(x));
p.addParameter('nWindows', 10, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('overlap', 0.2, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 1);
p.addParameter('imageWH', [1600, 900], ...
    @(x) isnumeric(x) && numel(x) == 2);
p.addParameter('outDir', '', @(x) ischar(x) || isstring(x));
p.addParameter('overwrite', false, @islogical);
parse(p, varargin{:});
P = p.Results;

if ischar(P.cells) || isstring(P.cells)
    P.cells = cellstr(P.cells);
end

% Default output: <spontCa>/llm/ (alongside auto/ and man/).
thisDir = fileparts(mfilename('fullpath'));
if isempty(P.outDir)
    P.outDir = fullfile(thisDir, 'llm');
end
P.outDir = char(P.outDir);
imgDir = fullfile(P.outDir, 'images');
if ~exist(imgDir, 'dir'), mkdir(imgDir); end

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
%  WINDOW GEOMETRY (equal-width, overlapping)
%  ========================================================================
% Solve W*(N - (N-1)*overlap) = recDur for window width W given a fixed
% N. Then step = W*(1-overlap); window k spans [(k-1)*step, (k-1)*step+W].
% Last window ends exactly at recDur.

N    = round(P.nWindows);
W    = recDur / (N - (N - 1) * P.overlap);
step = W * (1 - P.overlap);
fprintf('[spontCa_render] %d windows per cell (W=%.1fs, step=%.1fs, overlap=%.0f%%)\n', ...
    N, W, step, 100 * P.overlap);


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
    for iW = 1:N
        t0 = (iW - 1) * step;
        t1 = t0 + W;
        % Numerical safety: nudge last window to land exactly on recDur.
        if iW == N, t1 = recDur; end

        pngPath = fullfile(imgDir, sprintf('%s_w%02d.png', sName, iW));
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
        nRendered = nRendered + 1;
    end
    fprintf('[spontCa_render] %s : %d windows rendered (of %d)\n', ...
        sName, nRendered, N);
end

end     % EOF
