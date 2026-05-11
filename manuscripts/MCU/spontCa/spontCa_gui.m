function hFig = spontCa_gui(tbl, fs, varargin)
% SPONTCA_GUI Interactive per-cell QC for the SpontCa pipeline.
%
%   hFig = SPONTCA_GUI(TBL, FS, ...) opens a single-cell viewer for the
%   long-format table returned by SPONTCA_EVENTS. Layout (rows top-down):
%
%       row 1: mito trace + mito start/stop markers
%       row 2: cyto trace + cyto start/stop markers (x linked to row 1)
%       row 3: STA aligned to cyto starts (cyto + mito overlay)
%       row 4: cyto dur (hist) | mito dur (hist) | cyto-mito lag (hist)
%
%   Cyto is green, mito is red (cfg.clr.cmp). Event START markers are
%   drawn dashed and STOP markers dotted. Cells are listed in the side
%   dropdown as 'Ctrl_XX' / 'KO_YY' (genotype read from the name).
%
%   INPUTS:
%       tbl     - output of SPONTCA_EVENTS, including the per-row 'map'
%                 column and tbl.Properties.UserData.tWin.
%       fs      - sampling rate (Hz).
%
%   OPTIONAL (Name-Value):
%       'sbjID' - (char) initial cell to display.
%
%   See also: SPONTCA_LOAD, SPONTCA_EVENTS, SPONTCA_DETECT, PLOT_HIST,
%             RIPP_MAPS

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addRequired(p, 'fs',  @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'sbjID', '', @(x) ischar(x) || isstring(x) || isempty(x));
parse(p, tbl, fs, varargin{:});
initSbj = char(p.Results.sbjID);

cfg     = mcu_cfg;
clrCyto = cfg.clr.cmp(1, :);
clrMito = cfg.clr.cmp(2, :);
clrLag  = [0.40, 0.40, 0.40];

nT = size(tbl.trace, 2);
t  = (0:nT-1) / fs;

% STA window time axis (stored by spontCa_events)
if isfield(tbl.Properties.UserData, 'tWin')
    tWin = tbl.Properties.UserData.tWin;
else
    tWin = [];
end


%% ========================================================================
%  FIGURE LAYOUT
%  ========================================================================

hFig = figure('Name', 'SpontCa QC', 'NumberTitle', 'off', 'Color', 'w', ...
    'Units', 'pixels', 'Position', [80, 60, 1400, 900]);

pW = 0.12;
hSide = uipanel('Parent', hFig, 'Units', 'normalized', ...
    'Position', [0, 0, pW, 1], 'BorderType', 'none');
hMain = uipanel('Parent', hFig, 'Units', 'normalized', ...
    'Position', [pW, 0, 1 - pW, 1], 'BorderType', 'none');

tl = tiledlayout(hMain, 4, 3, 'TileSpacing', 'tight', 'Padding', 'compact');
axMito  = nexttile(tl, 1, [1, 3]);
axCyto  = nexttile(tl, 4, [1, 3]);
axSta   = nexttile(tl, 7, [1, 3]);
axDurC  = nexttile(tl, 10);
axDurM  = nexttile(tl, 11);
axLag   = nexttile(tl, 12);


%% ========================================================================
%  SIDE PANEL CONTROL
%  ========================================================================

uicontrol('Parent', hSide, 'Style', 'text', 'String', 'Cell:', ...
    'Units', 'normalized', 'Position', [0.05, 0.94, 0.9, 0.04], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');

cytoMask = tbl.compartment == 'Cyto';
cellList = cellstr(string(tbl.sbjID(cytoMask)));

ddCell = uicontrol('Parent', hSide, 'Style', 'popupmenu', ...
    'String', cellList, ...
    'Units', 'normalized', 'Position', [0.05, 0.90, 0.9, 0.04], ...
    'Callback', @onCellChange);

if ~isempty(initSbj)
    idx = find(strcmp(cellList, initSbj), 1);
    if ~isempty(idx), ddCell.Value = idx; end
end


%% ========================================================================
%  INIT
%  ========================================================================

onCellChange();


%% ========================================================================
%  CALLBACKS
%  ========================================================================

    function onCellChange(~, ~)
        items = get(ddCell, 'String');
        if isempty(items), return; end
        sName = items{get(ddCell, 'Value')};

        iC = find(tbl.sbjID == sName & tbl.compartment == 'Cyto');
        iM = find(tbl.sbjID == sName & tbl.compartment == 'Mito');
        if isempty(iC) || isempty(iM), return; end

        cyTrace = tbl.trace(iC, :);
        miTrace = tbl.trace(iM, :);
        cyStart = tbl.start{iC};
        cyStop  = tbl.stop{iC};
        miStart = tbl.start{iM};
        miStop  = tbl.stop{iM};
        miCoup  = tbl.coupled{iM};
        cyDur   = tbl.dur{iC};
        miDur   = tbl.dur{iM};
        miLag   = tbl.lag{iM};
        cyMap   = tbl.map{iC};
        miMap   = tbl.map{iM};

        % --- MITO TRACE (TOP) ---
        cla(axMito);
        hold(axMito, 'on');
        plot(axMito, t, miTrace, 'Color', clrMito, 'LineWidth', 0.7);
        axis(axMito, 'tight');
        yL = ylim(axMito);
        coupMask = miCoup & miDur > 0;
        plotMarks(axMito, miStart(coupMask), yL, clrMito, '--');
        plotMarks(axMito, miStop(coupMask),  yL, clrMito, ':');
        ylabel(axMito, 'Mito dF/F');
        set(axMito, 'XTickLabel', []);
        title(axMito, sprintf('%s | cyto %d events | mito %d / %d coupled', ...
            sName, length(cyStart), sum(miCoup), length(miStart)), ...
            'Interpreter', 'none');
        hold(axMito, 'off');

        % --- CYTO TRACE (BOTTOM) ---
        cla(axCyto);
        hold(axCyto, 'on');
        plot(axCyto, t, cyTrace, 'Color', clrCyto, 'LineWidth', 0.7);
        axis(axCyto, 'tight');
        yL = ylim(axCyto);
        plotMarks(axCyto, cyStart, yL, clrCyto, '--');
        plotMarks(axCyto, cyStop,  yL, clrCyto, ':');
        ylabel(axCyto, 'Cyto dF/F');
        xlabel(axCyto, 'Time (s)');
        hold(axCyto, 'off');

        linkaxes([axMito, axCyto], 'x');

        % --- STA (overlay cyto + mito averages, cyto-start aligned) ---
        cla(axSta);
        if ~isempty(tWin) && ~isempty(cyMap)
            hold(axSta, 'on');
            plotMeanSEM(axSta, tWin, cyMap, clrCyto, 'Cyto');
            plotMeanSEM(axSta, tWin, miMap, clrMito, 'Mito');
            xline(axSta, 0, 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
            xlabel(axSta, 'Time from cyto start (s)');
            ylabel(axSta, 'dF/F');
            legend(axSta, 'Location', 'northeast', 'Box', 'off');
            title(axSta, 'STA (mean \pm SEM across cyto events)', ...
                'FontWeight', 'normal');
            grid(axSta, 'on');
            hold(axSta, 'off');
        else
            text(axSta, 0.5, 0.5, '(no events for STA)', ...
                'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                'Color', [0.5 0.5 0.5]);
        end

        % --- HISTOGRAMS ---
        plotHistOrBlank(axDurC, cyDur(cyDur > 0), clrCyto, ...
            'Cyto event duration', 'Duration (s)');
        plotHistOrBlank(axDurM, miDur(miDur > 0 & miCoup), clrMito, ...
            'Mito event duration', 'Duration (s)');
        plotHistOrBlank(axLag,  miLag(~isnan(miLag)), clrLag, ...
            'Cyto-Mito lag', 'Lag (s)');
    end


    function plotMarks(ax, x, yL, clr, ls)
        % Batched vertical lines at x positions, drawn in a single PLOT call.
        if isempty(x), return; end
        x  = x(:);
        xx = reshape([x, x, nan(length(x), 1)]', [], 1);
        yy = repmat([yL(1); yL(2); NaN], length(x), 1);
        plot(ax, xx, yy, 'Color', [clr, 0.55], 'LineStyle', ls, ...
            'LineWidth', 0.7, 'HandleVisibility', 'off');
    end


    function plotMeanSEM(ax, tW, mat, clr, lbl)
        nn = sum(~isnan(mat), 1);
        mu = mean(mat, 1, 'omitnan');
        se = std(mat, 0, 1, 'omitnan') ./ sqrt(max(nn, 1));
        xConf = [tW, fliplr(tW)];
        yConf = [mu - se, fliplr(mu + se)];
        vld   = ~isnan(yConf);
        if any(vld)
            fill(ax, xConf(vld), yConf(vld), clr, 'FaceAlpha', 0.25, ...
                'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        plot(ax, tW, mu, 'Color', clr, 'LineWidth', 2, 'DisplayName', lbl);
    end


    function plotHistOrBlank(ax, vals, clr, ttl, xlbl)
        cla(ax);
        if ~isempty(vals) && numel(vals) > 1
            plot_hist([], vals, 'hAx', ax, 'c', clr, ...
                'flgKDE', true, 'flgStat', true);
        else
            text(ax, 0.5, 0.5, '(no events)', 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5]);
        end
        xlabel(ax, xlbl);
        title(ax, ttl, 'FontWeight', 'normal');
    end

end     % EOF
