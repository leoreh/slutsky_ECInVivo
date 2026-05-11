function hFig = spontCa_gui(tbl, fs, varargin)
% SPONTCA_GUI Interactive per-cell QC for the SpontCa pipeline.
%
%   hFig = SPONTCA_GUI(TBL, FS, ...) opens a single-cell viewer for the
%   long-format table returned by SPONTCA_EVENTS. Layout (top to bottom):
%
%       cyto trace + cyto start/stop markers
%       mito trace + mito start/stop markers  (x linked to cyto, near-zero gap)
%       ETA  |  cyto dur  |  mito dur  |  cyto-mito lag
%
%   Cyto and mito colors come from cfg.clr.cmp. Event START markers are
%   drawn dashed and STOP markers dotted. The dropdown lists every cell
%   as 'Ctrl_XX' / 'KO_YY' (genotype read from the name). Axes use manual
%   positioning so the cyto/mito gap is ~1 px while the gap above the
%   histogram row stays comfortable.
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

if isfield(tbl.Properties.UserData, 'tWin')
    tWin = tbl.Properties.UserData.tWin;
else
    tWin = [];
end


%% ========================================================================
%  FIGURE LAYOUT (manual axes positions)
%  ========================================================================

hFig = figure('Name', 'SpontCa QC', 'NumberTitle', 'off', 'Color', 'w', ...
    'Units', 'pixels', 'Position', [80, 60, 1400, 900]);

pW = 0.12;
hSide = uipanel('Parent', hFig, 'Units', 'normalized', ...
    'Position', [0, 0, pW, 1], 'BorderType', 'none');
hMain = uipanel('Parent', hFig, 'Units', 'normalized', ...
    'Position', [pW, 0, 1 - pW, 1], 'BorderType', 'none');

% Geometry inside hMain
margL = 0.07;  margR = 0.02;
margT = 0.06;  margB = 0.09;
gapBig = 0.07;     % between mito trace and bottom row
gapTiny = 0.005;   % between cyto and mito (near-zero)
gapH = 0.045;      % horizontal between bottom-row panels

availW = 1 - margL - margR;
hTrace = 0.20;     % each trace tile
hRow   = 0.30;     % bottom row (ETA + 3 hists)
nCols  = 4;
panelW = (availW - (nCols - 1) * gapH) / nCols;

yRow  = margB;
yMito = yRow  + hRow   + gapBig;
yCyto = yMito + hTrace + gapTiny;

axCyto = axes('Parent', hMain, 'Position', [margL, yCyto, availW, hTrace]);
axMito = axes('Parent', hMain, 'Position', [margL, yMito, availW, hTrace]);
axEta  = axes('Parent', hMain, 'Position', ...
    [margL + 0 * (panelW + gapH), yRow, panelW, hRow]);
axDurC = axes('Parent', hMain, 'Position', ...
    [margL + 1 * (panelW + gapH), yRow, panelW, hRow]);
axDurM = axes('Parent', hMain, 'Position', ...
    [margL + 2 * (panelW + gapH), yRow, panelW, hRow]);
axLag  = axes('Parent', hMain, 'Position', ...
    [margL + 3 * (panelW + gapH), yRow, panelW, hRow]);


%% ========================================================================
%  SIDE PANEL CONTROL
%  ========================================================================

uicontrol('Parent', hSide, 'Style', 'text', 'String', 'Cell:', ...
    'Units', 'normalized', 'Position', [0.05, 0.94, 0.9, 0.04], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');

cellList = cellstr(string(tbl.sbjID(tbl.compartment == 'Cyto')));

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

        % --- CYTO TRACE (TOP) ---
        cla(axCyto, 'reset');
        hold(axCyto, 'on');
        plot(axCyto, t, cyTrace, 'Color', clrCyto, 'LineWidth', 0.7);
        axis(axCyto, 'tight');
        yL = ylim(axCyto);
        plotMarks(axCyto, cyStart, yL, clrCyto, '--');
        plotMarks(axCyto, cyStop,  yL, clrCyto, ':');
        ylabel(axCyto, 'Cyto dF/F');
        set(axCyto, 'XTickLabel', []);
        title(axCyto, sprintf('%s | cyto %d events | mito %d / %d coupled', ...
            sName, length(cyStart), sum(miCoup), length(miStart)), ...
            'Interpreter', 'none');
        hold(axCyto, 'off');

        % --- MITO TRACE (BOTTOM, x-linked to cyto) ---
        cla(axMito, 'reset');
        hold(axMito, 'on');
        plot(axMito, t, miTrace, 'Color', clrMito, 'LineWidth', 0.7);
        axis(axMito, 'tight');
        yL = ylim(axMito);
        coupMask = miCoup & miDur > 0;
        plotMarks(axMito, miStart(coupMask), yL, clrMito, '--');
        plotMarks(axMito, miStop(coupMask),  yL, clrMito, ':');
        ylabel(axMito, 'Mito dF/F');
        xlabel(axMito, 'Time (s)');
        hold(axMito, 'off');

        linkaxes([axCyto, axMito], 'x');

        % --- ETA (overlay cyto + mito averages, cyto-start aligned) ---
        cla(axEta, 'reset');
        if ~isempty(tWin) && ~isempty(cyMap)
            hold(axEta, 'on');
            plotMeanSEM(axEta, tWin, cyMap, clrCyto, 'Cyto');
            plotMeanSEM(axEta, tWin, miMap, clrMito, 'Mito');
            yL = ylim(axEta);
            plot(axEta, [0 0], yL, 'k:', 'LineWidth', 0.8, ...
                'HandleVisibility', 'off');
            xlabel(axEta, 'Time from cyto start (s)');
            ylabel(axEta, 'dF/F');
            legend(axEta, 'Location', 'northeast', 'Box', 'off');
            title(axEta, 'ETA', 'FontWeight', 'normal');
            grid(axEta, 'on');
            hold(axEta, 'off');
        else
            text(axEta, 0.5, 0.5, '(no events for ETA)', ...
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
        % Batched vertical lines drawn in a single PLOT call.
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
        cla(ax, 'reset');
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
