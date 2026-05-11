function hFig = spontCa_gui(tbl, fs, varargin)
% SPONTCA_GUI Interactive per-cell QC for the SpontCa pipeline.
%
%   hFig = SPONTCA_GUI(TBL, FS, ...) opens a single-cell viewer for the
%   long-format table returned by SPONTCA_EVENTS + SPONTCA_COUPLE.
%
%   LAYOUT (top to bottom):
%       cyto trace + cyto start/stop markers
%       mito trace + mito start/stop markers  (x linked to cyto)
%       etaCyto | etaMito | durCyto | durMito | lag
%
%   ETA panels overlay cyto (red) and mito (green) ETAs, the first aligned
%   on cyto starts (mapCyto), the second aligned on mito starts (mapMito).
%   Cyto-independent mito events are drawn with reduced alpha so they stay
%   visually distinct without changing color.
%
%   INPUTS:
%       tbl - output of SPONTCA_EVENTS + SPONTCA_COUPLE, including the
%             per-row 'mapCyto' / 'mapMito' columns and the per-row
%             'cytoIndependent' / 'lag' columns on mito rows.
%       fs  - sampling rate (Hz).
%
%   OPTIONAL (Name-Value):
%       'sbjID' - (char) initial cell to display.
%
%   See also: SPONTCA_LOAD, SPONTCA_EVENTS, SPONTCA_COUPLE, SPONTCA_DETECT,
%             PLOT_HIST

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

alphaCoup  = 0.55;  % cyto-coupled mito events
alphaIndep = 0.25;  % cyto-independent mito events

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
    'Units', 'pixels', 'Position', [80, 60, 1500, 900]);

pW = 0.11;
hSide = uipanel('Parent', hFig, 'Units', 'normalized', ...
    'Position', [0, 0, pW, 1], 'BorderType', 'none');
hMain = uipanel('Parent', hFig, 'Units', 'normalized', ...
    'Position', [pW, 0, 1 - pW, 1], 'BorderType', 'none');

% Geometry inside hMain
margL = 0.07;  margR = 0.02;
margT = 0.06;  margB = 0.09;
gapBig = 0.07;     % between mito trace and bottom row
gapTiny = 0.005;   % between cyto and mito (near-zero)
gapH = 0.035;      % horizontal between bottom-row panels

availW = 1 - margL - margR;
hTrace = 0.20;     % each trace tile
hRow   = 0.30;     % bottom row (2 ETAs + 3 hists)
nCols  = 5;
panelW = (availW - (nCols - 1) * gapH) / nCols;

yRow  = margB;
yMito = yRow  + hRow   + gapBig;
yCyto = yMito + hTrace + gapTiny;

axCyto    = axes('Parent', hMain, 'Position', [margL, yCyto, availW, hTrace]);
axMito    = axes('Parent', hMain, 'Position', [margL, yMito, availW, hTrace]);
axEtaCyto = axes('Parent', hMain, 'Position', ...
    [margL + 0 * (panelW + gapH), yRow, panelW, hRow]);
axEtaMito = axes('Parent', hMain, 'Position', ...
    [margL + 1 * (panelW + gapH), yRow, panelW, hRow]);
axDurC    = axes('Parent', hMain, 'Position', ...
    [margL + 2 * (panelW + gapH), yRow, panelW, hRow]);
axDurM    = axes('Parent', hMain, 'Position', ...
    [margL + 3 * (panelW + gapH), yRow, panelW, hRow]);
axLag     = axes('Parent', hMain, 'Position', ...
    [margL + 4 * (panelW + gapH), yRow, panelW, hRow]);


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

        cyTrace   = tbl.trace(iC, :);
        miTrace   = tbl.trace(iM, :);
        cyStart   = tbl.start{iC};
        cyStop    = tbl.stop{iC};
        miStart   = tbl.start{iM};
        miStop    = tbl.stop{iM};
        miIndep   = tbl.cytoIndependent{iM};
        cyDur     = tbl.dur{iC};
        miDur     = tbl.dur{iM};
        miLag     = tbl.lag{iM};
        cyMapC    = tbl.mapCyto{iC};
        miMapC    = tbl.mapCyto{iM};
        cyMapM    = tbl.mapMito{iC};
        miMapM    = tbl.mapMito{iM};

        nIndep = sum(miIndep);

        % --- CYTO TRACE (TOP) ---
        cla(axCyto, 'reset');
        hold(axCyto, 'on');
        plot(axCyto, t, cyTrace, 'Color', clrCyto, 'LineWidth', 0.7);
        axis(axCyto, 'tight');
        yL = ylim(axCyto);
        plotMarks(axCyto, cyStart, yL, clrCyto, '--', alphaCoup);
        plotMarks(axCyto, cyStop,  yL, clrCyto, ':',  alphaCoup);
        ylabel(axCyto, 'Cyto dF/F');
        set(axCyto, 'XTickLabel', []);
        title(axCyto, sprintf('%s | cyto %d events | mito %d events (%d independent)', ...
            sName, length(cyStart), length(miStart), nIndep), ...
            'Interpreter', 'none');
        hold(axCyto, 'off');

        % --- MITO TRACE (BOTTOM, x-linked to cyto) ---
        cla(axMito, 'reset');
        hold(axMito, 'on');
        plot(axMito, t, miTrace, 'Color', clrMito, 'LineWidth', 0.7);
        axis(axMito, 'tight');
        yL = ylim(axMito);
        mCoup = ~miIndep;
        plotMarks(axMito, miStart(mCoup),  yL, clrMito, '--', alphaCoup);
        plotMarks(axMito, miStop(mCoup),   yL, clrMito, ':',  alphaCoup);
        plotMarks(axMito, miStart(miIndep), yL, clrMito, '--', alphaIndep);
        plotMarks(axMito, miStop(miIndep),  yL, clrMito, ':',  alphaIndep);
        ylabel(axMito, 'Mito dF/F');
        xlabel(axMito, 'Time (s)');
        hold(axMito, 'off');

        linkaxes([axCyto, axMito], 'x');

        % --- ETA cyto-aligned ---
        plotEta(axEtaCyto, tWin, cyMapC, miMapC, clrCyto, clrMito, ...
            'Cyto-aligned ETA', 'Time from cyto start (s)');

        % --- ETA mito-aligned ---
        plotEta(axEtaMito, tWin, cyMapM, miMapM, clrCyto, clrMito, ...
            'Mito-aligned ETA', 'Time from mito start (s)');

        % --- HISTOGRAMS ---
        plotHistOrBlank(axDurC, cyDur(cyDur > 0), clrCyto, ...
            'Cyto event duration', 'Duration (s)');
        plotHistOrBlank(axDurM, miDur(miDur > 0), clrMito, ...
            'Mito event duration', 'Duration (s)');
        plotHistOrBlank(axLag,  miLag(isfinite(miLag)), clrLag, ...
            'Cyto-Mito lag', 'Lag (s)');
    end


    function plotMarks(ax, x, yL, clr, ls, a)
        % Batched vertical lines drawn in a single PLOT call.
        if isempty(x), return; end
        x  = x(:);
        xx = reshape([x, x, nan(length(x), 1)]', [], 1);
        yy = repmat([yL(1); yL(2); NaN], length(x), 1);
        plot(ax, xx, yy, 'Color', [clr, a], 'LineStyle', ls, ...
            'LineWidth', 0.7, 'HandleVisibility', 'off');
    end


    function plotEta(ax, tW, mapC, mapM, clrC, clrM, ttl, xlbl)
        cla(ax, 'reset');
        if isempty(tW) || (isempty(mapC) && isempty(mapM))
            text(ax, 0.5, 0.5, '(no events for ETA)', ...
                'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                'Color', [0.5 0.5 0.5]);
            title(ax, ttl, 'FontWeight', 'normal');
            return;
        end
        hold(ax, 'on');
        if ~isempty(mapC)
            plotMeanSEM(ax, tW, mapC, clrC, 'Cyto');
        end
        if ~isempty(mapM)
            plotMeanSEM(ax, tW, mapM, clrM, 'Mito');
        end
        yL = ylim(ax);
        plot(ax, [0 0], yL, 'k:', 'LineWidth', 0.8, ...
            'HandleVisibility', 'off');
        xlabel(ax, xlbl);
        ylabel(ax, 'dF/F');
        legend(ax, 'Location', 'northeast', 'Box', 'off');
        title(ax, ttl, 'FontWeight', 'normal');
        grid(ax, 'on');
        hold(ax, 'off');
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
                'flgKDE', true, 'flgStat', false);
        else
            text(ax, 0.5, 0.5, '(no events)', 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5]);
        end
        xlabel(ax, xlbl);
        title(ax, ttl, 'FontWeight', 'normal');
    end

end     % EOF
