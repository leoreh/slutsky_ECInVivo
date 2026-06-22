function hFig = spk2ca_gui(tbl, Kd, n, Vmax)
% SPK2CA_GUI Launches the explorer for sweep results.
%
%   Controls the interactive scatter plot of Baseline vs Steady-State.
%   Kd / n / Vmax dropdowns switch the displayed variable pair. Built on the
%   shared graphics/+tblgui layer (uifigure); embeds tblGUI_scatHist.

% Figure
hFig = uifigure('Name', 'Sweep Explorer', 'Position', [100, 100, 1200, 800]);

% Layout: control strip on top, plot below
gMain = uigridlayout(hFig, [2, 1], 'RowHeight', {40, '1x'}, ...
    'Padding', [6 6 6 6], 'RowSpacing', 6);

gCtl = uigridlayout(gMain, [1, 7], ...
    'ColumnWidth', {35, 90, 25, 90, 50, 120, '1x'}, ...
    'Padding', [0 0 0 0], 'ColumnSpacing', 4);
gCtl.Layout.Row = 1;

uilabel(gCtl, 'Text', 'Kd:', 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
hPopKd = uidropdown(gCtl, 'Items', cellstr(string(Kd)), 'ItemsData', Kd, ...
    'Value', Kd(1), 'ValueChangedFcn', @update_view);

uilabel(gCtl, 'Text', 'n:', 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
hPopN = uidropdown(gCtl, 'Items', cellstr(string(n)), 'ItemsData', n, ...
    'Value', n(1), 'ValueChangedFcn', @update_view);

uilabel(gCtl, 'Text', 'Vmax:', 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
strVmax = cellfun(@v2str, Vmax, 'UniformOutput', false);
hPopVmax = uidropdown(gCtl, 'Items', strVmax, 'ItemsData', Vmax, ...
    'Value', Vmax{1}, 'ValueChangedFcn', @update_view);

% Plot panel hosts the embedded scatter GUI
plotPanel = uipanel(gMain, 'BorderType', 'none');
plotPanel.Layout.Row = 2;

% Initial variable pair (first parameter set)
[xVar, yVar] = get_varnames(Kd(1), n(1), Vmax{1});

% Launch embedded scatter GUI. Its state (setters, axes) lives on plotPanel.
tblGUI_scatHist(tbl, 'Parent', plotPanel, 'xVar', xVar, 'yVar', yVar, 'grpVar', 'Group');
hAxScat = plotPanel.UserData.hAxScatter;
update_title(Kd(1), n(1), Vmax{1});

%% ========================================================================
%  CALLBACKS
%  ========================================================================

    function update_view(~, ~)
        kV = hPopKd.Value;
        nV = hPopN.Value;
        vV = hPopVmax.Value;

        [xNew, yNew] = get_varnames(kV, nV, vV);
        plotPanel.UserData.setXYVarsFcn(xNew, yNew);
        sync_axes();
        update_title(kV, nV, vV);
    end

    function sync_axes()
        ax = plotPanel.UserData.hAxScatter;
        xl = xlim(ax); yl = ylim(ax);
        newLim = [min(xl(1), yl(1)), max(xl(2), yl(2))];
        xlim(ax, newLim);
        ylim(ax, newLim);
    end

    function [xv, yv] = get_varnames(k, n_hill, vmax_val)
        if isempty(vmax_val)
            strV = 'Auto';
        else
            strV = sprintf('%g', vmax_val);
            strV = strrep(strV, '.', 'p');
            strV = strrep(strV, '-', 'n');
        end
        suffix = sprintf('_K%g_n%g_V%s', k, n_hill, strV);
        suffix = strrep(suffix, '.', 'p');
        xv = ['mBsl' suffix];
        yv = ['mSs' suffix];
    end

    function update_title(k, n_hill, vmax_val)
        if isempty(vmax_val), vStr = 'Auto'; else, vStr = sprintf('%g', vmax_val); end
        title(hAxScat, sprintf('Sweep: Kd=%.2f, n=%g, Vmax=%s', k, n_hill, vStr));
    end

    function s = v2str(val)
        if isempty(val), s = 'Auto'; else, s = sprintf('%g', val); end
    end

end
