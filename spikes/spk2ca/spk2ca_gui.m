function spk2ca_gui(tbl, Kd, n, Vmax)
% SPK2CA_GUI Launches the explorer for sweep results.
%
%   Controls the interactive scatter plot of Baseline vs Steady-State.

% Create Figure
hFig = figure('Name', 'Sweep Explorer', 'Color', 'w', ...
    'Position', [100, 100, 1200, 800], 'NumberTitle', 'off');

% Layout: Top Controls, Bottom Plot
hPnlCtl = uipanel('Parent', hFig, 'Position', [0 0.94 1 0.06], ...
    'BorderType', 'none');
hPnlPlot = uipanel('Parent', hFig, 'Position', [0 0 1 0.94], ...
    'BorderType', 'none');

% Control: Kd
uicontrol('Parent', hPnlCtl, 'Style', 'text', 'String', 'Kd:', ...
    'Units', 'normalized', 'Position', [0.02 0.25 0.03 0.5], ...
    'HorizontalAlignment', 'right', ...
    'FontWeight', 'bold');
hPopKd = uicontrol('Parent', hPnlCtl, 'Style', 'popupmenu', ...
    'String', string(Kd), 'Units', 'normalized', ...
    'Position', [0.06 0.25 0.06 0.5], ...
    'Callback', @update_view);

% Control: n
uicontrol('Parent', hPnlCtl, 'Style', 'text', 'String', 'n:', ...
    'Units', 'normalized', 'Position', [0.14 0.25 0.02 0.5], ...
    'HorizontalAlignment', 'right', ...
    'FontWeight', 'bold');
hPopN = uicontrol('Parent', hPnlCtl, 'Style', 'popupmenu', ...
    'String', string(n), 'Units', 'normalized', ...
    'Position', [0.17 0.25 0.06 0.5], ...
    'Callback', @update_view);

% Control: Vmax
uicontrol('Parent', hPnlCtl, 'Style', 'text', 'String', 'Vmax:', ...
    'Units', 'normalized', 'Position', [0.25 0.25 0.04 0.5], ...
    'HorizontalAlignment', 'right', ...
    'FontWeight', 'bold');

% Create Vmax labels for dropdown
strVmax = cellfun(@(x) v2str(x), Vmax, 'UniformOutput', false);

hPopVmax = uicontrol('Parent', hPnlCtl, 'Style', 'popupmenu', ...
    'String', strVmax, 'Units', 'normalized', ...
    'Position', [0.30 0.25 0.10 0.5], ...
    'Callback', @update_view);


% Initial Vars (First parameter set)
kVal = Kd(1);
nVal = n(1);
vVal = Vmax{1};
[xVar, yVar] = get_varnames(kVal, nVal, vVal);

% Launch GUI Once
hG = tblGUI_scatHist(tbl, 'Parent', hPnlPlot, ...
    'xVar', xVar, 'yVar', yVar, ...
    'grpVar', 'Group');

% Store handle to update title later
hAx = findobj(hG, 'Type', 'axes', 'Tag', '');
update_title(kVal, nVal, vVal);


    function update_view(~, ~)

        % Get Indices
        idxKd = hPopKd.Value;
        idxN  = hPopN.Value;
        idxV  = hPopVmax.Value;

        % Check bounds
        if idxKd > length(Kd) || idxN > length(n) || idxV > length(Vmax), return; end

        kV = Kd(idxKd);
        nV = n(idxN);
        vV = Vmax{idxV};

        % Construct Variable Names
        [xNew, yNew] = get_varnames(kV, nV, vV);

        % Update Scatter Plot via Callback
        hG.UserData.setXYVarsFcn(xNew, yNew);

        % Sync Axes
        sync_axes();

        % Update Title
        update_title(kV, nV, vV);
    end

    function sync_axes()
        data = hG.UserData;
        ax = data.hAxScatter;

        xl = xlim(ax);
        yl = ylim(ax);

        minVal = min(xl(1), yl(1));
        maxVal = max(xl(2), yl(2));

        newLim = [minVal, maxVal];
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
        if ~isempty(hAx)
            if isempty(vmax_val)
                vStr = 'Auto';
            else
                vStr = sprintf('%g', vmax_val);
            end
            title(hAx, sprintf('Sweep: Kd=%.2f, n=%g, Vmax=%s', k, n_hill, vStr));
        end
    end

    function s = v2str(val)
        if isempty(val)
            s = 'Auto';
        else
            s = sprintf('%g', val);
        end
    end

end
