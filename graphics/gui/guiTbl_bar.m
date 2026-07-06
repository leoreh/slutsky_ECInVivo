function hFig = guiTbl_bar(tbl, varargin)
% GUITBL_BAR Interactive bar plot with grouping and error bars.
%
%   hFig = guiTbl_bar(tbl, ...) opens a GUI to visualize the table 'tbl'
%   as a bar chart.
%
%   INPUT:
%       tbl         (table) The data table to visualize.
%
%   OPTIONAL KEY-VALUE PAIRS:
%       'yVar'        (string) Initial Y variable name (Numeric)
%       'xVar'        (string) Initial X variable name (Categorical)
%       'grpVar'      (string) Initial Group variable name (Categorical)
%       'mode'        (string) Initial view mode: 'bar' (default) or
%                     'points'. Toggle at runtime via the Mode button.
%       'Parent'      (handle) Parent container (a uifigure or a uifigure
%                     container such as a uipanel / uigridlayout cell).
%
%   Built on the shared graphics/gui layer (uifigure + uigridlayout).
%
%   See also: GUITBL_SCATHIST, GUITBL_XY

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addParameter(p, 'yVar', '', @(x) ischar(x) || isstring(x) || isempty(x));
addParameter(p, 'xVar', '', @(x) ischar(x) || isstring(x) || isempty(x));
addParameter(p, 'grpVar', '', @(x) ischar(x) || isstring(x) || isempty(x));
addParameter(p, 'mode', 'bar', @(x) ischar(x) || isstring(x));
addParameter(p, 'Parent', [], @(x) isempty(x) || isgraphics(x));
parse(p, tbl, varargin{:});

tbl = p.Results.tbl;
yVarIn = p.Results.yVar;
xVarIn = p.Results.xVar;
grpVarIn = p.Results.grpVar;
hParent = p.Results.Parent;

modeIn = lower(char(string(p.Results.mode)));
if ~ismember(modeIn, {'bar', 'points'})
    warning('guiTbl_bar:invalidMode', ...
        'Invalid mode ''%s''; falling back to ''bar''.', modeIn);
    modeIn = 'bar';
end

%% ========================================================================
%  INITIALIZATION
%  ========================================================================

% Numeric (Y) and categorical (X / Group) variables, sorted for the menus.
% 'vector' keeps only scalar-per-row numerics (bars aggregate one value/row).
[numericVars, catVars] = gui_classifyVars(tbl, 'SortNames', true, 'Shape', 'vector');

if isempty(numericVars)
    error('Input table must contain at least one numeric variable.');
end
if isempty(catVars)
    error('Input table must contain at least one categorical variable.');
end

% Defaults
curY = numericVars{1};
curX = catVars{1};
curGrp = 'None';

% Apply overrides
if ~isempty(yVarIn) && ismember(yVarIn, numericVars), curY = char(yVarIn); end
if ~isempty(xVarIn) && ismember(xVarIn, catVars), curX = char(xVarIn); end
if ~isempty(grpVarIn) && ismember(grpVarIn, catVars), curGrp = char(grpVarIn); end

% Group dropdown also offers 'None'
grpVars = [{'None'}, catVars];

% Figure / parent
if isempty(hParent)
    hContainer = uifigure('Name', 'Bar Plot Visualizer', 'Position', [100, 100, 1000, 700]);
    hFig = hContainer;
else
    hContainer = hParent;
    hFig = ancestor(hContainer, 'figure');
end

% GUI Data
guiData = struct();
guiData.tbl = tbl;
guiData.numericVars = numericVars;
guiData.catVars = catVars;
guiData.grpVars = grpVars;
guiData.chkGrp = gobjects(0);

%% ========================================================================
%  LAYOUT
%  ========================================================================

[~, gPlot, gCtrl] = gui_layout(hContainer, 'CtrlWidth', 210);
guiData.hAx = uiaxes(gPlot);

% Controls (stacked in the scrollable control column)
guiData.ddY = gui_labeledControl(gCtrl, 'dropdown', 'Y Variable (Numeric):', ...
    'Items', numericVars, 'Value', curY, 'ValueChangedFcn', @onUpdatePlot);
guiData.ddX = gui_labeledControl(gCtrl, 'dropdown', 'X Variable (Cat):', ...
    'Items', catVars, 'Value', curX, 'ValueChangedFcn', @onUpdatePlot);
guiData.ddGrp = gui_labeledControl(gCtrl, 'dropdown', 'Group Variable (Cat):', ...
    'Items', grpVars, 'Value', curGrp, 'ValueChangedFcn', @onGrpChange);
guiData.ddStatType = gui_labeledControl(gCtrl, 'dropdown', 'Statistic:', ...
    'Items', {'Arithmetic', 'Geometric', 'Median'}, 'ValueChangedFcn', @onUpdatePlot);

modeIsPoints = strcmp(modeIn, 'points');
if modeIsPoints, tgLabel = 'Mode: Points'; else, tgLabel = 'Mode: Bars'; end
guiData.tgMode = gui_labeledControl(gCtrl, 'toggle', '', ...
    'Text', tgLabel, 'Value', modeIsPoints, 'ValueChangedFcn', @onModeToggle);

gui_labeledControl(gCtrl, 'label', 'Groups:');
guiData.pnlGrp = gui_labeledControl(gCtrl, 'panel', '', 'RowHeight', '1x');

hContainer.UserData = guiData;

% Initial population & plot
onGrpChange(hContainer, []);

%% ========================================================================
%  CALLBACKS
%  ========================================================================

    function onUpdatePlot(~, ~)
        data = hContainer.UserData;
        tblIn = data.tbl;

        yName = data.ddY.Value;
        xName = data.ddX.Value;
        grpName = data.ddGrp.Value;

        % Data Prep
        yData = tblIn.(yName);
        xData = tblIn.(xName);
        if ~iscategorical(xData), xData = categorical(xData); end
        xCats = categories(xData);
        % Only keep categories actually used (for non-NaN responses)
        xCats = xCats(ismember(xCats, unique(xData(~isnan(yData)))));

        hasGrp = ~strcmp(grpName, 'None');
        if hasGrp
            grpData = tblIn.(grpName);
            if ~iscategorical(grpData), grpData = categorical(grpData); end
            allGcats = gui_catList(grpData);     % full (unfiltered) order
            if ~isempty(data.chkGrp)
                active = gui_selectedCats(data.chkGrp);
                gCats = intersect(allGcats, active, 'stable');
            else
                gCats = allGcats;
            end
        else
            grpData = [];
            allGcats = {'All'};
            gCats = {'All'};
        end

        statType = data.ddStatType.Value;

        % Per-group colors, stable by full category list
        [fullClr, idxOf] = gui_groupColors(allGcats);
        clrMat = zeros(max(numel(gCats), 1), 3);
        for iG = 1:numel(gCats)
            clrMat(iG, :) = fullClr(idxOf(gCats{iG}), :);
        end

        % Aggregate Data (bar expects matrix nXCats x nGroups)
        meanMat = nan(length(xCats), length(gCats));
        errLMat = nan(length(xCats), length(gCats));
        errHMat = nan(length(xCats), length(gCats));
        valsByCell = cell(length(xCats), length(gCats));

        for iX = 1:length(xCats)
            for iG = 1:length(gCats)
                if hasGrp
                    idx = (xData == xCats{iX}) & (grpData == gCats{iG});
                else
                    idx = (xData == xCats{iX});
                end

                vals = yData(idx);
                vals = vals(~isnan(vals));

                if ~isempty(vals)
                    valsByCell{iX, iG} = vals;
                    [m, lo, hi] = gui_groupStat(vals, statType);
                    meanMat(iX, iG) = m;
                    errLMat(iX, iG) = m - lo;
                    errHMat(iX, iG) = hi - m;
                end
            end
        end

        % Plotting
        ax = data.hAx;
        % cla() skips children with HandleVisibility='off' (used below for
        % legend control), so use delete(allchild(...)) to force-clear all.
        delete(allchild(ax));
        hold(ax, 'on');

        % Common grouped-bar geometry (also used to place jittered points)
        nbars = length(gCats);
        groupwidth = min(0.8, nbars/(nbars + 1.5));

        useScatter = logical(data.tgMode.Value);
        if useScatter
            [hLegend, legendCats] = renderScatter(ax, valsByCell, meanMat, ...
                errLMat, errHMat, clrMat, groupwidth, gCats);
        else
            hLegend = renderBars(ax, meanMat, errLMat, errHMat, clrMat, groupwidth);
            legendCats = gCats;
        end

        hold(ax, 'off');

        % Aesthetics
        set(ax, 'XTick', 1:length(xCats), 'XTickLabel', xCats);
        ylabel(ax, yName, 'Interpreter', 'none');
        xlabel(ax, xName, 'Interpreter', 'none');

        if hasGrp
            if ~isempty(hLegend)
                legend(ax, hLegend, legendCats, 'Location', 'best', 'Interpreter', 'none');
            else
                legend(ax, 'off');
            end
            title(ax, sprintf('%s by %s (grouped by %s)', yName, xName, grpName), 'Interpreter', 'none');
        else
            legend(ax, 'off');
            title(ax, sprintf('%s by %s', yName, xName), 'Interpreter', 'none');
        end

        grid(ax, 'on');
    end

    function hBars = renderBars(ax, meanMat, errLMat, errHMat, clrMat, groupwidth)
        % Grouped bars with overlaid error bars.
        % Returns the bar handles for legend assembly.

        hBars = bar(ax, meanMat, 'grouped');
        nbars = size(meanMat, 2);
        ngroups = size(meanMat, 1);

        for iG = 1:nbars
            hBars(iG).FaceColor = clrMat(iG, :);

            % Center of each bar in this group (matches bar's internal layout)
            xCent = (1:ngroups) - groupwidth/2 + (2*iG-1) * groupwidth / (2*nbars);
            errorbar(ax, xCent, meanMat(:,iG), errLMat(:,iG), errHMat(:,iG), ...
                'k', 'LineStyle', 'none', 'HandleVisibility', 'off');
        end
    end

    function [hLegend, legendCats] = renderScatter(ax, valsByCell, meanMat, ...
            errLMat, errHMat, clrMat, groupwidth, gCats)
        % Jittered individual points per (xCat, gCat), with the central
        % tendency and error bars from the active Statistic overlaid as a
        % filled diamond. One legend handle per non-empty group.

        [nX, nG] = size(valsByCell);
        nbars = nG;
        barW = groupwidth / nbars;
        jitterHalf = 0.35 * barW;

        hLegend = gobjects(0);
        legendCats = {};

        for iG = 1:nG
            grpClr = clrMat(iG, :);
            grpHasData = false;
            grpLegendHandle = gobjects(1);

            for iX = 1:nX
                vals = valsByCell{iX, iG};
                if isempty(vals), continue; end

                % Bar-center x-position (same formula as renderBars)
                xCent = iX - groupwidth/2 + (2*iG-1) * groupwidth / (2*nbars);

                % Uniform jitter in [-jitterHalf, +jitterHalf]
                nv = numel(vals);
                jit = (rand(nv, 1) - 0.5) * 2 * jitterHalf;
                xJit = xCent + jit;

                % First non-empty cell per group contributes the legend handle;
                % the rest are hidden from the legend.
                if grpHasData
                    visFlag = 'off';
                else
                    visFlag = 'on';
                end

                hS = scatter(ax, xJit, vals, 20, grpClr, 'filled', ...
                    'MarkerFaceAlpha', 0.35, 'MarkerEdgeColor', 'none', ...
                    'HandleVisibility', visFlag);

                if ~grpHasData
                    grpLegendHandle = hS;
                    grpHasData = true;
                end

                % Summary marker (central tendency)
                if ~isnan(meanMat(iX, iG))
                    plot(ax, xCent, meanMat(iX, iG), 'kd', ...
                        'MarkerFaceColor', grpClr, 'MarkerSize', 8, ...
                        'HandleVisibility', 'off');
                    errorbar(ax, xCent, meanMat(iX, iG), ...
                        errLMat(iX, iG), errHMat(iX, iG), 'k', ...
                        'LineStyle', 'none', 'HandleVisibility', 'off');
                end
            end

            if grpHasData
                hLegend(end+1) = grpLegendHandle; %#ok<AGROW>
                legendCats{end+1} = gCats{iG};   %#ok<AGROW>
            end
        end
    end

    function onGrpChange(~, ~)
        data = hContainer.UserData;
        grpName = data.ddGrp.Value;
        if strcmp(grpName, 'None')
            delete(allchild(data.pnlGrp));
            data.chkGrp = gobjects(0);
        else
            cats = gui_catList(data.tbl.(grpName));
            data.chkGrp = gui_filterPanel(data.pnlGrp, cats, @onFilterChange);
        end
        hContainer.UserData = data;
        onUpdatePlot();
    end

    function onFilterChange(~, ~)
        onUpdatePlot();
    end

    function onModeToggle(src, ~)
        if src.Value
            src.Text = 'Mode: Points';
        else
            src.Text = 'Mode: Bars';
        end
        onUpdatePlot();
    end

end
