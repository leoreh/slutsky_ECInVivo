function hFig = guiTbl_xy(xVec, dataTbl, varargin)
% GUITBL_XY Interactive visualization of Table variables against a vector.
%
%   hFig = guiTbl_xy(xVec, dataTbl, varargin) plots column 'yVar' from 'dataTbl'
%   against 'xVec'.
%   - "Y Var": Select which variable to plot on Y-axis.
%   - "Plot By": Splits data into separate tiles (subplots).
%   - "Group By": Groups data within each tile by color.
%
%   INPUTS:
%       xVec    (numeric)  Vector for X-axis.
%       dataTbl (table)    Data table containing variables to plot.
%
%   OPTIONAL KEY-VALUE PAIRS:
%       'yVar'             (char/str) Initial variable to plot.
%       'tileFlow'         (char) 'flow', 'vertical' (1 col), 'horizontal' (1 row).
%       'tileVar'/'grpVar' (char/str) Initial tile / group variables.
%       'Parent'           (handle) uifigure or uifigure container.
%       'SelectionCallback'/'GroupByCallback' (function_handle) host coordination.
%       'xLbl'             (char/str) X-axis label.
%
%   Built on the shared graphics/gui layer (uifigure + uigridlayout).
%
%   See also: GUITBL_SCATHIST

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addParameter(p, 'yVar', [], @(x) ischar(x) || isstring(x) || isempty(x));
addParameter(p, 'tileFlow', 'vertical', @(x) any(strcmpi(x, {'flow', 'vertical', 'horizontal'})));
addParameter(p, 'tileVar', [], @(x) ischar(x) || isstring(x) || isempty(x));
addParameter(p, 'grpVar', [], @(x) ischar(x) || isstring(x) || isempty(x));
addParameter(p, 'SelectionCallback', [], @(x) isempty(x) || isa(x, 'function_handle'));
addParameter(p, 'Parent', [], @(x) isempty(x) || isgraphics(x));
addParameter(p, 'GroupByCallback', [], @(x) isempty(x) || isa(x, 'function_handle'));
addParameter(p, 'xLbl', 'Time / X', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});

initialYVar = p.Results.yVar;
initialTileVar = p.Results.tileVar;
initialGrpVar = p.Results.grpVar;
tileFlow = lower(char(p.Results.tileFlow));
hParent = p.Results.Parent;
selCbk = p.Results.SelectionCallback;
grpCbk = p.Results.GroupByCallback;
xLbl = p.Results.xLbl;

% Find all variables that match xVec dimensions (potential Y vars). This is
% specific to xy (matches a numeric matrix column or a cell of vectors to the
% length of xVec), so it stays local rather than using gui_classifyVars.
xLen = length(xVec);
yVars = {};
varNames = dataTbl.Properties.VariableNames;
for iVar = 1:length(varNames)
    raw = dataTbl.(varNames{iVar});
    if isnumeric(raw) && size(raw, 2) == xLen
        yVars{end+1} = varNames{iVar}; %#ok<AGROW>
    elseif iscell(raw) && ~isempty(raw) && length(raw{1}) == xLen
        yVars{end+1} = varNames{iVar}; %#ok<AGROW>
    end
end
if isempty(yVars)
    error('No suitable variable found in table matching xVec length.');
end

% Auto-select yVar if empty or invalid
if isempty(initialYVar) || ~ismember(initialYVar, yVars)
    yVar = yVars{1};
    if ~isempty(initialYVar)
        warning('Requested yVar "%s" not suitable. Defaulting to "%s".', initialYVar, yVar);
    end
else
    yVar = char(initialYVar);
end

%% ========================================================================
%  INITIALIZATION
%  ========================================================================

% Grouping/tiling variables (categorical / string / logical), 'None' first.
[~, catVars] = gui_classifyVars(dataTbl);
catVars = [{'None'}, catVars];

% Figure Setup
if isempty(hParent)
    hContainer = uifigure('Name', sprintf('XY Plot: %s', yVar), ...
        'Position', [100, 100, 1400, 800]);
    hFig = hContainer;
else
    hContainer = hParent;
    hFig = ancestor(hContainer, 'figure');
end

% Store GUI Data
guiData = struct();
guiData.xVec = xVec;
guiData.dataTbl = dataTbl;
guiData.yVar = yVar;
guiData.yVars = yVars;
guiData.catVars = catVars;
guiData.tileFlow = tileFlow;
guiData.chkPlotBy = gobjects(0);
guiData.chkGrpBy = gobjects(0);
guiData.tileInfo = [];
guiData.hlHandles = [];
guiData.highlightFcn = @highlightTraces;
guiData.setGroupVarFcn = @setGroupVar;
guiData.selCbk = selCbk;
guiData.grpCbk = grpCbk;
guiData.xLbl = xLbl;

%% ========================================================================
%  LAYOUT
%  ========================================================================

[~, gPlot, gCtrl] = gui_layout(hContainer, 'CtrlWidth', 190);

% Plot side: a panel hosts the tiledlayout (tiledlayout cannot parent directly
% into a uigridlayout cell).
guiData.hPanelRight = uipanel(gPlot, 'BorderType', 'none');
guiData.hLayout = tiledlayout(guiData.hPanelRight, 'flow', ...
    'TileSpacing', 'tight', 'Padding', 'compact');

% Controls
valY = find(strcmp(yVars, yVar), 1);
guiData.ddYVar = gui_labeledControl(gCtrl, 'dropdown', 'Y Variable:', ...
    'Items', yVars, 'Value', yVars{max(valY, 1)}, 'ValueChangedFcn', @onYVarChange);

guiData.ddPlotBy = gui_labeledControl(gCtrl, 'dropdown', 'Plot By (Tiles):', ...
    'Items', catVars, 'Value', pickCat(catVars, initialTileVar), ...
    'ValueChangedFcn', @onPlotByChange);
guiData.pnlPlotBy = gui_labeledControl(gCtrl, 'panel', '', 'RowHeight', '1x');

guiData.ddGrpBy = gui_labeledControl(gCtrl, 'dropdown', 'Group By (Colors):', ...
    'Items', catVars, 'Value', pickCat(catVars, initialGrpVar), ...
    'ValueChangedFcn', @onGrpByChange);
guiData.pnlGrpBy = gui_labeledControl(gCtrl, 'panel', '', 'RowHeight', '1x');

guiData.ddDispersion = gui_labeledControl(gCtrl, 'dropdown', 'Dispersion:', ...
    'Items', {'Traces', 'Spread', 'None'}, 'Value', 'Spread', 'ValueChangedFcn', @onUpdatePlot);
guiData.ddStatType = gui_labeledControl(gCtrl, 'dropdown', '', ...
    'Items', {'Arithmetic', 'Geometric', 'Median'}, 'ValueChangedFcn', @onUpdatePlot);

hContainer.UserData = guiData;

% Initial State
onPlotByChange(hContainer, []);
onGrpByChange(hContainer, []);
onUpdatePlot(hContainer, []);

%% ========================================================================
%  CALLBACKS
%  ========================================================================

    function onYVarChange(src, ~)
        data = hContainer.UserData;
        newVar = src.Value;
        if ~strcmp(data.yVar, newVar)
            data.yVar = newVar;
            hContainer.UserData = data;
            if strcmp(hContainer.Type, 'figure')
                hContainer.Name = sprintf('XY Plot: %s', newVar);
            end
            onUpdatePlot(hContainer, []);
        end
    end

    function onPlotByChange(~, ~)
        data = hContainer.UserData;
        data.chkPlotBy = populateFilter(data.ddPlotBy, data.pnlPlotBy);
        hContainer.UserData = data;
        onUpdatePlot(hContainer, []);
    end

    function onGrpByChange(~, ~)
        data = hContainer.UserData;
        data.chkGrpBy = populateFilter(data.ddGrpBy, data.pnlGrpBy);
        hContainer.UserData = data;

        if ~isempty(data.grpCbk)
            [~, allCats] = gui_selectedCats(data.chkGrpBy);
            data.grpCbk(data.ddGrpBy.Value, allCats, hContainer);
        end
        onUpdatePlot(hContainer, []);
    end

    function onFilterChange(~, ~)
        data = hContainer.UserData;
        onUpdatePlot(hContainer, []);
        if ~isempty(data.grpCbk)
            activeCats = gui_selectedCats(data.chkGrpBy);
            data.grpCbk(data.ddGrpBy.Value, activeCats, hContainer);
        end
    end

    function chk = populateFilter(dd, pnl)
        % Build the category checkboxes for a Plot By / Group By dropdown.
        varName = dd.Value;
        if strcmp(varName, 'None')
            delete(allchild(pnl));
            chk = gobjects(0);
        else
            cats = gui_catList(hContainer.UserData.dataTbl.(varName));
            chk = gui_filterPanel(pnl, cats, @onFilterChange);
        end
    end

    function onUpdatePlot(~, ~)
        data = hContainer.UserData;

        % Tiles (Plot By)
        varPB = data.ddPlotBy.Value;
        if strcmp(varPB, 'None')
            catsPB = {'All'};
        else
            catsPB = gui_selectedCats(data.chkPlotBy);
            if isempty(catsPB)
                gui_notify(hContainer, 'Select at least one "Plot By" category.', 'info');
                return;
            end
        end

        % Recreate tiled layout
        delete(data.hLayout);
        data.hLayout = tiledlayout(data.hPanelRight, data.tileFlow, ...
            'TileSpacing', 'tight', 'Padding', 'compact');
        data.tileInfo = struct('catName', {}, 'hAx', {}, 'indices', {});
        data.hlHandles = [];

        % Groups (Group By)
        varGB = data.ddGrpBy.Value;
        if strcmp(varGB, 'None')
            catsGB = {'All'};
            allCatsGB = {'All'};
        else
            [catsGB, allCatsGB] = gui_selectedCats(data.chkGrpBy);
            if isempty(catsGB)
                gui_notify(hContainer, 'Select at least one "Group By" category.', 'info');
                return;
            end
        end

        % Stable group colors over the full category list
        [fullClr, idxOf] = gui_groupColors(allCatsGB);

        if ~ismember(data.yVar, data.dataTbl.Properties.VariableNames)
            warning('Selected variable %s not in table. Resetting.', data.yVar);
            return;
        end
        yRaw = data.dataTbl.(data.yVar);
        isMatrix = isnumeric(yRaw);

        % Plot options
        dispMode = data.ddDispersion.Value;       % Traces / Spread / None
        showTraces = strcmp(dispMode, 'Traces');
        showShade  = strcmp(dispMode, 'Spread');
        method = data.ddStatType.Value;           % Arithmetic / Geometric / Median

        % Geometric floor: 1 event per recording duration (xVec assumed hr).
        totalRange = range(data.xVec);
        if totalRange == 0, totalRange = median(diff(data.xVec), 'omitnan'); end
        floorVal = 1 / (totalRange * 3600);

        axHandles = [];

        for iTile = 1:length(catsPB)
            catTile = catsPB{iTile};

            if strcmp(varPB, 'None')
                idxTile = true(height(data.dataTbl), 1);
            else
                rawCol = data.dataTbl.(varPB);
                if islogical(rawCol) || ~iscategorical(rawCol), rawCol = categorical(rawCol); end
                idxTile = (rawCol == catTile);
            end
            if sum(idxTile) == 0, continue; end

            hAx = nexttile(data.hLayout);
            data.tileInfo(end+1).catName = catTile;
            data.tileInfo(end).hAx = hAx;
            data.tileInfo(end).indices = idxTile;
            axHandles(end+1) = hAx; %#ok<AGROW>
            hold(hAx, 'on');

            tileMeanMin = inf; tileMeanMax = -inf; hasData = false;

            for iGrp = 1:length(catsGB)
                catGrp = catsGB{iGrp};

                if strcmp(varGB, 'None')
                    idxGrp = true(height(data.dataTbl), 1);
                    clr = [0, 0, 0];
                else
                    rawColG = data.dataTbl.(varGB);
                    if islogical(rawColG) || ~iscategorical(rawColG), rawColG = categorical(rawColG); end
                    idxGrp = (rawColG == catGrp);
                    clr = fullClr(idxOf(catGrp), :);
                end

                finalIdx = idxTile & idxGrp;
                if sum(finalIdx) == 0, continue; end

                if isMatrix
                    subY = yRaw(finalIdx, :);
                else
                    subY = cell2mat(yRaw(finalIdx));
                end

                % Individual traces
                if showTraces
                    hLines = plot(hAx, data.xVec, subY', 'Color', [clr, 0.05], ...
                        'LineWidth', 0.5, 'HandleVisibility', 'off');
                    globIndices = find(finalIdx);
                    for iL = 1:length(hLines)
                        hLines(iL).UserData = globIndices(iL);
                        hLines(iL).ButtonDownFcn = @(s, e) onLineClick(s, e, hContainer);
                    end
                end

                % Central tendency and bounds
                if strcmpi(method, 'Geometric')
                    [mData, lowerBound, upperBound] = gui_groupStat(subY, method, 'Floor', floorVal);
                else
                    [mData, lowerBound, upperBound] = gui_groupStat(subY, method);
                end

                % Shade (SEM / CI)
                if showShade && sum(finalIdx) > 1
                    xConf = [data.xVec(:); flipud(data.xVec(:))];
                    yConf = [upperBound(:); flipud(lowerBound(:))];
                    if all(~isnan(yConf))
                        fill(hAx, xConf, yConf, clr, 'FaceAlpha', 0.2, ...
                            'EdgeColor', 'none', 'HandleVisibility', 'off');
                    end
                end

                % Mean line
                plot(hAx, data.xVec, mData, 'Color', clr, 'LineWidth', 2, ...
                    'DisplayName', sprintf('%s (n=%d)', catGrp, sum(finalIdx)));

                tileMeanMin = min(tileMeanMin, min(mData));
                tileMeanMax = max(tileMeanMax, max(mData));
                hasData = true;
            end

            grid(hAx, 'on');
            title(hAx, catTile, 'Interpreter', 'none');
            axis(hAx, 'tight');
            if hasData && ~isinf(tileMeanMin) && ~isinf(tileMeanMax)
                yRange = tileMeanMax - tileMeanMin;
                if yRange == 0, yRange = 1; end
                ylim(hAx, [tileMeanMin - 0.1*yRange, tileMeanMax + 0.1*yRange]);
            end
            if ~strcmp(varGB, 'None')
                legend(hAx, 'Location', 'best', 'Interpreter', 'none');
            end
            hold(hAx, 'off');
        end

        xlabel(data.hLayout, data.xLbl);
        ylabel(data.hLayout, data.yVar, 'Interpreter', 'none');
        title(data.hLayout, sprintf('%s Plot By: %s | Group By: %s', data.yVar, varPB, varGB), 'Interpreter', 'none');

        % Link X axes (drawnow first: linkaxes on uiaxes needs a render pass)
        if ~isempty(axHandles)
            drawnow;
            linkaxes(axHandles, 'x');
        end

        hContainer.UserData = data;
    end

    function highlightTraces(indices)
        data = hContainer.UserData;
        if isnumeric(indices)
            tmp = false(height(data.dataTbl), 1);
            tmp(indices) = true;
            indices = tmp;
        end
        if isfield(data, 'hlHandles'), delete(data.hlHandles); end
        data.hlHandles = [];
        if ~any(indices), hContainer.UserData = data; return; end

        yRaw = data.dataTbl.(data.yVar);
        isMatrix = isnumeric(yRaw);

        for i = 1:length(data.tileInfo)
            ti = data.tileInfo(i);
            selInTile = ti.indices & indices;
            if ~any(selInTile), continue; end

            if isMatrix
                subY = yRaw(selInTile, :);
            else
                subY = cell2mat(yRaw(selInTile));
            end

            if ismember('UnitID', data.dataTbl.Properties.VariableNames)
                uid = data.dataTbl.UnitID(selInTile);
                if iscell(uid), uid = uid{1}; end
                if isnumeric(uid), uidStr = num2str(uid); else, uidStr = char(uid); end
                unitIDArg = {'DisplayName', sprintf('Unit # %s', uidStr)};
            else
                unitIDArg = {'DisplayName', 'Selected Trace'};
            end

            hold(ti.hAx, 'on');
            h = plot(ti.hAx, data.xVec, subY', 'Color', [1 1 0 0.7], 'LineWidth', 1.5, unitIDArg{:});
            data.hlHandles = [data.hlHandles; h];
            hold(ti.hAx, 'off');
            legend(ti.hAx, 'show');
        end
        hContainer.UserData = data;
    end

    function onLineClick(srcLine, ~, hFigContainer)
        try
            data = hFigContainer.UserData;
            idx = srcLine.UserData;
            highlightTraces(idx);
            if ~isempty(data.selCbk)
                inPoints = false(height(data.dataTbl), 1);
                inPoints(idx) = true;
                data.selCbk(inPoints);
            end
        catch ME
            warning(ME.identifier, '%s', ME.message);
        end
    end

    function setGroupVar(varName, activeCats)
        data = hContainer.UserData;
        if ~ismember(varName, data.ddGrpBy.Items), return; end

        needsUpdate = false;
        if ~strcmp(data.ddGrpBy.Value, varName)
            data.ddGrpBy.Value = varName;
            data.chkGrpBy = populateFilter(data.ddGrpBy, data.pnlGrpBy);
            hContainer.UserData = data;
            needsUpdate = true;
        end

        if nargin > 1 && ~isempty(data.chkGrpBy)
            for i = 1:numel(data.chkGrpBy)
                val = ismember(data.chkGrpBy(i).Text, activeCats);
                if data.chkGrpBy(i).Value ~= val
                    data.chkGrpBy(i).Value = val;
                    needsUpdate = true;
                end
            end
        end

        if needsUpdate, onUpdatePlot(hContainer, []); end
    end

end

function val = pickCat(catVars, requested)
% Resolve an initial dropdown value: the requested category if present, else
% 'None'.
val = 'None';
if ~isempty(requested) && ismember(char(requested), catVars)
    val = char(requested);
end
end
