function hFig = tblGUI_xy(xVec, dataTbl, varargin)
% TBLGUI_XY Interactive visualization of Table variables against a vector.
%
%   hFig = tblGUI_xy(xVec, dataTbl, varargin) plots column 'yVar' from 'dataTbl'
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
%                          If empty, auto-selects first suitable variable.
%       'tileFlow'         (char) Layout of tiles: 'flow',
%                          'vertical' (1 column), or 'horizontal' (1 row).
%
%   OUTPUT:
%       hFig    (handle)   Figure handle.
%
%   See also: TBLGUI_SCATHIST

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addParameter(p, 'yVar', [], @(x) ischar(x) || isstring(x) || isempty(x));
addParameter(p, 'tileFlow', 'vertical', @(x) permember(x, {'flow', 'vertical', 'horizontal'}));
addParameter(p, 'SelectionCallback', [], @(x) isempty(x) || isa(x, 'function_handle'));
addParameter(p, 'Parent', [], @(x) isempty(x) || isgraphics(x));
addParameter(p, 'GroupByCallback', [], @(x) isempty(x) || isa(x, 'function_handle'));
parse(p, varargin{:});

initialYVar = p.Results.yVar;
tileFlow = p.Results.tileFlow;
hParent = p.Results.Parent;
selCbk = p.Results.SelectionCallback;
grpCbk = p.Results.GroupByCallback;

% Find all variables that match xVec dimensions (Potential Y-Vars)
xLen = length(xVec);
yVars = {};
varNames = dataTbl.Properties.VariableNames;

for iVar = 1:length(varNames)
    raw = dataTbl.(varNames{iVar});
    if isnumeric(raw) && size(raw, 2) == xLen
        yVars{end+1} = varNames{iVar}; %#ok<AGROW>
    elseif iscell(raw) && length(raw) > 0 && length(raw{1}) == xLen
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
    yVar = initialYVar;
end


%% ========================================================================
%  INITIALIZATION
%  ========================================================================

% Identify potential grouping variables (Categorical/String/Logical)
allVars = dataTbl.Properties.VariableNames;
catVars = allVars(varfun(@(x) iscategorical(x) || isstring(x) || islogical(x),...
    dataTbl, 'OutputFormat', 'uniform'));
catVars = [{'None'}, catVars];

% Figure Setup
if isempty(hParent)
    hContainer = figure('Name', sprintf('XY Plot: %s', yVar), 'NumberTitle', 'off', ...
        'Units', 'pixels', 'Position', [100, 100, 1400, 800]);
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
guiData.yVars = yVars; % List of available Y variables
guiData.catVars = catVars;
guiData.tileFlow = tileFlow;
guiData.chkPlotBy = [];
guiData.chkGrpBy = [];
guiData.colors = lines(20);
guiData.tileInfo = [];
guiData.hlHandles = [];
guiData.highlightFcn = @highlightTraces;
guiData.setGroupVarFcn = @setGroupVar;
guiData.selCbk = selCbk;
guiData.grpCbk = grpCbk;


%% ========================================================================
%  LAYOUT
%  ========================================================================

% Side Panel for Controls (Left)
panelW = 0.1;
hPanel = uipanel('Parent', hContainer, 'Units', 'normalized', ...
    'Position', [0, 0, panelW, 1]);

% Main Plotting Area (Right)
hPanelRight = uipanel('Parent', hContainer, 'Units', 'normalized', ...
    'Position', [panelW, 0, 1-panelW, 1], 'BorderType', 'none');

guiData.hLayout = tiledlayout(hPanelRight, 'flow', ...
    'TileSpacing', 'tight', 'Padding', 'compact');
guiData.hPanelRight = hPanelRight;

% --- Controls in Side Panel ---
ctlH = 0.03;
ctlGap = 0.01;
currY = 0.95;

% 1. Y VARIABLE SELECTION
uicontrol('Parent', hPanel, 'Style', 'text', 'String', 'Y Variable:', ...
    'Units', 'normalized', 'Position', [0.05, currY, 0.9, ctlH], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');
currY = currY - ctlH;

valY = find(strcmp(yVars, yVar), 1);
if isempty(valY), valY = 1; end

guiData.ddYVar = uicontrol('Parent', hPanel, 'Style', 'popupmenu', ...
    'String', yVars, 'Value', valY, ...
    'Units', 'normalized', 'Position', [0.05, currY, 0.9, ctlH], ...
    'Callback', @onYVarChange);
currY = currY - ctlH - ctlGap*2; % Extra gap


% 2. PLOT BY (Tiles)
uicontrol('Parent', hPanel, 'Style', 'text', 'String', 'Plot By (Tiles):', ...
    'Units', 'normalized', 'Position', [0.05, currY, 0.9, ctlH], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');
currY = currY - ctlH;

guiData.ddPlotBy = uicontrol('Parent', hPanel, 'Style', 'popupmenu', ...
    'String', catVars, 'Units', 'normalized', ...
    'Position', [0.05, currY, 0.9, ctlH], ...
    'Callback', @onPlotByChange);
currY = currY - ctlH - ctlGap;

% Container for Plot By Checkboxes
guiData.pnlPlotBy = uipanel('Parent', hPanel, 'BorderType', 'none', ...
    'Units', 'normalized', 'Position', [0.05, 0.55, 0.9, currY - 0.55]);

% Separator
currY = 0.50;
uicontrol('Parent', hPanel, 'Style', 'text', 'String', '-----------------------------', ...
    'Units', 'normalized', 'Position', [0.05, currY, 0.9, ctlH], ...
    'HorizontalAlignment', 'center', 'ForegroundColor', [0.5 0.5 0.5]);
currY = currY - ctlH - ctlGap;

% 3. GROUP BY (Colors)
uicontrol('Parent', hPanel, 'Style', 'text', 'String', 'Group By (Colors):', ...
    'Units', 'normalized', 'Position', [0.05, currY, 0.9, ctlH], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');
currY = currY - ctlH;

guiData.ddGrpBy = uicontrol('Parent', hPanel, 'Style', 'popupmenu', ...
    'String', catVars, 'Units', 'normalized', ...
    'Position', [0.05, currY, 0.9, ctlH], ...
    'Callback', @onGrpByChange);
currY = currY - ctlH - ctlGap;

% Container for Group By Checkboxes
guiData.pnlGrpBy = uipanel('Parent', hPanel, 'BorderType', 'none', ...
    'Units', 'normalized', 'Position', [0.05, 0.10, 0.9, currY - 0.10]);


hContainer.UserData = guiData;
hFig = hContainer;

% Initial State
onPlotByChange(hFig, []);
onGrpByChange(hFig, []);

% Initial Plot
onUpdatePlot(hContainer, []);

%% ========================================================================
%  CALLBACKS
%  ========================================================================

    function onYVarChange(src, ~)
        data = hContainer.UserData;
        idx = get(src, 'Value');
        items = get(src, 'String');

        newVar = items{idx};

        if ~strcmp(data.yVar, newVar)
            data.yVar = newVar;
            hContainer.UserData = data;

            % Update Figure Name if standalone
            if strcmp(hContainer.Type, 'figure')
                set(hContainer, 'Name', sprintf('XY Plot: %s', newVar));
            end

            onUpdatePlot(hContainer, []);
        end
    end

    function onPlotByChange(src, ~)
        data = hContainer.UserData;
        populateCheckboxes(src, data, data.ddPlotBy, data.pnlPlotBy, 'chkPlotBy');
        onUpdatePlot(src, []);
    end

    function onGrpByChange(src, ~)
        data = hContainer.UserData;
        populateCheckboxes(src, data, data.ddGrpBy, data.pnlGrpBy, 'chkGrpBy');

        idx = get(data.ddGrpBy, 'Value');
        items = get(data.ddGrpBy, 'String');
        varName = items{idx};

        if ~isempty(data.grpCbk)
            % Gather all cats
            populateCheckboxes(src, data, data.ddGrpBy, data.pnlGrpBy, 'chkGrpBy'); % Refresh data
            data = hContainer.UserData; % Reload

            if isempty(data.chkGrpBy)
                activeCats = {};
            else
                allCats = arrayfun(@(x) string(get(x, 'String')), data.chkGrpBy);
                activeCats = cellstr(allCats);
            end

            data.grpCbk(varName, activeCats, hContainer);
        end

        onUpdatePlot(src, []);
    end

    function onFilterChange(~, ~)
        data = hContainer.UserData;
        idx = get(data.ddGrpBy, 'Value');
        items = get(data.ddGrpBy, 'String');
        varName = items{idx};

        onUpdatePlot(hContainer, []); % Update Local

        % Broadcast
        if ~isempty(data.grpCbk)
            if isempty(data.chkGrpBy)
                activeCats = {};
            else
                selectedIdx = arrayfun(@(x) get(x, 'Value'), data.chkGrpBy);
                allCats = arrayfun(@(x) string(get(x, 'String')), data.chkGrpBy);
                activeCats = cellstr(allCats(logical(selectedIdx)));
            end
            data.grpCbk(varName, activeCats, hContainer);
        end
    end

    function populateCheckboxes(src, data, hDD, hPanel, storeField)
        % Clear existing
        delete(hPanel.Children);

        idx = get(hDD, 'Value');
        varName = data.catVars{idx};

        chkHandles = [];

        if strcmp(varName, 'None')
            % No checkboxes needed
        else
            % Get Uniques
            raw = data.dataTbl.(varName);
            if islogical(raw), raw = categorical(raw); end
            if ~iscategorical(raw), raw = categorical(raw); end
            cats = categories(raw);
            cats = cats(ismember(cats, unique(raw))); % Only present categories

            nCats = length(cats);
            h = 1 / max(10, nCats + 1); % Simple spacing

            w = 0.9;
            chkHandles = gobjects(1, nCats); % Initialize
            for kVar = 1:nCats
                yPos = 1 - kVar*h;
                chkHandles(kVar) = uicontrol('Parent', hPanel, 'Style', 'checkbox', ...
                    'String', cats{kVar}, 'Units', 'normalized', ...
                    'Position', [0, yPos, w, h], 'Value', 1, ...
                    'Callback', @onFilterChange); % Auto-update & Sync
            end
        end

        data.(storeField) = chkHandles;
        hContainer.UserData = data;
    end

    function onUpdatePlot(~, ~)
        data = hContainer.UserData;

        % Determine Tile Splits (Plot By)
        idxPB = get(data.ddPlotBy, 'Value');
        varPB = data.catVars{idxPB};

        catsPB = {'All'};
        if ~strcmp(varPB, 'None')
            % Get selected categories from checkboxes
            hChk = data.chkPlotBy;
            if isempty(hChk)
                catsPB = {}; % Should not happen if non-None
            else
                selectedIdx = arrayfun(@(x) get(x, 'Value'), hChk);
                if sum(selectedIdx) == 0
                    msgbox('Please select at least one "Plot By" category.', 'Info');
                    return;
                end
                allCats = arrayfun(@(x) string(get(x, 'String')), hChk);
                catsPB = cellstr(allCats(logical(selectedIdx)));
            end
        end
        nTiles = length(catsPB);

        % Recreate layout
        delete(data.hLayout);
        data.hLayout = tiledlayout(data.hPanelRight, data.tileFlow, ...
            'TileSpacing', 'tight', 'Padding', 'compact');
        hContainer.UserData = data;

        % Reset Tile Info
        data.tileInfo = struct('catName', {}, 'hAx', {}, 'indices', {});
        data.hlHandles = [];

        % Determine Line Groups (Group By)
        idxGB = get(data.ddGrpBy, 'Value');
        varGB = data.catVars{idxGB};

        catsGB = {'All'};
        if ~strcmp(varGB, 'None')
            % Get selected categories from checkboxes
            hChk = data.chkGrpBy;
            if isempty(hChk)
                catsGB = {};
            else
                selectedIdx = arrayfun(@(x) get(x, 'Value'), hChk);
                if sum(selectedIdx) == 0
                    msgbox('Please select at least one "Group By" category.', 'Info');
                    return;
                end
                allCats = arrayfun(@(x) string(get(x, 'String')), hChk);
                catsGB = cellstr(allCats(logical(selectedIdx)));
            end
        end

        % Retrieve Data (Dynamic based on Selection)
        if ~ismember(data.yVar, data.dataTbl.Properties.VariableNames)
            % Fallback if var missing
            warning('Selected variable %s not in table. Resetting.', data.yVar);
            return;
        end
        yRaw = data.dataTbl.(data.yVar);
        isMatrix = isnumeric(yRaw);

        axHandles = []; % Store axes for linking

        % --- RENDER LOOP ---
        for iTile = 1:length(catsPB)
            catTile = catsPB{iTile};

            % Filter Data for Tile
            if strcmp(varPB, 'None')
                idxTile = true(height(data.dataTbl), 1);
            else
                rawCol = data.dataTbl.(varPB);
                if islogical(rawCol), rawCol = categorical(rawCol); end
                if ~iscategorical(rawCol), rawCol = categorical(rawCol); end
                idxTile = (rawCol == catTile);
            end

            if sum(idxTile) == 0, continue; end

            nexttile(data.hLayout);
            hAx = gca;

            % Store Tile Info
            data.tileInfo(end+1).catName = catTile;
            data.tileInfo(end).hAx = hAx;
            data.tileInfo(end).indices = idxTile;

            axHandles(end+1) = hAx; %#ok<AGROW>
            hold(hAx, 'on');

            tileMeanMin = inf;
            tileMeanMax = -inf;
            hasData = false;

            % Iterate groups within tile
            for iGrp = 1:length(catsGB)
                catGrp = catsGB{iGrp};

                % Filter Data for Group AND Tile
                if strcmp(varGB, 'None')
                    idxGrp = true(height(data.dataTbl), 1);
                    clr = [0, 0, 0];
                else
                    rawColG = data.dataTbl.(varGB);
                    if islogical(rawColG), rawColG = categorical(rawColG); end
                    if ~iscategorical(rawColG), rawColG = categorical(rawColG); end
                    idxGrp = (rawColG == catGrp);

                    % Cycle colors
                    cIdx = mod(iGrp-1, size(data.colors,1)) + 1;
                    clr = data.colors(cIdx, :);
                end

                finalIdx = idxTile & idxGrp;
                if sum(finalIdx) == 0, continue; end

                % Extract Data
                if isMatrix
                    subY = yRaw(finalIdx, :);
                else
                    subY = cell2mat(yRaw(finalIdx));
                end

                % Plot
                % Light individual lines (High transparency) 0.05
                hLines = plot(hAx, data.xVec, subY', 'Color', [clr, 0.05],...
                    'LineWidth', 0.5, 'HandleVisibility', 'off');

                % Assign global indices to UserData of lines
                globIndices = find(finalIdx);
                for iL = 1:length(hLines)
                    hLines(iL).UserData = globIndices(iL);
                    hLines(iL).ButtonDownFcn = @(s,e) onLineClick(s, e, hContainer);
                end

                % Bold mean line
                mfr = mean(subY, 1, 'omitnan');
                plot(hAx, data.xVec, mfr, 'Color', clr, 'LineWidth', 2, ...
                    'DisplayName', sprintf('%s (n=%d)', catGrp, sum(finalIdx)));

                % Update Ranges for Y Lim
                tileMeanMin = min(tileMeanMin, min(mfr));
                tileMeanMax = max(tileMeanMax, max(mfr));
                hasData = true;
            end

            % Aesthetics
            grid(hAx, 'on');
            title(hAx, catTile, 'Interpreter', 'none');
            axis(hAx, 'tight');

            % Apply Y-Limits based on Means
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

        xlabel(data.hLayout, 'Time / X');
        ylabel(data.hLayout, data.yVar, 'Interpreter', 'none');
        title(data.hLayout, sprintf('%s Plot By: %s | Group By: %s', data.yVar, varPB, varGB), 'Interpreter', 'none');

        % Link X axes
        if ~isempty(axHandles)
            linkaxes(axHandles, 'x');
        end

        % Save updated data with tile info
        hContainer.UserData = data;

    end

    function highlightTraces(indices)
        % Highlights specific traces (logical or linear indices)
        data = hContainer.UserData;

        % convert to logical if numeric
        if isnumeric(indices)
            tmp = false(height(data.dataTbl), 1);
            tmp(indices) = true;
            indices = tmp;
        end

        % Clear existing highlights
        if isfield(data, 'hlHandles')
            delete(data.hlHandles);
        end
        data.hlHandles = [];

        if ~any(indices), return; end

        yRaw = data.dataTbl.(data.yVar);
        isMatrix = isnumeric(yRaw);

        % Loop through tiles
        for i = 1:length(data.tileInfo)
            ti = data.tileInfo(i);

            % Intersection of tile indices and selected indices
            selInTile = ti.indices & indices;

            if ~any(selInTile), continue; end

            % Get Data
            if isMatrix
                subY = yRaw(selInTile, :);
            else
                subY = cell2mat(yRaw(selInTile));
            end

            % Get UnitID if available
            unitIDArg = {};
            if ismember('UnitID', data.dataTbl.Properties.VariableNames)
                uid = data.dataTbl.UnitID(selInTile);
                if iscell(uid), uid = uid{1}; end
                if isnumeric(uid), uidStr = num2str(uid); else, uidStr = char(uid); end
                unitIDArg = {'DisplayName', sprintf('Unit # %s', uidStr)};
            else
                unitIDArg = {'DisplayName', 'Selected Trace'};
            end

            % Plot Highlight
            hold(ti.hAx, 'on');
            % yellow bolder line for selection
            h = plot(ti.hAx, data.xVec, subY', 'Color', [1 1 0 0.7], 'LineWidth', 1.5, unitIDArg{:});
            data.hlHandles = [data.hlHandles; h];
            hold(ti.hAx, 'off');

            legend(ti.hAx, 'show'); % Ensure legend updates
        end

        hContainer.UserData = data;
    end

    function onLineClick(srcLine, ~, hFigContainer)
        % Callback for line click
        % Highlighting internally + Trigger external callback

        try
            data = hFigContainer.UserData;
            idx = srcLine.UserData;

            % 1. Internal Highlight
            highlightTraces(idx);

            % 2. External Callback
            if ~isempty(data.selCbk)
                % Pass logical array or index? Standardizing on logical usually safer
                % but index is fine too. Let's send logical to be consistent.
                inPoints = false(height(data.dataTbl), 1);
                inPoints(idx) = true;
                data.selCbk(inPoints);
            end

            fprintf('Selected Trace Index: %d\n', idx);

        catch ME
            warning('Error in line click: %s', ME.message);
        end
    end

    function setGroupVar(varName, activeCats)
        data = hContainer.UserData;
        items = get(data.ddGrpBy, 'String');
        idx = find(strcmp(items, varName));

        needsUpdate = false;

        % 1. Change Variable
        if ~isempty(idx) && idx ~= get(data.ddGrpBy, 'Value')
            set(data.ddGrpBy, 'Value', idx);
            % Update checkboxes Manually
            populateCheckboxes(hContainer, data, data.ddGrpBy, data.pnlGrpBy, 'chkGrpBy');
            data = hContainer.UserData;
            needsUpdate = true;
        end

        % 2. Apply Filters
        if exist('activeCats', 'var') && isfield(data, 'chkGrpBy') && ~isempty(data.chkGrpBy)
            for i = 1:length(data.chkGrpBy)
                catStr = get(data.chkGrpBy(i), 'String');
                val = ismember(catStr, activeCats);
                if get(data.chkGrpBy(i), 'Value') ~= val
                    set(data.chkGrpBy(i), 'Value', val);
                    needsUpdate = true;
                end
            end
        end

        if needsUpdate
            onUpdatePlot(hContainer, []);
        end
    end

end

function tf = permember(str, allowed)
tf = any(strcmpi(str, allowed));
end
