function hFig = tblGUI_xy(xVec, dataTbl, varargin)
% TBLGUI_XY Interactive visualization of Table variables against a vector.
%
%   hFig = tblGUI_xy(xVec, dataTbl, varargin) plots column 'yVar' from 'dataTbl'
%   against 'xVec'.
%   - "Plot By": Splits data into separate tiles (subplots).
%   - "Group By": Groups data within each tile by color.
%
%   INPUTS:
%       xVec    (numeric)  Vector for X-axis.
%       dataTbl (table)    Data table containing variables to plot.
%
%   OPTIONAL KEY-VALUE PAIRS:
%       'yVar'             (char/str) Name of the column in dataTbl to plot.
%                          If empty, auto-selects a numeric column matching
%                          xVec dimensions.
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
addRequired(p, 'xVec', @isnumeric);
addRequired(p, 'dataTbl', @istable);
addOptional(p, 'yVar', [], @(x) ischar(x) || isstring(x) || isempty(x));
addParameter(p, 'tileFlow', 'vertical', @(x) permember(x, {'flow', 'vertical', 'horizontal'}));
parse(p, xVec, dataTbl, varargin{:});

xVec = p.Results.xVec;
dataTbl = p.Results.dataTbl;
yVar = p.Results.yVar;
tileFlow = p.Results.tileFlow;

% Auto-select yVar if empty
if isempty(yVar)
    varNames = dataTbl.Properties.VariableNames;
    xLen = length(xVec);

    for iVar = 1:length(varNames)
        raw = dataTbl.(varNames{iVar});
        if isnumeric(raw) && size(raw, 2) == xLen
            yVar = varNames{iVar};
            break;
        elseif iscell(raw) && length(raw) > 0 && length(raw{1}) == xLen
            yVar = varNames{iVar};
            break;
        end
    end

    if isempty(yVar)
        error('No suitable variable found in table matching xVec length.');
    end
    fprintf('Auto-selected yVar: %s\n', yVar);
end

if ~ismember(yVar, dataTbl.Properties.VariableNames)
    error('Variable "%s" not found in table.', yVar);
end

%% ========================================================================
%  INITIALIZATION
%  ========================================================================

% Identify potential grouping variables
allVars = dataTbl.Properties.VariableNames;
catVars = allVars(varfun(@(x) iscategorical(x) || isstring(x) || islogical(x),...
    dataTbl, 'OutputFormat', 'uniform'));
catVars = [{'None'}, catVars];

% Figure Setup
% Use explicit pixels to act nicely on screens
hFig = figure('Name', sprintf('XY Plot: %s', yVar), 'NumberTitle', 'off', ...
    'Units', 'pixels', 'Position', [100, 100, 1400, 800]);

% Store GUI Data
guiData = struct();
guiData.xVec = xVec;
guiData.dataTbl = dataTbl;
guiData.yVar = yVar;
guiData.catVars = catVars;
guiData.tileFlow = tileFlow;
guiData.chkPlotBy = []; % Handles for checkboxes
guiData.chkGrpBy = [];  % Handles for checkboxes
guiData.colors = lines(20); % Pre-define colors

%% ========================================================================
%  LAYOUT
%  ========================================================================

% Side Panel for Controls (Left)
panelW = 0.1;
% Create Left Panel
hPanel = uipanel('Parent', hFig, 'Units', 'normalized', ...
    'Position', [0, 0, panelW, 1]);

% Create Right Panel for Plots
hPanelRight = uipanel('Parent', hFig, 'Units', 'normalized', ...
    'Position', [panelW, 0, 1-panelW, 1], 'BorderType', 'none');

% Main Plotting Area (Right)
% Parent tiledlayout to the RIGHT PANEL to avoid bleeding
guiData.hLayout = tiledlayout(hPanelRight, 'flow', ...
    'TileSpacing', 'tight', 'Padding', 'compact');
guiData.hPanelRight = hPanelRight; % Store parent panel for layout recreation

% --- Controls in Side Panel ---
ctlH = 0.03;
ctlGap = 0.01;
currY = 0.95;

% PLOT BY (Tiles)
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

% GROUP BY (Colors)
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

% Update Button
guiData.btnUpdate = uicontrol('Parent', hPanel, 'Style', 'pushbutton', ...
    'String', 'Update Plot', ...
    'Units', 'normalized', 'Position', [0.05, 0.02, 0.9, 0.05], ...
    'Callback', @onUpdatePlot, 'FontWeight', 'bold', 'FontSize', 10);


guidata(hFig, guiData);

% Initial State
onPlotByChange(hFig, []);
onGrpByChange(hFig, []);

% Initial Plot
onUpdatePlot(hFig, []);

%% ========================================================================
%  CALLBACKS
%  ========================================================================

    function onPlotByChange(src, ~)
        data = guidata(src);
        populateCheckboxes(src, data, data.ddPlotBy, data.pnlPlotBy, 'chkPlotBy');
    end

    function onGrpByChange(src, ~)
        data = guidata(src);
        populateCheckboxes(src, data, data.ddGrpBy, data.pnlGrpBy, 'chkGrpBy');
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
            for iVar = 1:nCats
                yPos = 1 - iVar*h;
                chkHandles(iVar) = uicontrol('Parent', hPanel, 'Style', 'checkbox', ...
                    'String', cats{iVar}, 'Units', 'normalized', ...
                    'Position', [0, yPos, w, h], 'Value', 1); % Default Select All
            end
        end

        data.(storeField) = chkHandles;
        guidata(src, data);
    end

    function onUpdatePlot(src, ~)
        data = guidata(src);

        % Manage Layout based on Arrangement
        % If layout type changes, we might need to be careful, but just
        % configuring it here is fine.
        nTiles = 1; % Temporary default

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

        % Recreate layout to support different arrangements
        delete(data.hLayout);
        data.hLayout = tiledlayout(data.hPanelRight, data.tileFlow, ...
            'TileSpacing', 'tight', 'Padding', 'compact');
        guidata(src, data); 

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

        % Retrieve Data
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
                plot(hAx, data.xVec, subY', 'Color', [clr, 0.05],...
                    'LineWidth', 0.5, 'HandleVisibility', 'off');

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
            axis(hAx, 'tight'); % Set X tightly

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

    end

end

function tf = permember(str, allowed)
tf = any(strcmpi(str, allowed));
end
