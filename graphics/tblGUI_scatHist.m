function hFig = tblGUI_scatHist(tbl, varargin)
% TBLGUI_SCATHIST Interactive scatter plot with marginal histograms and grouping.
%
%   tblGUI_scatHist(tbl, ...) opens a GUI to visualize the table 'tbl'.
%   Allows for dynamic variable selection, grouping, and interactive point
%   selection.
%
%   INPUTS:
%       tbl         - (table) The data table to visualize.
%
%   OPTIONAL KEY-VALUE PAIRS:
%       'xVar'      - (char) Initial X variable name.
%       'yVar'      - (char) Initial Y variable name.
%       'szVar'     - (char) Initial Size variable name.
%       'grpVar'    - (char) Initial Group variable name.
%       'clr'       - (m x 3) RGB color matrix for groups.
%       'alpha'     - (scalar) Marker transparency (0-1).
%       'varsExclude' - (cell) Variables to exclude from dropdowns.
%                     Default: {'UnitID', 'Name', 'Mouse', 'File'}.
%       'Parent'    - (handle) Parent container (figure/panel).
%       'SelectionCallback' - (func) Callback for point selection.
%       'GroupByCallback'   - (func) Callback for group selection.
%
%   EXAMPLE:
%       tblGUI_scatHist(fetTbl, 'xVar', 'tp', 'yVar', 'mfr', 'grpVar', 'Name');
%
%   See also: LME_ANALYSE

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addParameter(p, 'varsExclude', {'UnitID', 'Mouse', 'File'}, @iscell);
addParameter(p, 'xVar', '', @(x) ischar(x) || isstring(x) || isempty(x));
addParameter(p, 'yVar', '', @(x) ischar(x) || isstring(x) || isempty(x));
addParameter(p, 'szVar', '', @(x) ischar(x) || isstring(x) || isempty(x));
addParameter(p, 'grpVar', '', @(x) ischar(x) || isstring(x) || isempty(x));
addParameter(p, 'alpha', 0.6, @isnumeric);
addParameter(p, 'clr', [], @(x) isempty(x) || size(x,2)==3);
addParameter(p, 'Parent', [], @(x) isempty(x) || isgraphics(x));
addParameter(p, 'SelectionCallback', [], @(x) isempty(x) || isa(x, 'function_handle'));
addParameter(p, 'GroupByCallback', [], @(x) isempty(x) || isa(x, 'function_handle'));

parse(p, tbl, varargin{:});

% Unpack Results
tbl         = p.Results.tbl;
varsExclude = p.Results.varsExclude;
hParent     = p.Results.Parent;
selCbk      = p.Results.SelectionCallback;
grpCbk      = p.Results.GroupByCallback;
defAlpha    = p.Results.alpha;
defClr      = p.Results.clr;

% Initial Variable Inputs
xVarIn      = p.Results.xVar;
yVarIn      = p.Results.yVar;
szVarIn     = p.Results.szVar;
grpVarIn    = p.Results.grpVar;

%% ========================================================================
%  INITIALIZATION
%  ========================================================================

% --- Variable Classification ---
allVars = tbl.Properties.VariableNames;

% Numeric Variables (for X, Y, Size)
numericVars = allVars(varfun(@isnumeric, tbl, 'OutputFormat', 'uniform'));
numericVars = setdiff(numericVars, varsExclude);

if isempty(numericVars)
    error('Input table must contain at least one numeric variable.');
end

% Categorical Variables (for Grouping)
% We treat strings and logicals as potential categories too.
catVars = allVars(varfun(@(x) iscategorical(x) || isstring(x) || islogical(x),...
    tbl, 'OutputFormat', 'uniform'));
catVars = setdiff(catVars, varsExclude);

% --- Defaults ---
% Prioritize inputs, otherwise using simple heuristics (1st and 2nd vars)
curX = numericVars{1};
curY = numericVars{min(2, length(numericVars))};
curSize = 'None';
curGrp  = 'None';

if ~isempty(catVars), curGrp = catVars{1}; end

% Apply Overrides
if ~isempty(xVarIn) && ismember(xVarIn, numericVars),   curX = xVarIn;   end
if ~isempty(yVarIn) && ismember(yVarIn, numericVars),   curY = yVarIn;   end
if ~isempty(szVarIn) && ismember(szVarIn, numericVars), curSize = szVarIn; end
if ~isempty(grpVarIn) && ismember(grpVarIn, catVars),   curGrp = grpVarIn; end

% --- Figure Setup ---
if isempty(hParent)
    hContainer = figure('Name', 'Table Visualizer', 'NumberTitle', 'off', ...
        'Position', [100, 100, 1110, 900], 'MenuBar', 'none', 'ToolBar', 'figure');
    hFig = hContainer;
else
    hContainer = hParent;
    hFig = ancestor(hContainer, 'figure');
end

% --- Data Storage ---
% Store state in the container's UserData
guiData = struct();
guiData.tbl             = tbl;
guiData.numericVars     = numericVars;
guiData.catVars         = catVars;
guiData.pointHandles    = [];   % To store scatter handles
guiData.polyRoi         = [];   % To store polygon ROI
guiData.defAlpha        = defAlpha;
guiData.defClr          = defClr;
guiData.selCbk          = selCbk;
guiData.grpCbk          = grpCbk;
guiData.highlightFcn    = @highlightPoints; % Expose function
guiData.setGroupVarFcn  = @setGroupVar;     % Expose function
guiData.setXYVarsFcn    = @setXYVars;       % Expose function
guiData.hHighlight      = [];   % Store handle for external highlight
guiData.fitHandles      = [];   % Store handles for regression lines
guiData.chkGrp          = [];   % Checkbox handles
guiData.hEquality       = [];   % Store handle for equality line

%% ========================================================================
%  LAYOUT
%  ========================================================================

% Margins and spacing
panelW      = 0.2;  % Width of control panel
marg        = 0.05;
bottomMarg  = 0.1;

% Main Scatter: Bottom-Left (relative to plot area)
scatterW    = 0.55;
scatterH    = 0.68;
startX      = panelW + marg;
startY      = bottomMarg;

axScatterPos = [startX, startY, scatterW, scatterH];
% Top Histogram (X dist)
axHistXPos   = [startX, startY + scatterH + 0.02, scatterW, 0.15];
% Right Histogram (Y dist)
axHistYPos   = [startX + scatterW + 0.02, startY, 0.15, scatterH];

guiData.hAxScatter = axes('Parent', hContainer, 'Position', axScatterPos);
guiData.hAxHistX   = axes('Parent', hContainer, 'Position', axHistXPos);
guiData.hAxHistY   = axes('Parent', hContainer, 'Position', axHistYPos);

% Link axes for zooming
linkaxes([guiData.hAxScatter, guiData.hAxHistX], 'x');

%% ========================================================================
%  CONTROLS
%  ========================================================================

% Control Panel Positioning
ctlTop      = 0.96;
ctlH        = 0.04;
ctlGap      = 0.01;
ctlW        = panelW - 0.02;     % Full width
ctlW_Half   = ctlW / 2 - 0.005;  % Half width (2-col)
ctlW_Third  = (ctlW - 2*ctlGap) / 3; % Third width (3-col)

ctlX        = 0.01;
ctlX_Right  = ctlX + ctlW_Half + 0.01;

% 3-col positions
ctlX_1      = ctlX;
ctlX_2      = ctlX + ctlW_Third + ctlGap;
ctlX_3      = ctlX + 2*(ctlW_Third + ctlGap);

% --- X Variable ---
uicontrol('Parent', hContainer, 'Style', 'text', 'String', 'X Variable / Scale:', ...
    'Units', 'normalized', 'Position', [ctlX, ctlTop - ctlH, ctlW, ctlH*0.7], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');
guiData.ddX = uicontrol('Parent', hContainer, 'Style', 'popupmenu', ...
    'String', numericVars, 'Units', 'normalized', ...
    'Position', [ctlX, ctlTop - 2*ctlH, ctlW_Half, ctlH], ...
    'Callback', @onUpdatePlot);
guiData.ddXScale = uicontrol('Parent', hContainer, 'Style', 'popupmenu', ...
    'String', {'Linear', 'Log'}, 'Units', 'normalized', ...
    'Position', [ctlX_Right, ctlTop - 2*ctlH, ctlW_Half, ctlH], ...
    'Callback', @onUpdatePlot);

% --- Discretize / Adaptive / Prct (Row 3) ---
yPosChk = ctlTop - 3*ctlH - ctlGap;

guiData.chkDisc = uicontrol('Parent', hContainer, 'Style', 'checkbox', ...
    'String', 'Disc', 'Units', 'normalized', ...
    'Position', [ctlX_1, yPosChk, ctlW_Third, ctlH], ...
    'Callback', @onUpdatePlot);

guiData.chkAdapt = uicontrol('Parent', hContainer, 'Style', 'checkbox', ...
    'String', 'Adpt', 'Units', 'normalized', ...
    'Position', [ctlX_2, yPosChk, ctlW_Third, ctlH], ...
    'Callback', @onUpdatePlot);

guiData.chkPrct = uicontrol('Parent', hContainer, 'Style', 'checkbox', ...
    'String', 'Prct', 'Units', 'normalized', ...
    'Position', [ctlX_3, yPosChk, ctlW_Third, ctlH], ...
    'Callback', @onUpdatePlot);

% --- Y Variable ---
currYTop = ctlTop - 3*ctlH - 2*ctlGap - ctlH;
uicontrol('Parent', hContainer, 'Style', 'text', 'String', 'Y Variable / Scale:', ...
    'Units', 'normalized', 'Position', [ctlX, currYTop, ctlW, ctlH*0.7], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');
guiData.ddY = uicontrol('Parent', hContainer, 'Style', 'popupmenu', ...
    'String', numericVars, 'Units', 'normalized', ...
    'Position', [ctlX, currYTop - ctlH, ctlW_Half, ctlH], ...
    'Callback', @onUpdatePlot);
guiData.ddYScale = uicontrol('Parent', hContainer, 'Style', 'popupmenu', ...
    'String', {'Linear', 'Log', 'X-Linked'}, 'Units', 'normalized', ...
    'Position', [ctlX_Right, currYTop - ctlH, ctlW_Half, ctlH], ...
    'Callback', @onUpdatePlot);

% --- Fit Option (Moved up, Inline) ---
yFit = currYTop - 2*ctlH - ctlGap;
uicontrol('Parent', hContainer, 'Style', 'text', 'String', 'Fit:', ...
    'Units', 'normalized', 'Position', [ctlX, yFit + ctlH*0.1, ctlW*0.3, ctlH*0.6], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');
guiData.ddFit = uicontrol('Parent', hContainer, 'Style', 'popupmenu', ...
    'String', {'None', 'Linear', 'Ortho'}, 'Units', 'normalized', ...
    'Position', [ctlX + ctlW*0.35, yFit, ctlW*0.65, ctlH], ...
    'Callback', @onUpdatePlot);

% --- Size Variable (Inline) ---
ySize = yFit - ctlH - ctlGap;
uicontrol('Parent', hContainer, 'Style', 'text', 'String', 'Size:', ...
    'Units', 'normalized', 'Position', [ctlX, ySize + ctlH*0.1, ctlW*0.3, ctlH*0.6], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');
guiData.ddSize = uicontrol('Parent', hContainer, 'Style', 'popupmenu', ...
    'String', [{'None'}, numericVars], 'Units', 'normalized', ...
    'Position', [ctlX + ctlW*0.35, ySize, ctlW*0.65, ctlH], ...
    'Callback', @onUpdatePlot);

% --- Group Variable (Inline) ---
yGrp = ySize - ctlH - ctlGap;
uicontrol('Parent', hContainer, 'Style', 'text', 'String', 'Group:', ...
    'Units', 'normalized', 'Position', [ctlX, yGrp + ctlH*0.1, ctlW*0.3, ctlH*0.6], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');
guiData.ddGrp = uicontrol('Parent', hContainer, 'Style', 'popupmenu', ...
    'String', [{'None'}, catVars], 'Units', 'normalized', ...
    'Position', [ctlX + ctlW*0.35, yGrp, ctlW*0.65, ctlH], ...
    'Callback', @onGrpChange);

% Set Initial Values
set(guiData.ddX, 'Value', find(strcmp(numericVars, curX)));
set(guiData.ddY, 'Value', find(strcmp(numericVars, curY)));

% Size ('None' is 1)
idxSz = find(strcmp(numericVars, curSize));
if isempty(idxSz), set(guiData.ddSize, 'Value', 1);
else,              set(guiData.ddSize, 'Value', idxSz + 1); end

% Group ('None' is 1)
idxG = find(strcmp(catVars, curGrp));
if isempty(idxG), set(guiData.ddGrp, 'Value', 1);
else,             set(guiData.ddGrp, 'Value', idxG + 1); end

% Fit ('None' is 1)
set(guiData.ddFit, 'Value', 1);

% --- Filter Panel ---
grpDDBottom = yGrp;
panelTop    = grpDDBottom - 0.02;
panelBottom = 0.20; % Adjusted for buttons
panelH      = panelTop - panelBottom;

guiData.pnlGrp = uipanel('Parent', hContainer, 'BorderType', 'none', ...
    'Units', 'normalized', 'Position', [ctlX, panelBottom, ctlW, panelH]);

% --- Buttons ---
% Side-by-side
btnY = 0.13;
guiData.btnSelect = uicontrol('Parent', hContainer, 'Style', 'pushbutton', ...
    'String', 'Slct Grp', ...
    'Units', 'normalized', 'Position', [ctlX, btnY, ctlW_Half, 0.05], ...
    'Callback', @onSelectRegion, 'FontWeight', 'bold');

guiData.btnSelectDot = uicontrol('Parent', hContainer, 'Style', 'pushbutton', ...
    'String', 'Slct Dot', ...
    'Units', 'normalized', 'Position', [ctlX_Right, btnY, ctlW_Half, 0.05], ...
    'Callback', @onSelectDot, 'FontWeight', 'bold');

guiData.btnSave = uicontrol('Parent', hContainer, 'Style', 'pushbutton', ...
    'String', 'Save Table to Workspace', ...
    'Units', 'normalized', 'Position', [ctlX, 0.05, ctlW, 0.06], ...
    'Callback', @onSaveTable);

% Store guidata
hContainer.UserData = guiData;

% Finalize
hFig = hContainer;

% Trigger Initial Update
onGrpChange(hContainer);
onUpdatePlot(hContainer, []);


%% ========================================================================
%  CALLBACKS
%  ========================================================================

    function onGrpChange(src, ~)
        % Called when the Group Dropdown changes. Re-populates the checkbox list.
        data = hContainer.UserData;
        populateCheckboxes(data);

        % Trigger External Callback
        idxGrp   = get(data.ddGrp, 'Value');
        grpItems = get(data.ddGrp, 'String');
        grpName  = grpItems{idxGrp};

        if ~isempty(data.grpCbk)
            % New group var -> Enable ALL categories by default
            data = hContainer.UserData; % Reload
            if isempty(data.chkGrp)
                activeCats = {};
            else
                allCats = arrayfun(@(x) string(get(x, 'String')), data.chkGrp);
                activeCats = cellstr(allCats);
            end
            data.grpCbk(grpName, activeCats, hContainer);
        end

        onUpdatePlot(src, []);
    end

    function onFilterChange(~, ~)
        % Called when a single checkbox is toggled.
        data = hContainer.UserData;
        idxGrp   = get(data.ddGrp, 'Value');
        grpItems = get(data.ddGrp, 'String');
        grpName  = grpItems{idxGrp};

        fprintf('Filter Change: %s\n', grpName);
        onUpdatePlot(hContainer, []);

        % Trigger External Callback
        if ~isempty(data.grpCbk)
            if isempty(data.chkGrp)
                activeCats = {};
            else
                selectedIdx = arrayfun(@(x) get(x, 'Value'), data.chkGrp);
                allCats     = arrayfun(@(x) string(get(x, 'String')), data.chkGrp);
                activeCats  = cellstr(allCats(logical(selectedIdx)));
            end
            data.grpCbk(grpName, activeCats, hContainer);
        end
    end

    function populateCheckboxes(data)
        % Helper to fill the filter panel with checkboxes for each category.
        delete(data.pnlGrp.Children);
        idxGrp   = get(data.ddGrp, 'Value');
        grpItems = get(data.ddGrp, 'String');
        grpName  = grpItems{idxGrp};

        if strcmp(grpName, 'None')
            data.chkGrp = [];
        else
            raw = data.tbl.(grpName);
            if islogical(raw), raw = categorical(raw); end
            if ~iscategorical(raw), raw = categorical(raw); end

            % Get categories present in the data
            cats = categories(raw);
            cats = cats(ismember(cats, unique(raw)));

            nCats = length(cats);
            h = 1 / max(10, nCats + 1);
            w = 0.9;
            data.chkGrp = gobjects(1, nCats); % Initialize to clear old handles

            for iCat = 1:nCats
                yPos = 1 - iCat*h;
                data.chkGrp(iCat) = uicontrol('Parent', data.pnlGrp, 'Style', 'checkbox', ...
                    'String', cats{iCat}, 'Units', 'normalized', ...
                    'Position', [0, yPos, w, h], 'Value', 1, ...
                    'Callback', @onFilterChange);
            end
        end
        hContainer.UserData = data;
    end

    function onUpdatePlot(~, ~)
        % Main plotting function.
        data = hContainer.UserData;
        tbl  = data.tbl;

        % --- Retrieve Selections ---
        idxX    = get(data.ddX, 'Value');
        idxY    = get(data.ddY, 'Value');
        idxSize = get(data.ddSize, 'Value');
        idxGrp  = get(data.ddGrp, 'Value');
        idxFit  = get(data.ddFit, 'Value');

        xName = data.numericVars{idxX};
        yName = data.numericVars{idxY};

        % --- Size Logic (Quadratic Scaling) ---
        if idxSize > 1
            sizeVarItems = get(data.ddSize, 'String');
            sizeName     = sizeVarItems{idxSize};
            rawSz        = tbl.(sizeName);
            rawSz        = rawSz + abs(min(rawSz)); % Shift to positive

            % Robust normalization
            lims  = prctile(rawSz, [5 95]);
            minSz = lims(1);
            maxSz = lims(2);

            if maxSz > minSz
                clippedSz = max(min(rawSz, maxSz), minSz);
                normSz    = (clippedSz - minSz) / (maxSz - minSz);

                % Quadratic transform (x^5 for strong contrast)
                normSq = normSz .^ 5;

                % Map to range 20-100
                szData = 20 + 80 * normSq;
            else
                szData = repmat(20, length(rawSz), 1);
            end
        else
            szData = repmat(20, height(tbl), 1);
        end

        % --- Group Logic & Colors ---
        grpName     = '';
        groups      = ones(height(tbl), 1);
        grpLabels   = {'All'};

        % Initialize colors container
        grpColorsMap = [];
        fullCatList  = {};

        grpItems = get(data.ddGrp, 'String');
        if idxGrp > 1 % "None" is 1
            grpName = grpItems{idxGrp};
            rawGrp  = tbl.(grpName);

            if islogical(rawGrp)
                rawGrp = categorical(rawGrp);
            elseif ~iscategorical(rawGrp)
                rawGrp = categorical(rawGrp);
            end

            groups = rawGrp;
            grpLabels = categories(groups);

            % Generate colors based on ALL potential categories to ensure stability
            % (Logic: If we have 3 groups A,B,C, A should always be Blue, even if B is unchecked)
            fullCatList = categories(rawGrp);

            % Filter to only those present in data if desired, but for stability
            % 'categories' usually returns the full definition.
            % Let's intersect with present data to be safe against unused categories,
            % but still take the FULL set of PRESENT data, not the FILTERED set.
            fullCatList = fullCatList(ismember(fullCatList, unique(rawGrp)));
            nTotalCats  = length(fullCatList);

            if ~isempty(data.defClr) && size(data.defClr,1) >= nTotalCats
                baseColors = data.defClr;
            else
                baseColors = lines(nTotalCats);
            end

            % Adjust 'grpLabels' based on user checkboxes (Filtering)
            if ~isempty(data.chkGrp)
                validH = isgraphics(data.chkGrp);
                if any(validH)
                    selectedIdx = arrayfun(@(x) get(x, 'Value'), data.chkGrp(validH));
                    allCatsStr  = arrayfun(@(x) string(get(x, 'String')), data.chkGrp(validH));
                    activeCats  = cellstr(allCatsStr(logical(selectedIdx)));

                    % Filter the labels to iterate over
                    grpLabels = intersect(grpLabels, activeCats, 'stable');
                end
            end

            isGrpActive = true;
        else
            isGrpActive = false;
            baseColors  = lines(1);
        end

        % Get Data Vectors
        xData = tbl.(xName);
        yData = tbl.(yName);

        % Discretization Flag
        isDisc = get(data.chkDisc, 'Value');

        % --- Scales & Limits ---
        scaleX = 'linear';
        if get(data.ddXScale, 'Value') == 2, scaleX = 'log'; end

        yScaleVal = get(data.ddYScale, 'Value');
        isLinked  = (yScaleVal == 3);

        if isLinked
            scaleY = scaleX;
        else
            scaleY = 'linear';
            if yScaleVal == 2, scaleY = 'log'; end
        end

        % Calculate Limits (for Bins)
        xLim = calcLimits(xData, scaleX);
        if xLim(1) >= xLim(2)
            if strcmp(scaleX, 'log'), xLim = [xLim(1)/1.1, xLim(1)*1.1];
            else,                     xLim = [xLim(1)-0.5, xLim(1)+0.5]; end
        end

        yLim = calcLimits(yData, scaleY);
        if yLim(1) >= yLim(2)
            if strcmp(scaleY, 'log'), yLim = [yLim(1)/1.1, yLim(1)*1.1];
            else,                     yLim = [yLim(1)-0.5, yLim(1)+0.5]; end
        end

        if isLinked
            jointMin = min(xLim(1), yLim(1));
            jointMax = max(xLim(2), yLim(2));
            xLim = [jointMin, jointMax];
            yLim = [jointMin, jointMax];
        end

        % Calculate Global Bin Edges
        nBins = 30;
        if strcmp(scaleX, 'log')
            % Ensure positive for log
            if xLim(1) <= 0, xLim(1) = min(xData(xData>0)); if isempty(xLim(1)), xLim(1)=0.1; end; end
            xEdges = logspace(log10(xLim(1)), log10(xLim(2)), nBins);
        else
            xEdges = linspace(xLim(1), xLim(2), nBins);
        end

        if strcmp(scaleY, 'log')
            if yLim(1) <= 0, yLim(1) = min(yData(yData>0)); if isempty(yLim(1)), yLim(1)=0.1; end; end
            yEdges = logspace(log10(yLim(1)), log10(yLim(2)), nBins);
        else
            yEdges = linspace(yLim(1), yLim(2), nBins);
        end

        % Set Axis Scales Immediately
        set(data.hAxScatter, 'XScale', scaleX, 'YScale', scaleY);
        set(data.hAxHistX,   'XScale', scaleX, 'YScale', 'linear');
        set(data.hAxHistY,   'XScale', 'linear', 'YScale', scaleY);

        % --- Drawing ---
        % Clean up previous fit lines (handled explicitly due to HandleVisibility=off)
        if ~isempty(data.fitHandles)
            delete(data.fitHandles(isgraphics(data.fitHandles)));
        end
        data.fitHandles = [];

        if ~isempty(data.hEquality)
            delete(data.hEquality(isgraphics(data.hEquality)));
        end
        data.hEquality = [];

        delete(allchild(data.hAxScatter));
        delete(allchild(data.hAxHistX));
        delete(allchild(data.hAxHistY));

        hold(data.hAxScatter, 'on');
        hold(data.hAxHistX, 'on');
        hold(data.hAxHistY, 'on');

        % Plotting Loop
        % Store stats for second pass
        statsData = {};

        % Plotting Loop (Pass 1: Scatter & Histograms)
        for iG = 1:length(grpLabels)
            if isGrpActive
                currCat = grpLabels{iG};
                idx     = (groups == currCat);

                % Color Stability Logic:
                % Find the index of this category in the FULL list of categories
                [~, globalIdx] = ismember(currCat, fullCatList);
                if globalIdx > 0
                    cIdx = mod(globalIdx-1, size(baseColors,1)) + 1;
                    cG   = baseColors(cIdx, :);
                else
                    cG = [0 0 0]; % Fallback
                end
            else
                idx     = true(height(tbl), 1);
                currCat = 'All';
                cG      = baseColors(1, :);
            end

            if ~any(idx), continue; end

            % Extract Subset
            xG  = xData(idx);
            yG  = yData(idx);
            szG = szData(idx);

            nPoints = sum(idx);

            % Scatter Plot Params
            scatAlpha = data.defAlpha;
            if isDisc
                scatAlpha = scatAlpha * 0.1;
                szG = szG * 0.3;
            end

            % Get Fit Type
            fitTypes = {'None', 'Linear', 'Ortho'};
            curFit   = fitTypes{idxFit};

            % Prepare Label
            lbl = sprintf('%s (n=%d)', string(currCat), nPoints);

            % Plot using plot_scat
            [hS, hF] = plot_scat([], xG, yG, ...
                'hAx', data.hAxScatter, ...
                'sz', szG, ...
                'c', cG, ...
                'alpha', scatAlpha, ...
                'marker', 'o', ...
                'fitType', curFit, ...
                'flgStats', true, ... % Calculate and append stats
                'dispName', lbl);
            
            % Post-process Scatter (Interaction)
            set(hS, 'HitTest', 'off', 'PickableParts', 'none');
            
            % Store fit handles
            if isgraphics(hF)
               data.fitHandles = [data.fitHandles; hF];
            end


            % Binned Stats Calculation
            if isDisc
                isAdapt = get(data.chkAdapt, 'Value');
                isPrct  = get(data.chkPrct, 'Value');

                sStruct = calcBinnedStats(xG, yG, cG, isAdapt, isPrct, xEdges, scaleX);

                if ~isempty(sStruct)
                    statsData{end+1} = sStruct; %#ok<AGROW>
                end
            end

            % Histogram X
            if strcmp(scaleX, 'log'), normX = 'probability'; else, normX = 'pdf'; end
            plot_hist([], xG, 'hAx', data.hAxHistX, 'bins', xEdges, 'c', cG, ...
                'scale', scaleX, 'orient', 'vertical', 'norm', normX, ...
                'flgKDE', true, 'flgStat', false);

            % Histogram Y (Horizontal)
            if strcmp(scaleY, 'log'), normY = 'probability'; else, normY = 'pdf'; end
            plot_hist([], yG, 'hAx', data.hAxHistY, 'bins', yEdges, 'c', cG, ...
                'scale', scaleY, 'orient', 'horizontal', 'norm', normY, ...
                'flgKDE', true, 'flgStat', false);
        end

        % Plotting Loop (Pass 2: Error Bars on TOP)
        if isDisc
            for k = 1:length(statsData)
                s = statsData{k};

                % Filter NaNs
                valid = ~isnan(s.meds);
                if ~any(valid), continue; end

                % Err bar color (darker)
                errColor = max(0, s.color * 0.7);

                hErr = errorbar(data.hAxScatter, s.centers(valid), s.meds(valid), ...
                    s.neg(valid), s.pos(valid), ...
                    'Color', errColor, 'LineWidth', 2, 'LineStyle', '-', ...
                    'Marker', 'o', 'MarkerSize', 6, 'MarkerFaceColor', errColor, ...
                    'CapSize', 0);

                % Hide from legend
                hErr.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end
        end

        % Labels & Aesthetics
        xlabel(data.hAxScatter, xName, 'Interpreter', 'none');
        ylabel(data.hAxScatter, yName, 'Interpreter', 'none');

        if isGrpActive
            legend(data.hAxScatter, 'Location', 'best', 'Interpreter', 'tex');
        end

        % Equality Line
        if isLinked
            % Re-plot equality line since we cleared axes
            data.hEquality = plot(data.hAxScatter, xLim, xLim, 'k--', 'LineWidth', 1, ...
                'HitTest', 'off', 'PickableParts', 'none', 'HandleVisibility', 'off');
        end

        % Apply Padding and Set Limits
        applyPaddedLimits(data.hAxScatter, xLim, yLim, scaleX, scaleY);

        grid(data.hAxScatter, 'on');
        data.hAxHistX.XAxis.Visible = 'off';
        data.hAxHistX.YAxis.Visible = 'off';
        data.hAxHistY.XAxis.Visible = 'off';
        data.hAxHistY.YAxis.Visible = 'off';

        linkaxes([data.hAxScatter, data.hAxHistX], 'x');

        hold(data.hAxScatter, 'off');
        hold(data.hAxHistX, 'off');
        hold(data.hAxHistY, 'off');

        % Update UI State
        if isGrpActive
            set(data.btnSelect, 'Enable', 'on', 'String', 'Select Group');
            set(data.btnSelectDot, 'Enable', 'on', 'String', 'Select Dot');
        else
            set(data.btnSelect, 'Enable', 'off', 'String', 'Select Group');
            set(data.btnSelectDot, 'Enable', 'off', 'String', 'Select Dot');
        end
        
        % Save State (including fitHandles)
        hContainer.UserData = data;
    end

% --- Helper: Limit Calculation ---
    function lims = calcLimits(vec, scaleType)
        if strcmp(scaleType, 'log')
            v = vec(~isnan(vec) & ~isinf(vec) & vec > 0);
            if isempty(v), lims = [0.1 1]; else, lims = [min(v), max(v)]; end
        else
            v = vec(~isnan(vec) & ~isinf(vec));
            if isempty(v), lims = [0 1]; else, lims = [min(v), max(v)]; end
        end
    end % EOF calcLimits

% --- Helper: Apply Padding ---
    function applyPaddedLimits(ax, xL, yL, sX, sY)
        % X
        if strcmp(sX, 'log')
            logMin = log10(xL(1)); logMax = log10(xL(2));
            span = logMax - logMin; if span==0, span=1; end
            xlim(ax, [10^(logMin - 0.05*span), 10^(logMax + 0.05*span)]);
        else
            span = diff(xL); if span==0, span=1; end
            xlim(ax, [xL(1) - 0.05*span, xL(2) + 0.05*span]);
        end
        % Y
        if strcmp(sY, 'log')
            logMin = log10(yL(1)); logMax = log10(yL(2));
            span = logMax - logMin; if span==0, span=1; end
            ylim(ax, [10^(logMin - 0.05*span), 10^(logMax + 0.05*span)]);
        else
            span = diff(yL); if span==0, span=1; end
            ylim(ax, [yL(1) - 0.05*span, yL(2) + 0.05*span]);
        end
    end % EOF applyPaddedLimits

% --- Helper: Binned Stats Calculation ---
    function sStruct = calcBinnedStats(xG, yG, cG, isAdapt, isPrct, xEdges, scaleX)
        % CALCBINNEDSTATS Calculates binned statistics (median, error bars).
        %
        %   Can use fixed edges (xEdges) or adaptive binning based on
        %   percentiles.

        sStruct = [];
        if isempty(xG) || isempty(yG), return; end

        % Data filtering
        % We already filtered loop-wise, but good to be safe if reusing
        valid = ~isnan(xG) & ~isnan(yG);
        xG = xG(valid);
        yG = yG(valid);
        if isempty(xG), return; end

        % 1. Strategy Setup
        if isAdapt
            % --- Adaptive: Edges based on quantiles ---
            % 16 bins for adaptive (0:100 linspace)
            pcts  = linspace(0, 100, 16);
            edges = unique(prctile(xG, pcts));

            if length(edges) < 2, return; end

            % Discretize
            binIdx = discretize(xG, edges);

            % We will iterate over the bins that actually exist
            uBins      = unique(binIdx(~isnan(binIdx)))';
            nBinSlots  = length(uBins);

            % For adaptive, we must calculate centers from data
            useDataCtr = true;
            preCenters = [];

        else
            % --- Fixed: Uniform/Log Edges ---
            edges = xEdges;

            % Discretize
            binIdx = discretize(xG, edges);

            % Pre-calculate centers
            if strcmp(scaleX, 'log')
                logEdges   = log10(edges);
                logCenters = (logEdges(1:end-1) + logEdges(2:end)) / 2;
                preCenters = 10.^logCenters;
            else
                preCenters = (edges(1:end-1) + edges(2:end)) / 2;
            end

            % We iterate over ALL defined bins for fixed timeline
            uBins      = 1:length(preCenters);
            nBinSlots  = length(uBins);
            useDataCtr = false;
        end

        % 2. Calculate Stats Loop
        bMeds = nan(1, nBinSlots);
        bLow  = nan(1, nBinSlots);
        bHigh = nan(1, nBinSlots);
        bCtrs = nan(1, nBinSlots);

        hasData = false;

        for iB = 1:nBinSlots
            currBinIdx = uBins(iB);

            % Identify points in bin
            inBin = (binIdx == currBinIdx);
            nIn   = sum(inBin);

            if nIn >= 5
                hasData = true;

                % Y-Stats
                % Y-Stats
                ySub = yG(inBin);
                
                if isPrct
                    % Request: If Pressed (Checked) -> Mean +/- SEM
                    mu  = mean(ySub, 'omitnan');
                    sd  = std(ySub, 'omitnan');
                    nn  = sum(~isnan(ySub));
                    sem = sd / sqrt(nn);
                    
                    bLow(iB)  = mu - sem;
                    bMeds(iB) = mu;
                    bHigh(iB) = mu + sem;
                else
                    % Default (Unchecked) -> Percentiles
                    p = prctile(ySub, [10, 50, 90]);
    
                    bLow(iB)  = p(1);
                    bMeds(iB) = p(2);
                    bHigh(iB) = p(3);
                end

                % X-Center
                if useDataCtr
                    xSub      = xG(inBin);
                    bCtrs(iB) = median(xSub, 'omitnan');
                else
                    bCtrs(iB) = preCenters(currBinIdx);
                end
            end
        end

        % 3. Package Result
        if hasData
            sStruct.centers = bCtrs;
            sStruct.meds    = bMeds;
            sStruct.neg     = bMeds - bLow;
            sStruct.pos     = bHigh - bMeds;
            sStruct.color   = cG;
        end

    end % EOF calcBinnedStats

    function onSelectRegion(~, ~)
        % Polygon selection tool
        data = hContainer.UserData;
        ax   = data.hAxScatter;
        roi  = drawpolygon(ax);

        if isempty(roi.Position)
            delete(roi); return;
        end

        % Process
        processSelection(hContainer, data, roi, false);
    end

    function onSelectDot(~, ~)
        % Single dot selection tool
        data = hContainer.UserData;
        ax   = data.hAxScatter;

        % Clear previous
        delete(findobj(ax, 'Type', 'images.roi.Point'));

        hPoint = drawpoint(ax, 'Label', 'Target');
        if isempty(hPoint.Position), delete(hPoint); return; end

        % Snap & Listen
        snapToData(hPoint, data);

        addlistener(hPoint, 'ROIClicked', ...
            @(roiSrc, evt) onDotDoubleClick(roiSrc, evt, hContainer));
        addlistener(hPoint, 'MovingROI', ...
            @(roiSrc, evt) snapToData(roiSrc, data));

        title(ax, 'Double-click to assign group. Drag to browse traces.', 'Color', 'r');
    end

    function snapToData(hPoint, data)
        % Snaps the ROI point to the nearest data point.
        pos = hPoint.Position;

        idxX  = get(data.ddX, 'Value');
        idxY  = get(data.ddY, 'Value');
        xName = data.numericVars{idxX};
        yName = data.numericVars{idxY};

        xData = data.tbl.(xName);
        yData = data.tbl.(yName);

        % Normalized Distance Calculation (to handle different axis scales)
        ax    = data.hAxScatter;
        xlims = ax.XLim;
        ylims = ax.YLim;

        xNorm = (xData - xlims(1)) / diff(xlims);
        yNorm = (yData - ylims(1)) / diff(ylims);
        pNorm = [(pos(1)-xlims(1))/diff(xlims), (pos(2)-ylims(1))/diff(ylims)];

        distSq      = (xNorm - pNorm(1)).^2 + (yNorm - pNorm(2)).^2;
        [~, minIdx] = min(distSq);

        newPos = [xData(minIdx), yData(minIdx)];

        % Update Position
        if sum((hPoint.Position - newPos).^2) > 1e-10
            hPoint.Position = newPos;
        end

        hPoint.UserData = minIdx;

        % Trigger Selection Callback (for browsing)
        inPoints = false(height(data.tbl), 1);
        inPoints(minIdx) = true;

        if isfield(data, 'selCbk') && ~isempty(data.selCbk)
            data.selCbk(inPoints);
        end
    end

    function onDotDoubleClick(hPoint, evt, ~)
        if strcmp(evt.SelectionType, 'double')
            data = hContainer.UserData;
            idx  = hPoint.UserData;

            % Select this point and open assignment dialog
            processSelectionIdx(hContainer, data, idx);

            delete(hPoint);
            title(data.hAxScatter, '');
        end
    end

    function processSelection(src, data, roi, isSinglePoint)
        idxX  = get(data.ddX, 'Value');
        idxY  = get(data.ddY, 'Value');
        xName = data.numericVars{idxX};
        yName = data.numericVars{idxY};

        xData = data.tbl.(xName);
        yData = data.tbl.(yName);

        inPoints  = inpolygon(xData, yData, roi.Position(:,1), roi.Position(:,2));
        nSelected = sum(inPoints);

        if nSelected == 0
            msgbox('No points selected.', 'Selection');
            delete(roi);
            return;
        end

        if isSinglePoint && ~isempty(data.selCbk)
            data.selCbk(inPoints);
        end

        assignGroup(src, data, inPoints, nSelected);
        delete(roi);
    end

    function processSelectionIdx(src, data, idx)
        inPoints      = false(height(data.tbl), 1);
        inPoints(idx) = true;

        if ~isempty(data.selCbk)
            data.selCbk(inPoints);
        end

        assignGroup(src, data, inPoints, 1);
    end

    function assignGroup(src, data, inPoints, nSelected)
        % Opens a dialog to assign selected points to a group.
        idxGrp   = get(data.ddGrp, 'Value');
        grpItems = get(data.ddGrp, 'String');
        grpName  = grpItems{idxGrp};

        currentGrpCol = data.tbl.(grpName);

        if iscategorical(currentGrpCol)
            cats = categories(currentGrpCol);
        elseif islogical(currentGrpCol)
            cats = {'false', 'true'};
        else
            cats = unique(string(currentGrpCol));
        end

        [indx, tf] = listdlg('PromptString', sprintf('Assign %d points to:', nSelected), ...
            'SelectionMode', 'single', ...
            'ListString', cats);

        if tf
            selectedCat = cats{indx};

            % Update Table
            if iscategorical(currentGrpCol)
                data.tbl.(grpName)(inPoints) = selectedCat;
            elseif islogical(currentGrpCol)
                val = strcmpi(selectedCat, 'true');
                data.tbl.(grpName)(inPoints) = val;
            else
                data.tbl.(grpName)(inPoints) = selectedCat;
            end

            hContainer.UserData = data;
            onUpdatePlot(src, []);
            fprintf('Updated %d points.\n', nSelected);
        end
    end

    function onSaveTable(~, ~)
        data = hContainer.UserData;
        assignin('base', 'fetTbl_mod', data.tbl);
        msgbox('Table saved to workspace as "fetTbl_mod".', 'Saved');
    end

    function setGroupVar(varName, activeCats)
        % External setter for group variable and filters
        data = hContainer.UserData;
        grpItems = get(data.ddGrp, 'String');
        idx = find(strcmp(grpItems, varName));

        needsUpdate = false;

        % 1. Change Variable
        if ~isempty(idx) && idx ~= get(data.ddGrp, 'Value')
            set(data.ddGrp, 'Value', idx);
            populateCheckboxes(data);
            data = hContainer.UserData; % Reload
            needsUpdate = true;
        end

        % 2. Apply Filters
        if exist('activeCats', 'var') && ~isempty(data.chkGrp)
            for i = 1:length(data.chkGrp)
                catStr = get(data.chkGrp(i), 'String');
                val    = ismember(catStr, activeCats);
                if get(data.chkGrp(i), 'Value') ~= val
                    set(data.chkGrp(i), 'Value', val);
                    needsUpdate = true;
                end
            end
        end

        if needsUpdate
            onUpdatePlot(hContainer, []);
        end
    end

    function setXYVars(xName, yName)
        % External setter for X and Y variables
        data = hContainer.UserData;

        idxX = find(strcmp(data.numericVars, xName));
        idxY = find(strcmp(data.numericVars, yName));

        if isempty(idxX) || isempty(idxY)
            warning('Variables %s or %s not found in numericVars.', xName, yName);
            return;
        end

        currX = get(data.ddX, 'Value');
        currY = get(data.ddY, 'Value');
        needsUpdate = false;

        if idxX ~= currX
            set(data.ddX, 'Value', idxX);
            needsUpdate = true;
        end
        if idxY ~= currY
            set(data.ddY, 'Value', idxY);
            needsUpdate = true;
        end

        if needsUpdate
            onUpdatePlot(hContainer, []);
        end
    end

    function highlightPoints(indices)
        % External function to highlight points with a halo
        data = hContainer.UserData;

        if isfield(data, 'hHighlight') && ~isempty(data.hHighlight)
            delete(data.hHighlight);
        end
        data.hHighlight = [];

        idxX  = get(data.ddX, 'Value');
        idxY  = get(data.ddY, 'Value');
        xName = data.numericVars{idxX};
        yName = data.numericVars{idxY};

        xData = data.tbl.(xName);
        yData = data.tbl.(yName);

        selX = xData(indices);
        selY = yData(indices);

        if isempty(selX), return; end

        hold(data.hAxScatter, 'on');
        data.hAxHighlight = plot(data.hAxScatter, selX, selY, 'o', ...
            'Color', 'r', 'LineWidth', 2, 'MarkerSize', 10, ...
            'PickableParts', 'none', 'HitTest', 'off', ...
            'DisplayName', 'Selected');
        hold(data.hAxScatter, 'off');

        data.hHighlight = data.hAxHighlight; % Fix naming
        hContainer.UserData = data;
    end

end
