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
        'Position', [100, 100, 1000, 700], 'MenuBar', 'none', 'ToolBar', 'figure');
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
guiData.chkGrp          = [];   % Checkbox handles

%% ========================================================================
%  LAYOUT
%  ========================================================================

% Margins and spacing
panelW      = 0.2;  % Width of control panel
marg        = 0.05;
bottomMarg  = 0.1;

% Main Scatter: Bottom-Left (relative to plot area)
scatterW    = 0.55;
scatterH    = 0.55;
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
ctlTop      = 0.95;
ctlH        = 0.04;
ctlGap      = 0.01;
ctlW        = panelW - 0.02;     % Full width
ctlW_Half   = ctlW / 2 - 0.005;  % Half width

ctlX        = 0.01;
ctlX_Right  = ctlX + ctlW_Half + 0.01;

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

% --- Y Variable ---
currYTop = ctlTop - 2*ctlH - ctlGap - ctlH;
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

% --- Size Variable ---
currSizeTop = currYTop - ctlH - ctlGap - ctlH;
uicontrol('Parent', hContainer, 'Style', 'text', 'String', 'Size Variable:', ...
    'Units', 'normalized', 'Position', [ctlX, currSizeTop, ctlW, ctlH*0.7], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');
guiData.ddSize = uicontrol('Parent', hContainer, 'Style', 'popupmenu', ...
    'String', [{'None'}, numericVars], 'Units', 'normalized', ...
    'Position', [ctlX, currSizeTop - ctlH, ctlW, ctlH], ...
    'Callback', @onUpdatePlot);

% --- Group Variable ---
currGrpTop = currSizeTop - ctlH - ctlGap - ctlH;
uicontrol('Parent', hContainer, 'Style', 'text', 'String', 'Group Variable:', ...
    'Units', 'normalized', 'Position', [ctlX, currGrpTop, ctlW, ctlH*0.7], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');
guiData.ddGrp = uicontrol('Parent', hContainer, 'Style', 'popupmenu', ...
    'String', [{'None'}, catVars], 'Units', 'normalized', ...
    'Position', [ctlX, currGrpTop - ctlH, ctlW, ctlH], ...
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

% --- Filter Panel ---
grpDDBottom = currGrpTop - ctlH;
panelTop    = grpDDBottom - 0.02;
panelBottom = 0.30;
panelH      = panelTop - panelBottom;

guiData.pnlGrp = uipanel('Parent', hContainer, 'BorderType', 'none', ...
    'Units', 'normalized', 'Position', [ctlX, panelBottom, ctlW, panelH]);

% --- Buttons ---
guiData.btnSelect = uicontrol('Parent', hContainer, 'Style', 'pushbutton', ...
    'String', 'Select Group', ...
    'Units', 'normalized', 'Position', [ctlX, 0.22, ctlW, 0.05], ...
    'Callback', @onSelectRegion, 'FontWeight', 'bold');

guiData.btnSelectDot = uicontrol('Parent', hContainer, 'Style', 'pushbutton', ...
    'String', 'Select Dot', ...
    'Units', 'normalized', 'Position', [ctlX, 0.16, ctlW, 0.05], ...
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

            for i = 1:nCats
                yPos = 1 - i*h;
                data.chkGrp(i) = uicontrol('Parent', data.pnlGrp, 'Style', 'checkbox', ...
                    'String', cats{i}, 'Units', 'normalized', ...
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

        xName = data.numericVars{idxX};
        yName = data.numericVars{idxY};

        % --- Size Logic (Quadratic Scaling) ---
        szData = [];
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
        isGrpActive = false;

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
        cla(data.hAxScatter);
        cla(data.hAxHistX);
        cla(data.hAxHistY);

        hold(data.hAxScatter, 'on');
        hold(data.hAxHistX, 'on');
        hold(data.hAxHistY, 'on');

        % Plotting Loop
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

            % Correlation Stats
            [rho, pval] = corr(xG, yG, 'Type', 'Spearman', 'Rows', 'complete');
            if pval < 0.0001, pStr = 'p<0.0001';
            else,             pStr = sprintf('p=%.4f', pval); end
            statsStr = sprintf('\\rho=%.2f, %s', rho, pStr);

            lbl = sprintf('%s (n=%d, %s)', string(currCat), nPoints, statsStr);

            % Scatter Plot
            scatter(data.hAxScatter, xG, yG, szG, cG, 'filled', ...
                'DisplayName', lbl, ...
                'MarkerFaceAlpha', data.defAlpha, ...
                'HitTest', 'off', 'PickableParts', 'none');

            % Histogram X
            plot_histPdf(data.hAxHistX, xG, xEdges, cG, scaleX, 'vertical');

            % Histogram Y (Horizontal)
            plot_histPdf(data.hAxHistY, yG, yEdges, cG, scaleY, 'horizontal');
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
            plot(data.hAxScatter, xLim, xLim, 'k--', 'LineWidth', 1, ...
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
    end

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
    end

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

    function plot_histPdf(ax, val, edges, c, scale, orient)
        % Helper to plot histogram with KDE overlay

        % Filter data
        val = val(~isnan(val) & ~isinf(val));
        if isempty(val), return; end

        % Histogram
        histogram(ax, val, 'BinEdges', edges, 'FaceColor', c, ...
            'EdgeColor', 'none', 'FaceAlpha', 0.2, 'Normalization', 'pdf', ...
            'DisplayStyle', 'bar', 'Orientation', orient);

        % KDE Overlay
        if length(val) < 2, return; end

        pts = linspace(min(edges), max(edges), 200);

        if strcmp(scale, 'log')
            val = val(val > 0);
            if length(val) < 2, return; end

            % Log-KDE: Estimate on log10(x), then transform PDF density back
            [f_log, ~] = ksdensity(log10(val), log10(pts));
            f = f_log ./ (pts * log(10));
        else
            % Linear-KDE
            [f, ~] = ksdensity(val, pts);
        end

        % Check for NaNs
        if all(isnan(f)), return; end

        % Plot Curve
        if strcmp(orient, 'vertical')
            plot(ax, pts, f, 'Color', c, 'LineWidth', 2);
        else
            plot(ax, f, pts, 'Color', c, 'LineWidth', 2);
        end
    end

end
