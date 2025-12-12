function hFig = tblGUI_scatHist(tbl, varargin)
% TBLGUI_SCATHIST Interactive scatter plot with marginal histograms and grouping.
%
%   tblGUI_scatHist(tbl, ...) opens a GUI to visualize the table 'tbl'.
%   Features include variable selection, quadratic size scaling, dynamic
%   grouping, and region-based group assignment.
%
%   INPUT:
%       tbl         (table) The data table to visualize.
%
%   OPTIONAL KEY-VALUE PAIRS:
%       'cfg'         (struct) Initial configuration structure with fields:
%                     - xVar:   (string) Initial X variable name
%                     - yVar:   (string) Initial Y variable name
%                     - szVar:  (string) Initial Size variable name
%                     - grpVar: (string) Initial Group variable name
%                     - clr:    (m x 3)  RGB color matrix for groups
%                     - alpha:  (scalar) Marker transparency (0-1)
%
%       'varsExclude' (cell) List of variable names to exclude from dropdowns.
%                     Default: {'UnitID', 'Name', 'Mouse', 'File'}.
%
%   EXAMPLE:
%       cfg = struct('xVar','tp', 'yVar','mfr', 'szVar','pVal', 'grpVar','initType');
%       tblGUI_scatHist(fetTbl, 'cfg', cfg);
%

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'tbl', @istable);
addOptional(p, 'varsExclude', {'UnitID', 'Name', 'Mouse', 'File'}, @iscell);
addOptional(p, 'cfg', struct(), @isstruct);
addParameter(p, 'Parent', [], @(x) isempty(x) || isgraphics(x));
addParameter(p, 'SelectionCallback', [], @(x) isempty(x) || isa(x, 'function_handle'));
parse(p, tbl, varargin{:});

tbl = p.Results.tbl;
varsExclude = p.Results.varsExclude;
cfg = p.Results.cfg;
hParent = p.Results.Parent;
selCbk = p.Results.SelectionCallback;

%% ========================================================================
%  INITIALIZATION
%  ========================================================================

% Valid variables for X/Y/Size (Numeric) and Group (Categorical/String)
allVars = tbl.Properties.VariableNames;
numericVars = allVars(varfun(@isnumeric, tbl, 'OutputFormat', 'uniform'));
numericVars = setdiff(numericVars, varsExclude);

% Ensure Group variables are categorical. If logical or string, let's treat
% them as potential grouping vars, but for the dropdown we prefer
% categorical.
catVars = allVars(varfun(@(x) iscategorical(x) || isstring(x) || islogical(x),...
    tbl, 'OutputFormat', 'uniform'));
catVars = setdiff(catVars, varsExclude);

if isempty(numericVars)
    error('Input table must contain at least one numeric variable.');
end

% Defaults (Prioritize cfg, then heuristics)
curX = numericVars{1};
curY = numericVars{min(2, length(numericVars))};
curSize = 'None';
curGrp = '';
if ~isempty(catVars), curGrp = catVars{1}; end
defAlpha = 0.6;
defClr = [];

% Apply cfg overrides if present
if isfield(cfg, 'xVar') && ismember(cfg.xVar, numericVars), curX = cfg.xVar; end
if isfield(cfg, 'yVar') && ismember(cfg.yVar, numericVars), curY = cfg.yVar; end
if isfield(cfg, 'szVar') && ismember(cfg.szVar, numericVars), curSize = cfg.szVar; end
if isfield(cfg, 'grpVar') && ismember(cfg.grpVar, catVars), curGrp = cfg.grpVar; end
if isfield(cfg, 'alpha') && isnumeric(cfg.alpha), defAlpha = cfg.alpha; end
if isfield(cfg, 'clr') && size(cfg.clr,2)==3, defClr = cfg.clr; end

% Figure Setup
% Figure Setup
if isempty(hParent)
    hContainer = figure('Name', 'Table Visualizer', 'NumberTitle', 'off', ...
        'Position', [100, 100, 1000, 700], 'MenuBar', 'none', 'ToolBar', 'figure');
    hFig = hContainer;
else
    hContainer = hParent;
    hFig = ancestor(hContainer, 'figure');
end

% Data storage in figure handle for callbacks
guiData = struct();
guiData.tbl = tbl;
guiData.numericVars = numericVars;
guiData.catVars = catVars;
guiData.pointHandles = []; % To store scatter handles
guiData.polyRoi = [];      % To store polygon ROI
guiData.defAlpha = defAlpha;
guiData.defClr = defClr;
guiData.selCbk = selCbk;

%% ========================================================================
%  LAYOUT
%  ========================================================================

% Margins and spacing
panelW = 0.2; % Width of control panel
marg = 0.05;
bottomMarg = 0.1; % Increased for X label

% Axes positions (Scatter, HistX, HistY)
% [left bottom width height]
% Main Scatter: Bottom-Left (relative to plot area)
% Adjusted height to match control panel top alignment more closely
scatterW = 0.55;
scatterH = 0.55;
startX   = panelW + marg;
startY   = bottomMarg;

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
linkaxes([guiData.hAxScatter, guiData.hAxHistY], 'y');

%% ========================================================================
%  CONTROLS
%  ========================================================================

% Controls start from top aligned with HistX top approximately
% HistX top is at: startY + scatterH + 0.02 + 0.15 = ~0.82 if default values
% Let's position controls from top down.
ctlTop = axHistXPos(2) + axHistXPos(4); % Top of Hist X
ctlH = 0.04;
ctlGap = 0.01;
ctlW = panelW - 0.02;
ctlX = 0.01;

% X Var
uicontrol('Parent', hContainer, 'Style', 'text', 'String', 'X Variable:', ...
    'Units', 'normalized', 'Position', [ctlX, ctlTop - ctlH, ctlW, ctlH*0.7], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');
guiData.ddX = uicontrol('Parent', hContainer, 'Style', 'popupmenu', ...
    'String', numericVars, 'Units', 'normalized', ...
    'Position', [ctlX, ctlTop - 2*ctlH, ctlW, ctlH], ...
    'Callback', @onUpdatePlot);

% Y Var
currYTop = ctlTop - 2*ctlH - ctlGap - ctlH;
uicontrol('Parent', hContainer, 'Style', 'text', 'String', 'Y Variable:', ...
    'Units', 'normalized', 'Position', [ctlX, currYTop, ctlW, ctlH*0.7], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');
guiData.ddY = uicontrol('Parent', hContainer, 'Style', 'popupmenu', ...
    'String', numericVars, 'Units', 'normalized', ...
    'Position', [ctlX, currYTop - ctlH, ctlW, ctlH], ...
    'Callback', @onUpdatePlot);

% Size Var
currSizeTop = currYTop - ctlH - ctlGap - ctlH;
uicontrol('Parent', hContainer, 'Style', 'text', 'String', 'Size Variable:', ...
    'Units', 'normalized', 'Position', [ctlX, currSizeTop, ctlW, ctlH*0.7], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');
guiData.ddSize = uicontrol('Parent', hContainer, 'Style', 'popupmenu', ...
    'String', [{'None'}, numericVars], 'Units', 'normalized', ...
    'Position', [ctlX, currSizeTop - ctlH, ctlW, ctlH], ...
    'Callback', @onUpdatePlot);

% Group Var
currGrpTop = currSizeTop - ctlH - ctlGap - ctlH;
uicontrol('Parent', hContainer, 'Style', 'text', 'String', 'Group Variable:', ...
    'Units', 'normalized', 'Position', [ctlX, currGrpTop, ctlW, ctlH*0.7], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');
guiData.ddGrp = uicontrol('Parent', hContainer, 'Style', 'popupmenu', ...
    'String', [{'None'}, catVars], 'Units', 'normalized', ...
    'Position', [ctlX, currGrpTop - ctlH, ctlW, ctlH], ...
    'Callback', @onUpdatePlot);

% Set initial dropdown values
set(guiData.ddX, 'Value', find(strcmp(numericVars, curX)));
set(guiData.ddY, 'Value', find(strcmp(numericVars, curY)));
% Set Size (Handle 'None' at index 1)
idxSz = find(strcmp(numericVars, curSize));
if isempty(idxSz), set(guiData.ddSize, 'Value', 1); % None
else, set(guiData.ddSize, 'Value', idxSz + 1); end

% Set Group (Handle 'None' at index 1)
idxG = find(strcmp(catVars, curGrp));
if isempty(idxG), set(guiData.ddGrp, 'Value', 1); % None
else, set(guiData.ddGrp, 'Value', idxG + 1); end


% -- Buttons --

% Position interactive buttons lower down
guiData.btnSelect = uicontrol('Parent', hContainer, 'Style', 'pushbutton', ...
    'String', 'Select Group', ...
    'Units', 'normalized', 'Position', [ctlX, 0.40, ctlW, 0.05], ...
    'Callback', @onSelectRegion, 'FontWeight', 'bold');

guiData.btnSelectDot = uicontrol('Parent', hContainer, 'Style', 'pushbutton', ...
    'String', 'Select Dot', ...
    'Units', 'normalized', 'Position', [ctlX, 0.34, ctlW, 0.05], ...
    'Callback', @onSelectDot, 'FontWeight', 'bold');

uicontrol('Parent', hContainer, 'Style', 'text', ...
    'String', 'Instruction: Use "Select Group" (Polygon) or "Select Dot" (Double-click to assign) to update groups.', ...
    'Units', 'normalized', 'Position', [ctlX, 0.25, ctlW, 0.08], ...
    'HorizontalAlignment', 'left', 'ForegroundColor', [0.4 0.4 0.4], ...
    'FontSize', 8);

guiData.btnSave = uicontrol('Parent', hContainer, 'Style', 'pushbutton', ...
    'String', 'Save Table to Workspace', ...
    'Units', 'normalized', 'Position', [ctlX, 0.10, ctlW, 0.06], ...
    'Callback', @onSaveTable);

% Store guidata
hContainer.UserData = guiData;
% Return handle to container
hFig = hContainer;

% Initial Plot
onUpdatePlot(hContainer, []);

%% ====================================================================
%  CALLBACKS
%  ====================================================================

    function onUpdatePlot(~, ~)
        data = hContainer.UserData;
        tbl = data.tbl;

        % Get selections
        idxX = get(data.ddX, 'Value');
        idxY = get(data.ddY, 'Value');
        idxSize = get(data.ddSize, 'Value');
        idxGrp = get(data.ddGrp, 'Value');

        xName = data.numericVars{idxX};
        yName = data.numericVars{idxY};

        % Size logic
        sizeName = '';
        szData = [];
        if idxSize > 1
            sizeVarItems = get(data.ddSize, 'String');
            sizeName = sizeVarItems{idxSize};
            rawSz = tbl.(sizeName);
            rawSz = rawSz + abs(min(rawSz));

            % Robust normalization + Quadratic Skew
            lims = prctile(rawSz, [5 95]);
            minSz = lims(1);
            maxSz = lims(2);

            if maxSz > minSz
                % Clip values
                clippedSz = max(min(rawSz, maxSz), minSz);

                % Linear 0-1
                normSz = (clippedSz - minSz) / (maxSz - minSz);

                % Quadratic transform: (norm^2)
                % This compresses small values towards 0, keeping them small.
                % Large values stay large.
                normSq = normSz .^ 5;

                % Map to Size Range: 10 to 150
                % Small dots (background) will be ~10-20
                szData = 20 + 80 * normSq;
            else
                szData = repmat(20, length(rawSz), 1);
            end
        else
            szData = repmat(20, height(tbl), 1);
        end

        % Group logic
        grpName = '';
        groups = ones(height(tbl), 1); % Default single group
        grpLabels = {'All'};
        isGrpActive = false;

        grpItems = get(data.ddGrp, 'String');
        if idxGrp > 1 % "None" is 1
            grpName = grpItems{idxGrp};
            rawGrp = tbl.(grpName);

            % Handle different group types (categorical, logical, etc)
            if islogical(rawGrp)
                rawGrp = categorical(rawGrp);
            elseif ~iscategorical(rawGrp)
                rawGrp = categorical(rawGrp); % Force categorical for plotting
            end

            groups = rawGrp;
            grpLabels = categories(groups);
            isGrpActive = true;
        else
            isGrpActive = false;
        end

        % Get Data
        xData = tbl.(xName);
        yData = tbl.(yName);

        % Clear Axes
        cla(data.hAxScatter);
        cla(data.hAxHistX);
        cla(data.hAxHistY);

        hold(data.hAxScatter, 'on');
        hold(data.hAxHistX, 'on');
        hold(data.hAxHistY, 'on');

        % Plotting Loop

        % Determine Colors
        if ~isempty(data.defClr) && size(data.defClr,1) >= length(grpLabels)
            colors = data.defClr;
        else
            colors = lines(length(grpLabels));
        end

        for iG = 1:length(grpLabels)
            if isGrpActive
                currCat = grpLabels{iG};
                idx = (groups == currCat);
            else
                idx = true(height(tbl), 1);
            end

            if ~any(idx), continue; end

            xG = xData(idx);
            yG = yData(idx);
            szG = szData(idx);

            % Safety check for colors index
            cIdx = mod(iG-1, size(colors,1)) + 1;
            cG = colors(cIdx, :);

            % Count points
            nPoints = sum(idx);
            lbl = sprintf('%s (n=%d)', string(grpLabels{iG}), nPoints);

            % Scatter
            % Storing handle for legend
            scatter(data.hAxScatter, xG, yG, szG, cG, 'filled', ...
                'DisplayName', lbl, ...
                'MarkerFaceAlpha', data.defAlpha, ...
                'HitTest', 'off', 'PickableParts', 'none'); % Optimization

            % Hist X
            histogram(data.hAxHistX, xG, 30, 'FaceColor', cG, ...
                'EdgeColor', 'none', 'Normalization', 'count', 'FaceAlpha', 0.5);

            % Hist Y (Orientation horizontal)
            histogram(data.hAxHistY, yG, 30, 'FaceColor', cG, ...
                'EdgeColor', 'none', 'Normalization', 'count', 'FaceAlpha', 0.5, ...
                'Orientation', 'horizontal');
        end

        % Labels & Aesthetics
        xlabel(data.hAxScatter, xName, 'Interpreter', 'none');
        ylabel(data.hAxScatter, yName, 'Interpreter', 'none');

        if isGrpActive
            legend(data.hAxScatter, 'Location', 'best', 'Interpreter', 'none');
        end

        % Update Axis Limits (Fit to data with padding)
        validX = xData(~isnan(xData) & ~isinf(xData));
        validY = yData(~isnan(yData) & ~isinf(yData));

        if isempty(validX), xLim = [0 1]; else, xLim = [min(validX), max(validX)]; end
        if isempty(validY), yLim = [0 1]; else, yLim = [min(validY), max(validY)]; end

        % Add 5% padding
        dx = diff(xLim); if dx==0, dx=1; end
        dy = diff(yLim); if dy==0, dy=1; end

        xlim(data.hAxScatter, [xLim(1)-0.05*dx, xLim(2)+0.05*dx]);
        ylim(data.hAxScatter, [yLim(1)-0.05*dy, yLim(2)+0.05*dy]);

        grid(data.hAxScatter, 'on');
        data.hAxHistX.XAxis.Visible = 'off';
        data.hAxHistX.YAxis.Visible = 'off';
        data.hAxHistY.XAxis.Visible = 'off';
        data.hAxHistY.YAxis.Visible = 'off';

        linkaxes([data.hAxScatter, data.hAxHistX], 'x');
        linkaxes([data.hAxScatter, data.hAxHistY], 'y');

        hold(data.hAxScatter, 'off');
        hold(data.hAxHistX, 'off');
        hold(data.hAxHistY, 'off');

        % Update Enable state of "Assign Group"
        isReady = isGrpActive;
        if isReady
            set(data.btnSelect, 'Enable', 'on', 'String', 'Select Group');
            set(data.btnSelectDot, 'Enable', 'on', 'String', 'Select Dot');
        else
            set(data.btnSelect, 'Enable', 'off', 'String', 'Select Group');
            set(data.btnSelectDot, 'Enable', 'off', 'String', 'Select Dot');
        end
    end

    function onSelectRegion(~, ~)
        data = hContainer.UserData;
        ax = data.hAxScatter;
        roi = drawpolygon(ax);

        if isempty(roi.Position)
            delete(roi); return;
        end

        % Standard processing...
        % Standard processing...
        processSelection(hContainer, data, roi, false);
    end

    function onSelectDot(~, ~)
        data = hContainer.UserData;
        ax = data.hAxScatter;

        % Wait for user to create a point
        delete(findobj(ax, 'Type', 'images.roi.Point')); % Clear previous points
        hPoint = drawpoint(ax, 'Label', 'Target');

        if isempty(hPoint.Position)
            delete(hPoint); return;
        end

        % Snap to nearest point & Highlight
        snapToData(hPoint, data);

        % Add listener for double click
        addlistener(hPoint, 'ROIClicked', @(roiSrc, evt) onDotDoubleClick(roiSrc, evt, hContainer));

        % Add listener for move to re-snap
        addlistener(hPoint, 'MovingROI', @(roiSrc, evt) snapToData(roiSrc, data));

        % Instruction
        title(ax, 'Double-click to assign group. Drag to browse traces.', 'Color', 'r');
    end

    function snapToData(hPoint, data)
        pos = hPoint.Position;

        idxX = get(data.ddX, 'Value');
        idxY = get(data.ddY, 'Value');
        xName = data.numericVars{idxX};
        yName = data.numericVars{idxY};

        xData = data.tbl.(xName);
        yData = data.tbl.(yName);

        % Simple Euclidean distance (Consider normalizing if axes differ widely)
        % For better UX, normalize by axis range
        ax = data.hAxScatter;
        xlims = ax.XLim; ylims = ax.YLim;

        xNorm = (xData - xlims(1)) / diff(xlims);
        yNorm = (yData - ylims(1)) / diff(ylims);
        pNorm = [(pos(1)-xlims(1))/diff(xlims), (pos(2)-ylims(1))/diff(ylims)];

        distSq = (xNorm - pNorm(1)).^2 + (yNorm - pNorm(2)).^2;
        [~, minIdx] = min(distSq);

        newPos = [xData(minIdx), yData(minIdx)];

        % Update position only if effectively different to avoid loop
        if sum((hPoint.Position - newPos).^2) > 1e-10
            hPoint.Position = newPos;
        end

        % Store index in UserData of the ROI
        hPoint.UserData = minIdx;

        % Trigger Callback immediately for single click/drag highlighting
        inPoints = false(height(data.tbl), 1);
        inPoints(minIdx) = true;
        if isfield(data, 'selCbk') && ~isempty(data.selCbk)
            data.selCbk(inPoints);
        end
    end

    function onDotDoubleClick(hPoint, evt, ~)
        if strcmp(evt.SelectionType, 'double')
            data = hContainer.UserData;
            idx = hPoint.UserData;

            % Process single point selection logic
            processSelectionIdx(hContainer, data, idx);

            delete(hPoint);
            % Restore title
            title(data.hAxScatter, '');
        end
    end

    function processSelection(src, data, roi, isSinglePoint)
        idxX = get(data.ddX, 'Value');
        idxY = get(data.ddY, 'Value');
        xName = data.numericVars{idxX};
        yName = data.numericVars{idxY};

        xData = data.tbl.(xName);
        yData = data.tbl.(yName);

        inPoints = inpolygon(xData, yData, roi.Position(:,1), roi.Position(:,2));
        nSelected = sum(inPoints);

        if nSelected == 0
            msgbox('No points selected.', 'Selection');
            delete(roi);
            return;
        end

        % Only trigger callback if NOT polygon selection (as per user request)
        if isSinglePoint
            % Trigger selCbk logic is removed from here and handled in snapToData or processSelectionIdx
            % Actually, user said polygon shouldn't affect traces.
            % So we simply do NOT call selCbk here.
            if ~isempty(data.selCbk)
                data.selCbk(inPoints);
            end
        end

        assignGroup(src, data, inPoints, nSelected);
        delete(roi);
    end

    function processSelectionIdx(src, data, idx)
        inPoints = false(height(data.tbl), 1);
        inPoints(idx) = true;

        if ~isempty(data.selCbk)
            data.selCbk(inPoints);
        end

        assignGroup(src, data, inPoints, 1);
    end

    function assignGroup(src, data, inPoints, nSelected)
        % Prompt for Group Assignment
        idxGrp = get(data.ddGrp, 'Value');
        grpItems = get(data.ddGrp, 'String');
        grpName = grpItems{idxGrp};

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
                if strcmpi(selectedCat, 'true'), val = true; else, val = false; end
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

end
