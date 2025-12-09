function plot_tblGUI(tbl, varargin)
% PLOT_TBLGUI Interactive scatter plot with marginal histograms and grouping.
%
%   plot_tblGUI(tbl, ...) opens a GUI to visualize the table 'tbl'.
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
%       plot_tblGUI(fetTbl, 'cfg', cfg);
%

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
p = inputParser;
addRequired(p, 'tbl', @istable);
addOptional(p, 'varsExclude', {'UnitID', 'Name', 'Mouse', 'File'}, @iscell);
addOptional(p, 'cfg', struct(), @isstruct);
parse(p, tbl, varargin{:});

tbl = p.Results.tbl;
varsExclude = p.Results.varsExclude;
cfg = p.Results.cfg;

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
hFig = figure('Name', 'Table Visualizer', 'NumberTitle', 'off', ...
    'Position', [100, 100, 1000, 700], 'MenuBar', 'none', 'ToolBar', 'figure');

% Data storage in figure handle for callbacks
guiData = struct();
guiData.tbl = tbl;
guiData.numericVars = numericVars;
guiData.catVars = catVars;
guiData.pointHandles = []; % To store scatter handles
guiData.polyRoi = [];      % To store polygon ROI
guiData.defAlpha = defAlpha;
guiData.defClr = defClr;

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

guiData.hAxScatter = axes('Parent', hFig, 'Position', axScatterPos);
guiData.hAxHistX   = axes('Parent', hFig, 'Position', axHistXPos);
guiData.hAxHistY   = axes('Parent', hFig, 'Position', axHistYPos);

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
uicontrol('Parent', hFig, 'Style', 'text', 'String', 'X Variable:', ...
    'Units', 'normalized', 'Position', [ctlX, ctlTop - ctlH, ctlW, ctlH*0.7], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');
guiData.ddX = uicontrol('Parent', hFig, 'Style', 'popupmenu', ...
    'String', numericVars, 'Units', 'normalized', ...
    'Position', [ctlX, ctlTop - 2*ctlH, ctlW, ctlH], ...
    'Callback', @onUpdatePlot);

% Y Var
currYTop = ctlTop - 2*ctlH - ctlGap - ctlH;
uicontrol('Parent', hFig, 'Style', 'text', 'String', 'Y Variable:', ...
    'Units', 'normalized', 'Position', [ctlX, currYTop, ctlW, ctlH*0.7], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');
guiData.ddY = uicontrol('Parent', hFig, 'Style', 'popupmenu', ...
    'String', numericVars, 'Units', 'normalized', ...
    'Position', [ctlX, currYTop - ctlH, ctlW, ctlH], ...
    'Callback', @onUpdatePlot);

% Size Var
currSizeTop = currYTop - ctlH - ctlGap - ctlH;
uicontrol('Parent', hFig, 'Style', 'text', 'String', 'Size Variable:', ...
    'Units', 'normalized', 'Position', [ctlX, currSizeTop, ctlW, ctlH*0.7], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');
guiData.ddSize = uicontrol('Parent', hFig, 'Style', 'popupmenu', ...
    'String', [{'None'}, numericVars], 'Units', 'normalized', ...
    'Position', [ctlX, currSizeTop - ctlH, ctlW, ctlH], ...
    'Callback', @onUpdatePlot);

% Group Var
currGrpTop = currSizeTop - ctlH - ctlGap - ctlH;
uicontrol('Parent', hFig, 'Style', 'text', 'String', 'Group Variable:', ...
    'Units', 'normalized', 'Position', [ctlX, currGrpTop, ctlW, ctlH*0.7], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');
guiData.ddGrp = uicontrol('Parent', hFig, 'Style', 'popupmenu', ...
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
guiData.btnSelect = uicontrol('Parent', hFig, 'Style', 'pushbutton', ...
    'String', 'Select & Assign Group', ...
    'Units', 'normalized', 'Position', [ctlX, 0.40, ctlW, 0.06], ...
    'Callback', @onSelectRegion, 'FontWeight', 'bold');

uicontrol('Parent', hFig, 'Style', 'text', ...
    'String', 'Instruction: Click "Select", draw polygon, double-click to finish.', ...
    'Units', 'normalized', 'Position', [ctlX, 0.34, ctlW, 0.05], ...
    'HorizontalAlignment', 'left', 'ForegroundColor', [0.4 0.4 0.4], ...
    'FontSize', 8);

guiData.btnSave = uicontrol('Parent', hFig, 'Style', 'pushbutton', ...
    'String', 'Save Table to Workspace', ...
    'Units', 'normalized', 'Position', [ctlX, 0.10, ctlW, 0.06], ...
    'Callback', @onSaveTable);

% Store guidata
guidata(hFig, guiData);

% Initial Plot
onUpdatePlot(hFig, []);

%% ====================================================================
%  CALLBACKS
%  ====================================================================

    function onUpdatePlot(src, ~)
        data = guidata(src);
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

            % Scatter
            % Storing handle for legend
            scatter(data.hAxScatter, xG, yG, szG, cG, 'filled', ...
                'DisplayName', string(grpLabels{iG}), ...
                'MarkerFaceAlpha', data.defAlpha);

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
        if isGrpActive
            set(data.btnSelect, 'Enable', 'on', 'String', 'Select & Assign Group');
        else
            set(data.btnSelect, 'Enable', 'off', 'String', 'Select (Select Group Var First)');
        end
    end

    function onSelectRegion(src, ~)
        data = guidata(src);

        % Ensure we are targeting the scatter plot
        ax = data.hAxScatter;

        % Use drawpolygon (R2018b+)
        roi = drawpolygon(ax);

        if isempty(roi.Position)
            delete(roi);
            return;
        end

        % Find points inside
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

        % Prompt for Group Assignment
        idxGrp = get(data.ddGrp, 'Value');
        grpItems = get(data.ddGrp, 'String');
        grpName = grpItems{idxGrp};

        % We need the categories of the target group
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
                % Convert string selection back to logical
                if strcmpi(selectedCat, 'true'), val = true; else, val = false; end
                data.tbl.(grpName)(inPoints) = val;
            else
                data.tbl.(grpName)(inPoints) = selectedCat;
            end

            % Update GUI Data
            guidata(src, data);

            % Refresh Plot
            onUpdatePlot(src, []);

            fprintf('Updated %d points in variable "%s" to "%s".\n', ...
                nSelected, grpName, string(selectedCat));
        end

        delete(roi);
    end

    function onSaveTable(src, ~)
        data = guidata(src);
        % Assign to base workspace
        assignin('base', 'fetTbl_mod', data.tbl);
        msgbox('Table saved to workspace as "fetTbl_mod".', 'Saved');
    end

end
