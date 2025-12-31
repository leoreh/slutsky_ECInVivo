function hFig = tblGUI_bar(tbl, varargin)
% TBLGUI_BAR Interactive bar plot with grouping and erro bars.
%
%   hFig = tblGUI_bar(tbl, ...) opens a GUI to visualize the table 'tbl'
%   as a bar chart.
%
%   INPUT:
%       tbl         (table) The data table to visualize.
%
%   OPTIONAL KEY-VALUE PAIRS:
%   OPTIONAL KEY-VALUE PAIRS:
%       'yVar'        (string) Initial Y variable name (Numeric)
%       'xVar'        (string) Initial X variable name (Categorical)
%       'grpVar'      (string) Initial Group variable name (Categorical)
%       'Parent'      (handle) Parent container.
%
%   See also: TBLGUI_SCATHIST, TBLGUI_XY

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addParameter(p, 'yVar', '', @(x) ischar(x) || isstring(x) || isempty(x));
addParameter(p, 'xVar', '', @(x) ischar(x) || isstring(x) || isempty(x));
addParameter(p, 'grpVar', '', @(x) ischar(x) || isstring(x) || isempty(x));
addParameter(p, 'Parent', [], @(x) isempty(x) || isgraphics(x));
parse(p, tbl, varargin{:});

tbl = p.Results.tbl;
yVarIn = p.Results.yVar;
xVarIn = p.Results.xVar;
grpVarIn = p.Results.grpVar;
hParent = p.Results.Parent;

%% ========================================================================
%  INITIALIZATION
%  ========================================================================

% Valid variables
allVars = tbl.Properties.VariableNames;

% Numeric Vars for Y-Axis
numericVars = allVars(varfun(@isnumeric, tbl, 'OutputFormat', 'uniform'));

% Categorical Vars for X-Axis and Group
% Include logical and string and categorical
catVars = allVars(varfun(@(x) iscategorical(x) || isstring(x) || islogical(x),...
    tbl, 'OutputFormat', 'uniform'));

% Sort for niceness
numericVars = sort(numericVars);
catVars = sort(catVars);

if isempty(numericVars)
    error('Input table must contain at least one numeric variable.');
end

% Defaults
curY = numericVars{1};
curX = '';
if ~isempty(catVars), curX = catVars{1}; end
curGrp = 'None';

% Apply overrides
if ~isempty(yVarIn) && ismember(yVarIn, numericVars), curY = yVarIn; end
if ~isempty(xVarIn) && ismember(xVarIn, catVars), curX = xVarIn; end
if ~isempty(grpVarIn) && ismember(grpVarIn, catVars), curGrp = grpVarIn; end

% Figure Setup
if isempty(hParent)
    hContainer = figure('Name', 'Bar Plot Visualizer', 'NumberTitle', 'off', ...
        'Units', 'pixels', 'Position', [100, 100, 1000, 700]);
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
% Prepend 'None' for Group
guiData.grpVars = [{'None'}, catVars];
guiData.chkGrp = [];     % Checkbox handles

%% ========================================================================
%  LAYOUT
%  ========================================================================

% Side Panel for Controls (Left) like tblGUI_xy
panelW = 0.2; % Slightly wider for readability
hPanelControl = uipanel('Parent', hContainer, 'Units', 'normalized', ...
    'Position', [0, 0, panelW, 1]);

% Main Plotting Area (Right)
hPanelPlot = uipanel('Parent', hContainer, 'Units', 'normalized', ...
    'Position', [panelW, 0, 1-panelW, 1], 'BorderType', 'none');

guiData.hAx = axes('Parent', hPanelPlot, 'Position', [0.1, 0.1, 0.85, 0.8]);

% --- Controls ---
ctlH = 0.04;
ctlGap = 0.01;
currY = 0.95;
ctlW = 0.9;
ctlX = 0.05;

% Y Variable (Response)
uicontrol('Parent', hPanelControl, 'Style', 'text', 'String', 'Y Variable (Numeric):', ...
    'Units', 'normalized', 'Position', [ctlX, currY, ctlW, ctlH], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');
currY = currY - ctlH;

guiData.ddY = uicontrol('Parent', hPanelControl, 'Style', 'popupmenu', ...
    'String', numericVars, 'Units', 'normalized', ...
    'Position', [ctlX, currY, ctlW, ctlH], ...
    'Callback', @onUpdatePlot);
currY = currY - ctlH - ctlGap*2;

% X Variable (Category)
uicontrol('Parent', hPanelControl, 'Style', 'text', 'String', 'X Variable (Cat):', ...
    'Units', 'normalized', 'Position', [ctlX, currY, ctlW, ctlH], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');
currY = currY - ctlH;

guiData.ddX = uicontrol('Parent', hPanelControl, 'Style', 'popupmenu', ...
    'String', catVars, 'Units', 'normalized', ...
    'Position', [ctlX, currY, ctlW, ctlH], ...
    'Callback', @onUpdatePlot);
currY = currY - ctlH - ctlGap*2;

% Group Variable (Splits bars)
uicontrol('Parent', hPanelControl, 'Style', 'text', 'String', 'Group Variable (Cat):', ...
    'Units', 'normalized', 'Position', [ctlX, currY, ctlW, ctlH], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');
currY = currY - ctlH;

guiData.ddGrp = uicontrol('Parent', hPanelControl, 'Style', 'popupmenu', ...
    'String', guiData.grpVars, 'Units', 'normalized', ...
    'Position', [ctlX, currY, ctlW, ctlH], ...
    'Position', [ctlX, currY, ctlW, ctlH], ...
    'Callback', @onGrpChange);

% Panel for Checkboxes
guiData.pnlGrp = uipanel('Parent', hPanelControl, 'BorderType', 'none', ...
    'Units', 'normalized', 'Position', [ctlX, 0.05, ctlW, currY - 0.06]);



% Set Init Values
set(guiData.ddY, 'Value', find(strcmp(numericVars, curY)));
idxX = find(strcmp(catVars, curX));
if isempty(idxX) && ~isempty(catVars), idxX = 1; end
if ~isempty(idxX), set(guiData.ddX, 'Value', idxX); end

idxGrp = find(strcmp(guiData.grpVars, curGrp));
if isempty(idxGrp), idxGrp = 1; end % Default None
set(guiData.ddGrp, 'Value', idxGrp);

hContainer.UserData = guiData;
guiData.pnlGrp.Parent = hPanelControl; % Ensure parentage just in case, though handled in constructor

% Initial Population & Plot
onGrpChange(hContainer, []);

%% ========================================================================
%  CALLBACKS
%  ========================================================================

    function onUpdatePlot(~, ~)
        data = hContainer.UserData;
        tblIn = data.tbl;

        idxY = get(data.ddY, 'Value');
        idxX = get(data.ddX, 'Value');
        idxGrp = get(data.ddGrp, 'Value');

        yName = data.numericVars{idxY};

        if isempty(data.catVars)
            % Should not happen given init check, but safety
            xlabel(data.hAx, 'No Categorical Vars');
            return;
        end
        xName = data.catVars{idxX};

        grpName = data.grpVars{idxGrp};

        % Data Prep
        yData = tblIn.(yName);
        xData = tblIn.(xName);

        % Ensure Categorical
        if ~iscategorical(xData), xData = categorical(xData); end
        xCats = categories(xData);
        % Only keep used categories
        xCats = xCats(ismember(xCats, unique(xData(~isnan(yData)))));

        hasGrp = ~strcmp(grpName, 'None');

        if hasGrp
            grpData = tblIn.(grpName);
            if ~iscategorical(grpData), grpData = categorical(grpData); end

            % Initial categories from data
            gCats = categories(grpData);
            gCats = gCats(ismember(gCats, unique(grpData(~isnan(yData)))));

            % Filter by Checkboxes
            if ~isempty(data.chkGrp)
                validH = isgraphics(data.chkGrp);
                if any(validH)
                    selectedIdx = arrayfun(@(x) get(x, 'Value'), data.chkGrp(validH));
                    allCats = arrayfun(@(x) string(get(x, 'String')), data.chkGrp(validH));
                    activeCats = cellstr(allCats(logical(selectedIdx)));

                    gCats = intersect(gCats, activeCats, 'stable');
                end
            end
        else
            gCats = {'All'};
        end

        % Aggregate Data for Bar Plot
        % Bar expects matrix nXCats x nGroups
        meanMat = nan(length(xCats), length(gCats));
        semMat = nan(length(xCats), length(gCats));

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
                    meanMat(iX, iG) = mean(vals);
                    semMat(iX, iG) = std(vals) / sqrt(length(vals));
                end
            end
        end

        % Plotting
        ax = data.hAx;
        cla(ax);
        hold(ax, 'on');

        b = bar(ax, meanMat, 'grouped');

        % Error Bars
        % Calculation of error bar positions for grouped bars
        ngroups = size(meanMat, 1);
        nbars = size(meanMat, 2);

        % Calculate the center of each bar
        groupwidth = min(0.8, nbars/(nbars + 1.5));

        for i = 1:nbars
            % Based on bar documentation logic for centers
            x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);

            errorbar(ax, x, meanMat(:,i), semMat(:,i), 'k', 'linestyle', 'none');
        end

        hold(ax, 'off');

        % Aesthetics
        set(ax, 'XTick', 1:length(xCats), 'XTickLabel', xCats);
        ylabel(ax, yName, 'Interpreter', 'none');
        xlabel(ax, xName, 'Interpreter', 'none');

        if hasGrp
            legend(ax, gCats, 'Location', 'best', 'Interpreter', 'none');
            title(ax, sprintf('%s by %s (grouped by %s)', yName, xName, grpName), 'Interpreter', 'none');
        else
            title(ax, sprintf('%s by %s', yName, xName), 'Interpreter', 'none');
        end

        grid(ax, 'on');
    end

    function onGrpChange(src, ~)
        data = hContainer.UserData;
        populateCheckboxes(data);
        onUpdatePlot(src, []);
    end

    function onFilterChange(~, ~)
        onUpdatePlot(hContainer, []);
    end

    function populateCheckboxes(data)
        delete(data.pnlGrp.Children);
        idxGrp = get(data.ddGrp, 'Value');
        grpItems = get(data.ddGrp, 'String');
        grpName = grpItems{idxGrp};

        if strcmp(grpName, 'None')
            data.chkGrp = [];
        else
            raw = data.tbl.(grpName);
            if islogical(raw), raw = categorical(raw); end
            if ~iscategorical(raw), raw = categorical(raw); end
            cats = categories(raw);
            cats = cats(ismember(cats, unique(raw))); % Show only existing

            nCats = length(cats);
            h = 1 / max(10, nCats + 1);
            w = 1;
            data.chkGrp = gobjects(1, nCats);
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

end
