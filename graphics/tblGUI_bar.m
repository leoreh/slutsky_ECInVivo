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
%       'yVar'        (string) Initial Y variable name (Numeric)
%       'xVar'        (string) Initial X variable name (Categorical)
%       'grpVar'      (string) Initial Group variable name (Categorical)
%       'mode'        (string) Initial view mode: 'bar' (default) or
%                     'points'. Toggle at runtime via the Mode button.
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
    warning('tblGUI_bar:invalidMode', ...
        'Invalid mode ''%s''; falling back to ''bar''.', modeIn);
    modeIn = 'bar';
end

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
guiData.tgMode = [];     % Mode toggle button handle

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
    'Callback', @onGrpChange);
currY = currY - ctlH - ctlGap*2;

% Statistic
uicontrol('Parent', hPanelControl, 'Style', 'text', 'String', 'Statistic:', ...
    'Units', 'normalized', 'Position', [ctlX, currY, ctlW, ctlH], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');
currY = currY - ctlH;

guiData.ddStatType = uicontrol('Parent', hPanelControl, 'Style', 'popupmenu', ...
    'String', {'Arithmetic', 'Geometric', 'Median'}, ...
    'Units', 'normalized', ...
    'Position', [ctlX, currY, ctlW, ctlH], ...
    'Value', 1, ...
    'Callback', @onUpdatePlot);
currY = currY - ctlH - ctlGap*2;

% Mode Toggle (Bars vs Points)
modeIsPoints = strcmp(modeIn, 'points');
if modeIsPoints
    tgLabel = 'Mode: Points';
else
    tgLabel = 'Mode: Bars';
end
guiData.tgMode = uicontrol('Parent', hPanelControl, 'Style', 'togglebutton', ...
    'String', tgLabel, 'Units', 'normalized', ...
    'Position', [ctlX, currY, ctlW, ctlH], ...
    'Value', double(modeIsPoints), ...
    'Callback', @onModeToggle);
currY = currY - ctlH - ctlGap*2;

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

        % Get Statistic
        idxStat = get(data.ddStatType, 'Value');
        statItems = get(data.ddStatType, 'String');
        statType = statItems{idxStat};

        % Aggregate Data for Bar Plot
        % Bar expects matrix nXCats x nGroups
        meanMat = nan(length(xCats), length(gCats));
        errLMat = nan(length(xCats), length(gCats));
        errHMat = nan(length(xCats), length(gCats));
        valsByCell = cell(length(xCats), length(gCats));    % Raw points per cell, for scatter mode

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
                    n = length(vals);
                    switch statType
                        case 'Arithmetic'
                            m = mean(vals);
                            s = std(vals) / sqrt(n);
                            meanMat(iX, iG) = m;
                            errLMat(iX, iG) = s;
                            errHMat(iX, iG) = s;

                        case 'Geometric'
                            vals(vals <= 0) = [];
                            if ~isempty(vals)
                                n = length(vals);
                                mLog = mean(log(vals));
                                sLog = std(log(vals)) / sqrt(n);
                                mGeo = exp(mLog);
                                meanMat(iX, iG) = mGeo;
                                errLMat(iX, iG) = mGeo - exp(mLog - sLog);
                                errHMat(iX, iG) = exp(mLog + sLog) - mGeo;
                            end

                        case 'Median'
                            m = median(vals);
                            q1 = prctile(vals, 25);
                            q3 = prctile(vals, 75);
                            iqrVal = q3 - q1;
                            notch = 1.57 * iqrVal / sqrt(n);

                            meanMat(iX, iG) = m;
                            errLMat(iX, iG) = notch;
                            errHMat(iX, iG) = notch;
                    end
                end
            end
        end

        % Plotting
        ax = data.hAx;
        % cla() skips children with HandleVisibility='off' (used below for
        % legend control), so use delete(allchild(...)) to force-clear all.
        delete(allchild(ax));
        hold(ax, 'on');

        % Per-group colors, shared by both render modes
        clrMat = lines(length(gCats));

        % Common grouped-bar geometry (also used to place jittered points)
        nbars = length(gCats);
        groupwidth = min(0.8, nbars/(nbars + 1.5));

        % Dispatch to the active render mode
        useScatter = (get(data.tgMode, 'Value') == 1);
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

    function onGrpChange(src, ~)
        data = hContainer.UserData;
        populateCheckboxes(data);
        onUpdatePlot(src, []);
    end

    function onFilterChange(~, ~)
        onUpdatePlot(hContainer, []);
    end

    function onModeToggle(src, ~)
        if get(src, 'Value') == 1
            set(src, 'String', 'Mode: Points');
        else
            set(src, 'String', 'Mode: Bars');
        end
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
