function hFig = tblGUI_raster(dataTbl, varargin)
% TBLGUI_RASTER Interactive raster plot for Table variables.
%
%   hFig = TBLGUI_RASTER(dataTbl, ...) plots spike rasters from 'dataTbl'.
%   - Top Panel: Overview (Time in Hours). Click to zoom.
%   - Bottom Panel: Zoomed view (Time in Seconds).
%   - Side Panel: Filter by categorical variables.
%
%   INPUTS:
%       dataTbl   (table)  Table containing spike times and metadata.
%
%   OPTIONAL KEY-VALUE PAIRS:
%       'timesVar'         (char) Column name for spike times {'spktimes'}.
%       'grpVar'           (char) Initial grouping/filtering variable.
%       'grpVal'           (char/cell) Initial value to filter by.
%       'Parent'           (handle) Parent container.
%
%   OUTPUT:
%       hFig      (handle) Figure handle.
%
%   See also: PLOT_RASTER, TBLGUI_XY

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'dataTbl', @istable);
addParameter(p, 'timesVar', 'spktimes', @ischar);
addParameter(p, 'grpVar', [], @(x) ischar(x) || isempty(x));
addParameter(p, 'grpVal', [], @(x) ischar(x) || iscell(x) || isstring(x) || isempty(x));
addParameter(p, 'Parent', [], @(x) isempty(x) || isgraphics(x));
parse(p, dataTbl, varargin{:});

timesVar = p.Results.timesVar;
initialGrpVar = p.Results.grpVar;
initialGrpVal = p.Results.grpVal;
hParent = p.Results.Parent;

if ~ismember(timesVar, dataTbl.Properties.VariableNames)
    error('Variable "%s" not found in table.', timesVar);
end


%% ========================================================================
%  INITIALIZATION
%  ========================================================================

% Identify potential grouping variables
allVars = dataTbl.Properties.VariableNames;
catVars = allVars(varfun(@(x) iscategorical(x) || isstring(x) || islogical(x),...
    dataTbl, 'OutputFormat', 'uniform'));
catVars = [{'None'}, catVars];

% Auto-select grpVar if valid
if ~isempty(initialGrpVar) && ismember(initialGrpVar, catVars)
    grpVar = initialGrpVar;
else
    grpVar = 'None';
end

% Figure Setup
if isempty(hParent)
    hContainer = figure('Name', 'Table Raster GUI', 'NumberTitle', 'off', ...
        'Units', 'pixels', 'Position', [100, 100, 1200, 800], 'Color', 'w');
    hFig = hContainer;
else
    hContainer = hParent;
    hFig = ancestor(hContainer, 'figure');
end

% Store GUI Data
guiData = struct();
guiData.dataTbl = dataTbl;
guiData.timesVar = timesVar;
guiData.catVars = catVars;
guiData.activeIndices = true(height(dataTbl), 1); % Currently visible rows
guiData.currentCenter = 0; % Center of zoomed view (seconds)
guiData.winWidth = 60;      % Default zoom width (seconds)
guiData.initialGrpVal = initialGrpVal; % Store for first filter

guiData.hAx = gobjects(2, 1);
guiData.chkGrpBy = [];


%% ========================================================================
%  LAYOUT
%  ========================================================================

% Side Panel for Controls (Left)
panelW = 0.15;
hPanel = uipanel('Parent', hContainer, 'Units', 'normalized', ...
    'Position', [0, 0, panelW, 1]);

% Main Plotting Area (Right)
hPanelRight = uipanel('Parent', hContainer, 'Units', 'normalized', ...
    'Position', [panelW, 0, 1-panelW, 1], 'BorderType', 'none', ...
    'BackgroundColor', 'w');

guiData.tlo = tiledlayout(hPanelRight, 2, 1, ...
    'TileSpacing', 'compact', 'Padding', 'compact');


% --- Controls in Side Panel ---
ctlH = 0.03;
ctlGap = 0.01;
currY = 0.95;

% 1. FILTER BY GROUP
uicontrol('Parent', hPanel, 'Style', 'text', 'String', 'Filter By:', ...
    'Units', 'normalized', 'Position', [0.05, currY, 0.9, ctlH], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');
currY = currY - ctlH;

valGrp = find(strcmp(catVars, grpVar), 1);
if isempty(valGrp), valGrp = 1; end

guiData.ddGrpBy = uicontrol('Parent', hPanel, 'Style', 'popupmenu', ...
    'String', catVars, 'Value', valGrp, ...
    'Units', 'normalized', 'Position', [0.05, currY, 0.9, ctlH], ...
    'Callback', @onGrpByChange);
currY = currY - ctlH - ctlGap;

% Container for Checkboxes
guiData.pnlGrpBy = uipanel('Parent', hPanel, 'BorderType', 'none', ...
    'Units', 'normalized', 'Position', [0.05, 0.4, 0.9, currY - 0.4]);


% 2. ZOOM WINDOW
currY = 0.35;
uicontrol('Parent', hPanel, 'Style', 'text', 'String', 'Zoom Window:', ...
    'Units', 'normalized', 'Position', [0.05, currY, 0.9, ctlH], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');
currY = currY - ctlH;

winOpts = [10, 30, 60, 120, 300, 600, 1800, 3600];
winLbls = {'10 s', '30 s', '1 min', '2 min', '5 min', '10 min', '30 min', '1 hr'};
guiData.winOpts = winOpts;

guiData.ddWin = uicontrol('Parent', hPanel, 'Style', 'popupmenu', ...
    'String', winLbls, 'Value', 3, ... % Default 1 min
    'Units', 'normalized', 'Position', [0.05, currY, 0.9, ctlH], ...
    'Callback', @onWinChange);


% Store Data
hContainer.UserData = guiData;

% Initial State
onGrpByChange(hContainer, []); % Populates checkboxes and plots


%% ========================================================================
%  CALLBACKS
%  ========================================================================

    function onGrpByChange(~, ~)
        data = hContainer.UserData;
        idx = get(data.ddGrpBy, 'Value');
        varName = data.catVars{idx};

        % Populate Checkboxes
        delete(data.pnlGrpBy.Children);
        data.chkGrpBy = [];

        if ~strcmp(varName, 'None')
            raw = data.dataTbl.(varName);
            if islogical(raw), raw = categorical(raw); end
            if ~iscategorical(raw), raw = categorical(raw); end
            cats = categories(raw);
            cats = cats(ismember(cats, unique(raw)));

            % Determine initial selection
            if ~isempty(data.initialGrpVal)
                % If user provided a specific value, select only that
                target = string(data.initialGrpVal);
                initVal = ismember(string(cats), target);
                % If no match found, default to all
                if sum(initVal) == 0
                    warning('grpVal "%s" not found in %s. Selecting all.', target, varName);
                    initVal = true(size(cats));
                end
                % Clear it so it only applies once
                data.initialGrpVal = [];
            else
                initVal = true(size(cats));
            end

            nCats = length(cats);
            h = 1 / max(10, nCats + 1);
            w = 0.9;
            for k = 1:nCats
                yPos = 1 - k*h;
                data.chkGrpBy(k) = uicontrol('Parent', data.pnlGrpBy, ...
                    'Style', 'checkbox', 'String', cats{k}, ...
                    'Units', 'normalized', 'Position', [0, yPos, w, h], ...
                    'Value', initVal(k), 'Callback', @onFilterChange);
            end
        end

        hContainer.UserData = data;
        updateActiveIndices();
        plotAll();
    end

    function onFilterChange(~, ~)
        updateActiveIndices();
        plotAll();
    end

    function onWinChange(~, ~)
        updateBottomPanel();
    end

    function onClickTop(~, event)
        % Update center position
        pt = event.IntersectionPoint;
        % Top plot is in Hours, convert to Seconds
        data = hContainer.UserData;
        data.currentCenter = pt(1) * 3600;
        hContainer.UserData = data;

        updateBottomPanel();
        drawFocusLine();
    end

    function updateActiveIndices()
        data = hContainer.UserData;
        idxVal = get(data.ddGrpBy, 'Value');
        varName = data.catVars{idxVal};

        if strcmp(varName, 'None')
            data.activeIndices = true(height(data.dataTbl), 1);
        else
            % Get allowed categories
            chk = data.chkGrpBy;
            if isempty(chk)
                % If no checkboxes (None), all active? No, checked above.
                % This block is reached if categories exist but checkboxes missing?
                data.activeIndices = true(height(data.dataTbl), 1);
            else
                areSel = arrayfun(@(x) get(x, 'Value'), chk);
                allCats = arrayfun(@(x) string(get(x, 'String')), chk);
                selCats = allCats(logical(areSel));

                raw = data.dataTbl.(varName);
                if islogical(raw), raw = categorical(raw); end
                if ~iscategorical(raw), raw = categorical(raw); end

                data.activeIndices = ismember(string(raw), selCats);
            end
        end
        hContainer.UserData = data;
    end

    function plotAll()
        data = hContainer.UserData;

        % Prepare Data
        spikes = data.dataTbl.(data.timesVar)(data.activeIndices);
        if ~iscell(spikes), spikes = {spikes}; end

        % Remove empty cells to avoid gaps in raster (User requested compact view)
        spikes = spikes(~cellfun('isempty', spikes));

        % --- Top Panel: Hours ---
        nexttile(data.tlo, 1);
        data.hAx(1) = gca;
        cla(data.hAx(1));

        % Convert spikes to hours for display
        spikesHrs = cellfun(@(x) x / 3600, spikes, 'UniformOutput', false);

        [~, hPlt] = plot_raster(spikesHrs, 'hAx', data.hAx(1), 'plotType', 'vertline');
        title(data.hAx(1), 'Overview (Hours) - Click to Focus');
        xlabel(data.hAx(1), 'Time (hr)');

        % Add interaction
        set(data.hAx(1), 'ButtonDownFcn', @onClickTop);
        set(hPlt, 'HitTest', 'off'); % Click passes to Axes
        set(hPlt, 'Color', [0.2 0.2 0.2 0.5]); % Slight transparency

        % --- Bottom Panel: Seconds ---
        nexttile(data.tlo, 2);
        data.hAx(2) = gca;
        cla(data.hAx(2));
        plot_raster(spikes, 'hAx', data.hAx(2), 'plotType', 'vertline');
        title(data.hAx(2), 'Zoomed View (Seconds)');
        xlabel(data.hAx(2), 'Time (s)');

        linkaxes(data.hAx, 'y');

        hContainer.UserData = data;
        updateBottomPanel();
        drawFocusLine();
    end

    function updateBottomPanel()
        data = hContainer.UserData;

        idxWin = get(data.ddWin, 'Value');
        width = data.winOpts(idxWin);

        center = data.currentCenter;

        xlim(data.hAx(2), [center - width/2, center + width/2]);
        drawFocusLine();
    end

    function drawFocusLine()
        data = hContainer.UserData;
        centerHr = data.currentCenter / 3600;

        hLine = findobj(data.hAx(1), 'Tag', 'FocusLine');
        if isempty(hLine)
            xline(data.hAx(1), centerHr, 'r-', 'LineWidth', 1.5, ...
                'Tag', 'FocusLine', 'HitTest', 'off');
        else
            set(hLine, 'Value', centerHr);
        end
    end

end
