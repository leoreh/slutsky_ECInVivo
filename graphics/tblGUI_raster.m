function hFig = tblGUI_raster(dataTbl, varargin)
% TBLGUI_RASTER Interactive raster plot for Table variables.
%
%   hFig = TBLGUI_RASTER(dataTbl, ...) plots spike rasters from 'dataTbl'.
%   A single raster panel shows all spikes (black) with optional burst
%   spike overlay (red). Navigation uses editable text boxes for the
%   display window center and width.
%
%   - Side Panel: Filter by categorical variables, adjust view window.
%   - Main Panel: Raster with fixed typography (Arial, 10/12 pt).
%
%   INPUTS:
%       dataTbl   (table)  Table containing spike times and metadata.
%
%   OPTIONAL KEY-VALUE PAIRS:
%       'timesVar'  (char) Column name for spike times {'spktimes'}.
%       'brstVar'   (char) Column name for burst spike times {''}.
%                   If empty or not found, burst overlay is skipped.
%       'grpVar'    (char) Initial grouping/filtering variable.
%       'grpVal'    (char/cell) Initial value to filter by.
%       'timeLim'   (1x2 numeric) Clip spike times to [start, end] {[]}.
%                   Applied once at initialization for speed.
%       'Parent'    (handle) Parent container.
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
addParameter(p, 'brstVar', '', @ischar);
addParameter(p, 'grpVar', [], @(x) ischar(x) || isempty(x));
addParameter(p, 'grpVal', [], @(x) ischar(x) || iscell(x) || isstring(x) || isempty(x));
addParameter(p, 'timeLim', [], @isnumeric);
addParameter(p, 'Parent', [], @(x) isempty(x) || isgraphics(x));
parse(p, dataTbl, varargin{:});

timesVar       = p.Results.timesVar;
brstVar        = p.Results.brstVar;
initialGrpVar  = p.Results.grpVar;
initialGrpVal  = p.Results.grpVal;
timeLim        = p.Results.timeLim;
hParent        = p.Results.Parent;

if ~ismember(timesVar, dataTbl.Properties.VariableNames)
    error('Variable "%s" not found in table.', timesVar);
end

% Resolve burst variable
flgBrst = ~isempty(brstVar) && ...
    ismember(brstVar, dataTbl.Properties.VariableNames);

% Clip spike times to timeLim (applied once for speed)
if ~isempty(timeLim)
    dataTbl.(timesVar) = cellfun(@(x) x(x >= timeLim(1) & x <= timeLim(2)), ...
        dataTbl.(timesVar), 'UniformOutput', false);
    if flgBrst
        dataTbl.(brstVar) = cellfun(@(x) x(x >= timeLim(1) & x <= timeLim(2)), ...
            dataTbl.(brstVar), 'UniformOutput', false);
    end
end


%% ========================================================================
%  INITIALIZATION
%  ========================================================================

% Identify categorical variables for filtering
allVars = dataTbl.Properties.VariableNames;
catVars = allVars(varfun(@(x) iscategorical(x) || isstring(x) || islogical(x), ...
    dataTbl, 'OutputFormat', 'uniform'));
catVars = [{'None'}, catVars];

% Auto-select grpVar if valid
if ~isempty(initialGrpVar) && ismember(initialGrpVar, catVars)
    grpVar = initialGrpVar;
else
    grpVar = 'None';
end

% Data extent (used for navigation defaults)
allSpks = vertcat(dataTbl.(timesVar){:});
if isempty(allSpks)
    tMin = 0;  tMax = 1;
else
    tMin = min(allSpks);  tMax = max(allSpks);
end


%% ========================================================================
%  GUI SETUP
%  ========================================================================

if isempty(hParent)
    hContainer = figure('Name', 'Table Raster GUI', 'NumberTitle', 'off', ...
        'Units', 'pixels', 'Position', [100, 100, 1200, 700], 'Color', 'w');
    hFig = hContainer;
else
    hContainer = hParent;
    hFig = ancestor(hContainer, 'figure');
end

% GUI State
guiData = struct();
guiData.dataTbl        = dataTbl;
guiData.timesVar       = timesVar;
guiData.brstVar        = brstVar;
guiData.flgBrst        = flgBrst;
guiData.catVars        = catVars;
guiData.activeIndices  = true(height(dataTbl), 1);
guiData.tMin           = tMin;
guiData.tMax           = tMax;
guiData.winCenter      = (tMin + tMax) / 2;
guiData.winWidth       = tMax - tMin;
guiData.initialGrpVal  = initialGrpVal;
guiData.chkGrpBy       = [];
guiData.renderData     = struct();


%% ========================================================================
%  LAYOUT
%  ========================================================================

% --- Side Panel (Left, 15%) ---
panelW = 0.15;
hPanel = uipanel('Parent', hContainer, 'Units', 'normalized', ...
    'Position', [0, 0, panelW, 1]);

% --- Main Plot Panel (Right, 85%) ---
hPanelRight = uipanel('Parent', hContainer, 'Units', 'normalized', ...
    'Position', [panelW, 0, 1 - panelW, 1], 'BorderType', 'none', ...
    'BackgroundColor', 'w');

guiData.hAx = axes('Parent', hPanelRight, 'Units', 'normalized', ...
    'Position', [0.06, 0.1, 0.9, 0.85]);


% --- Side Panel Controls ---
ctlH   = 0.03;
ctlGap = 0.01;
currY  = 0.95;

% Section 1: Filter By
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

% Checkbox container (dynamic height based on categories)
guiData.pnlGrpBy = uipanel('Parent', hPanel, 'BorderType', 'none', ...
    'Units', 'normalized', 'Position', [0.05, 0.45, 0.9, currY - 0.45]);

% Section 2: Navigation
currY = 0.40;
uicontrol('Parent', hPanel, 'Style', 'text', 'String', 'Navigation:', ...
    'Units', 'normalized', 'Position', [0.05, currY, 0.9, ctlH], ...
    'HorizontalAlignment', 'left', 'FontWeight', 'bold');
currY = currY - ctlH;

% Center (s)
uicontrol('Parent', hPanel, 'Style', 'text', 'String', 'Center (s):', ...
    'Units', 'normalized', 'Position', [0.05, currY, 0.9, ctlH], ...
    'HorizontalAlignment', 'left');
currY = currY - ctlH;

guiData.edCenter = uicontrol('Parent', hPanel, 'Style', 'edit', ...
    'String', num2str(round(guiData.winCenter, 1)), ...
    'Units', 'normalized', 'Position', [0.05, currY, 0.9, ctlH], ...
    'BackgroundColor', 'w', 'Callback', @onNavChange);
currY = currY - ctlH - ctlGap;

% Window (s)
uicontrol('Parent', hPanel, 'Style', 'text', 'String', 'Window (s):', ...
    'Units', 'normalized', 'Position', [0.05, currY, 0.9, ctlH], ...
    'HorizontalAlignment', 'left');
currY = currY - ctlH;

guiData.edWindow = uicontrol('Parent', hPanel, 'Style', 'edit', ...
    'String', num2str(round(guiData.winWidth, 1)), ...
    'Units', 'normalized', 'Position', [0.05, currY, 0.9, ctlH], ...
    'BackgroundColor', 'w', 'Callback', @onNavChange);
currY = currY - ctlH - ctlGap * 2;

% Step Buttons [<] [>]
btnW = 0.43;
uicontrol('Parent', hPanel, 'Style', 'pushbutton', 'String', '<', ...
    'Units', 'normalized', 'Position', [0.05, currY, btnW, ctlH], ...
    'Callback', @(~, ~) onStep(-1));
uicontrol('Parent', hPanel, 'Style', 'pushbutton', 'String', '>', ...
    'Units', 'normalized', 'Position', [0.05 + btnW + 0.04, currY, btnW, ctlH], ...
    'Callback', @(~, ~) onStep(1));
currY = currY - ctlH - ctlGap;

% Show All Button
uicontrol('Parent', hPanel, 'Style', 'pushbutton', 'String', 'Show All', ...
    'Units', 'normalized', 'Position', [0.05, currY, 0.9, ctlH], ...
    'Callback', @onShowAll);
currY = currY - ctlH - ctlGap * 3;

% Export Button
uicontrol('Parent', hPanel, 'Style', 'pushbutton', 'String', 'Export', ...
    'Units', 'normalized', 'Position', [0.05, currY, 0.9, ctlH], ...
    'FontWeight', 'bold', 'Callback', @(~, ~) tblGUI_raster_export(hFig));


% --- Store State & Initialize ---
hContainer.UserData = guiData;
onGrpByChange(hContainer, []);


%% ========================================================================
%  CALLBACKS
%  ========================================================================

    % -----------------------------------------------------------------
    % Group-By Dropdown: populate checkboxes and replot
    % -----------------------------------------------------------------
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
                target = string(data.initialGrpVal);
                initVal = ismember(string(cats), target);
                if sum(initVal) == 0
                    warning('grpVal "%s" not found in %s. Selecting all.', ...
                        target, varName);
                    initVal = true(size(cats));
                end
                data.initialGrpVal = [];
            else
                initVal = true(size(cats));
            end

            nCats = length(cats);
            h = 1 / max(10, nCats + 1);
            for k = 1 : nCats
                yPos = 1 - k * h;
                data.chkGrpBy(k) = uicontrol('Parent', data.pnlGrpBy, ...
                    'Style', 'checkbox', 'String', cats{k}, ...
                    'Units', 'normalized', 'Position', [0, yPos, 0.9, h], ...
                    'Value', initVal(k), 'Callback', @onFilterChange);
            end
        end

        hContainer.UserData = data;
        updateActiveIndices();
        updatePlot();
    end

    % -----------------------------------------------------------------
    % Filter Checkbox: update active indices and replot
    % -----------------------------------------------------------------
    function onFilterChange(~, ~)
        updateActiveIndices();
        updatePlot();
    end

    % -----------------------------------------------------------------
    % Navigation Text Boxes: parse center and window, then replot
    % -----------------------------------------------------------------
    function onNavChange(~, ~)
        data = hContainer.UserData;

        val = str2double(get(data.edCenter, 'String'));
        if ~isnan(val)
            data.winCenter = val;
        end

        val = str2double(get(data.edWindow, 'String'));
        if ~isnan(val) && val > 0
            data.winWidth = val;
        end

        hContainer.UserData = data;
        updatePlot();
    end

    % -----------------------------------------------------------------
    % Step Buttons: shift center by half the window width
    % -----------------------------------------------------------------
    function onStep(direction)
        data = hContainer.UserData;
        data.winCenter = data.winCenter + direction * data.winWidth / 2;
        set(data.edCenter, 'String', num2str(round(data.winCenter, 1)));
        hContainer.UserData = data;
        updatePlot();
    end

    % -----------------------------------------------------------------
    % Show All: reset navigation to full filtered-data extent
    % -----------------------------------------------------------------
    function onShowAll(~, ~)
        resetNavigation();
        updatePlot();
    end

    % -----------------------------------------------------------------
    % Update Active Indices from checkbox selection
    % -----------------------------------------------------------------
    function updateActiveIndices()
        data = hContainer.UserData;
        idxVal = get(data.ddGrpBy, 'Value');
        varName = data.catVars{idxVal};

        if strcmp(varName, 'None')
            data.activeIndices = true(height(data.dataTbl), 1);
        else
            chk = data.chkGrpBy;
            if isempty(chk)
                data.activeIndices = true(height(data.dataTbl), 1);
            else
                areSel  = arrayfun(@(x) get(x, 'Value'), chk);
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

    % -----------------------------------------------------------------
    % Reset Navigation: center and window to match filtered data extent
    % -----------------------------------------------------------------
    function resetNavigation()
        data = hContainer.UserData;
        spks = data.dataTbl.(data.timesVar)(data.activeIndices);
        allT = vertcat(spks{:});

        if isempty(allT)
            data.winCenter = 0;
            data.winWidth  = 1;
        else
            tLo = min(allT);
            tHi = max(allT);
            data.winCenter = (tLo + tHi) / 2;
            data.winWidth  = tHi - tLo;
        end

        set(data.edCenter, 'String', num2str(round(data.winCenter, 1)));
        set(data.edWindow, 'String', num2str(round(data.winWidth, 1)));
        hContainer.UserData = data;
    end

    % -----------------------------------------------------------------
    % Update Plot: draw raster for active units
    % -----------------------------------------------------------------
    function updatePlot()
        data = hContainer.UserData;

        % Get filtered spike times
        spikes = data.dataTbl.(data.timesVar)(data.activeIndices);
        if data.flgBrst
            brstSpks = data.dataTbl.(data.brstVar)(data.activeIndices);
        end

        % Remove empty units (keep indexing aligned for burst overlay)
        nonEmpty = ~cellfun('isempty', spikes);
        spikes = spikes(nonEmpty);
        if data.flgBrst
            brstSpks = brstSpks(nonEmpty);
        end

        % Clear and draw
        cla(data.hAx);
        hold(data.hAx, 'on');

        % All spikes (black)
        if ~isempty(spikes)
            plot_raster(spikes, 'hAx', data.hAx, ...
                'plotType', 'vertline', 'clr', [0 0 0]);
        end

        % Burst spikes overlay (red)
        if data.flgBrst && any(~cellfun('isempty', brstSpks))
            plot_raster(brstSpks, 'hAx', data.hAx, ...
                'plotType', 'vertline', 'clr', [1 0 0]);
        end

        hold(data.hAx, 'off');

        % X-limits from navigation
        xLo = data.winCenter - data.winWidth / 2;
        xHi = data.winCenter + data.winWidth / 2;
        xlim(data.hAx, [xLo, xHi]);

        % Typography (Arial, 10 pt ticks, 12 pt labels)
        set(data.hAx, 'FontName', 'Arial', 'FontSize', 10, 'YDir', 'normal');
        xlabel(data.hAx, 'Time (s)', 'FontName', 'Arial', 'FontSize', 12);
        ylabel(data.hAx, 'Unit No.', 'FontName', 'Arial', 'FontSize', 12);

        % Cache render data for tblGUI_raster_export
        rd.spikes   = spikes;
        rd.xLo      = xLo;
        rd.xHi      = xHi;
        rd.flgBrst  = data.flgBrst;
        if data.flgBrst
            rd.brstSpks = brstSpks;
        else
            rd.brstSpks = {};
        end
        data.renderData = rd;
        hContainer.UserData = data;
    end

end     % EOF
