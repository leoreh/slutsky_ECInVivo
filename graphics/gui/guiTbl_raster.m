function hFig = guiTbl_raster(dataTbl, varargin)
% GUITBL_RASTER Interactive raster plot for Table variables.
%
%   hFig = GUITBL_RASTER(dataTbl, ...) plots spike rasters from 'dataTbl'.
%   A single raster panel shows all spikes (black) with optional burst
%   spike overlay (red). Navigation uses editable boxes for the display
%   window center and width.
%
%   INPUTS:
%       dataTbl   (table)  Table containing spike times and metadata.
%
%   OPTIONAL KEY-VALUE PAIRS:
%       'timesVar'  (char) Column name for spike times {'spktimes'}.
%       'brstVar'   (char) Column name for burst spike times {''}.
%       'grpVar'    (char) Initial grouping/filtering variable.
%       'grpVal'    (char/cell) Initial value(s) to filter by.
%       'timeLim'   (1x2 numeric) Clip spike times to [start, end] {[]}.
%       'Parent'    (handle) uifigure or uifigure container.
%
%   OUTPUT:
%       hFig      (handle) Figure handle.
%
%   Built on the shared graphics/gui layer (uifigure + uigridlayout).
%
%   See also: PLOT_RASTER, GUITBL_XY, GUITBL_RASTEREXPORT

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

flgBrst = ~isempty(brstVar) && ismember(brstVar, dataTbl.Properties.VariableNames);

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

[~, catVars] = gui_classifyVars(dataTbl);
catVars = [{'None'}, catVars];

if ~isempty(initialGrpVar) && ismember(initialGrpVar, catVars)
    grpVar = initialGrpVar;
else
    grpVar = 'None';
end

% Data extent (navigation defaults)
allSpks = vertcat(dataTbl.(timesVar){:});
if isempty(allSpks)
    tMin = 0; tMax = 1;
else
    tMin = min(allSpks); tMax = max(allSpks);
end

if isempty(hParent)
    hContainer = uifigure('Name', 'Table Raster GUI', 'Position', [100, 100, 1200, 700]);
    hFig = hContainer;
else
    hContainer = hParent;
    hFig = ancestor(hContainer, 'figure');
end

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
guiData.chkGrpBy       = gobjects(0);
guiData.renderData     = struct();

%% ========================================================================
%  LAYOUT
%  ========================================================================

[~, gPlot, gCtrl, gActions] = gui_layout(hContainer, 'CtrlWidth', 200);
guiData.hAx = uiaxes(gPlot);

% Navigation (kept compact at the top)
gui_labeledControl(gCtrl, 'label', 'Navigation:');
guiData.edCenter = gui_labeledControl(gCtrl, 'editnum', 'Center (s):', ...
    'Value', round(guiData.winCenter, 1), 'ValueChangedFcn', @onNavChange);
guiData.edWindow = gui_labeledControl(gCtrl, 'editnum', 'Window (s):', ...
    'Value', round(guiData.winWidth, 1), 'ValueChangedFcn', @onNavChange);

stepPanel = gui_labeledControl(gCtrl, 'panel', '', 'RowHeight', 32);
gStep = uigridlayout(stepPanel, [1, 2], 'Padding', [0 0 0 0], 'ColumnSpacing', 4);
uibutton(gStep, 'Text', '<', 'ButtonPushedFcn', @(~, ~) onStep(-1));
uibutton(gStep, 'Text', '>', 'ButtonPushedFcn', @(~, ~) onStep(1));

gui_labeledControl(gCtrl, 'button', '', 'Text', 'Show All', 'ButtonPushedFcn', @onShowAll);

% Filter (expands to fill remaining control-column space)
guiData.ddGrpBy = gui_labeledControl(gCtrl, 'dropdown', 'Filter By:', ...
    'Items', catVars, 'Value', grpVar, 'ValueChangedFcn', @onGrpByChange);
guiData.pnlGrpBy = gui_labeledControl(gCtrl, 'panel', '', 'RowHeight', '1x');

% Export pinned at the bottom
uibutton(gActions, 'Text', 'Export', 'FontWeight', 'bold', ...
    'ButtonPushedFcn', @(~, ~) guiTbl_rasterExport(hContainer));

hContainer.UserData = guiData;
onGrpByChange(hContainer, []);

%% ========================================================================
%  CALLBACKS
%  ========================================================================

    function onGrpByChange(~, ~)
        data = hContainer.UserData;
        varName = data.ddGrpBy.Value;
        if strcmp(varName, 'None')
            delete(allchild(data.pnlGrpBy));
            data.chkGrpBy = gobjects(0);
        else
            cats = gui_catList(data.dataTbl.(varName));
            % Initial selection from grpVal (once)
            if ~isempty(data.initialGrpVal)
                initVal = ismember(string(cats), string(data.initialGrpVal));
                if ~any(initVal)
                    warning('grpVal not found in %s. Selecting all.', varName);
                    initVal = true(size(cats));
                end
                data.initialGrpVal = [];
            else
                initVal = true(size(cats));
            end
            data.chkGrpBy = gui_filterPanel(data.pnlGrpBy, cats, @onFilterChange, ...
                'InitVal', initVal);
        end
        hContainer.UserData = data;
        updateActiveIndices();
        updatePlot();
    end

    function onFilterChange(~, ~)
        updateActiveIndices();
        updatePlot();
    end

    function onNavChange(~, ~)
        data = hContainer.UserData;
        v = data.edCenter.Value;
        if ~isnan(v), data.winCenter = v; end
        v = data.edWindow.Value;
        if ~isnan(v) && v > 0, data.winWidth = v; end
        hContainer.UserData = data;
        updatePlot();
    end

    function onStep(direction)
        data = hContainer.UserData;
        data.winCenter = data.winCenter + direction * data.winWidth / 2;
        data.edCenter.Value = round(data.winCenter, 1);
        hContainer.UserData = data;
        updatePlot();
    end

    function onShowAll(~, ~)
        resetNavigation();
        updatePlot();
    end

    function updateActiveIndices()
        data = hContainer.UserData;
        varName = data.ddGrpBy.Value;
        if strcmp(varName, 'None') || isempty(data.chkGrpBy)
            data.activeIndices = true(height(data.dataTbl), 1);
        else
            selCats = gui_selectedCats(data.chkGrpBy);
            raw = data.dataTbl.(varName);
            if islogical(raw) || ~iscategorical(raw), raw = categorical(raw); end
            data.activeIndices = ismember(string(raw), selCats);
        end
        hContainer.UserData = data;
    end

    function resetNavigation()
        data = hContainer.UserData;
        spks = data.dataTbl.(data.timesVar)(data.activeIndices);
        allT = vertcat(spks{:});
        if isempty(allT)
            data.winCenter = 0; data.winWidth = 1;
        else
            tLo = min(allT); tHi = max(allT);
            data.winCenter = (tLo + tHi) / 2;
            data.winWidth  = tHi - tLo;
        end
        data.edCenter.Value = round(data.winCenter, 1);
        data.edWindow.Value = round(data.winWidth, 1);
        hContainer.UserData = data;
    end

    function updatePlot()
        data = hContainer.UserData;

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

        cla(data.hAx);
        hold(data.hAx, 'on');

        if ~isempty(spikes)
            plot_raster(spikes, 'hAx', data.hAx, 'plotType', 'vertline', 'clr', [0 0 0]);
        end
        if data.flgBrst && any(~cellfun('isempty', brstSpks))
            plot_raster(brstSpks, 'hAx', data.hAx, 'plotType', 'vertline', 'clr', [1 0 0]);
        end

        hold(data.hAx, 'off');

        xLo = data.winCenter - data.winWidth / 2;
        xHi = data.winCenter + data.winWidth / 2;
        xlim(data.hAx, [xLo, xHi]);

        % Typography (Arial, 10 pt ticks, 12 pt labels)
        set(data.hAx, 'FontName', 'Arial', 'FontSize', 10, 'YDir', 'normal');
        xlabel(data.hAx, 'Time (s)', 'FontName', 'Arial', 'FontSize', 12);
        ylabel(data.hAx, 'Unit No.', 'FontName', 'Arial', 'FontSize', 12);

        % Cache render data for guiTbl_rasterExport
        rd.spikes  = spikes;
        rd.xLo     = xLo;
        rd.xHi     = xHi;
        rd.flgBrst = data.flgBrst;
        if data.flgBrst
            rd.brstSpks = brstSpks;
        else
            rd.brstSpks = {};
        end
        data.renderData = rd;
        hContainer.UserData = data;
    end

end     % EOF
