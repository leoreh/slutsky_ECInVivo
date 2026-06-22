function hFig = tblGUI_scatHist(tbl, varargin)
% TBLGUI_SCATHIST Interactive scatter plot with marginal histograms and grouping.
%
%   tblGUI_scatHist(tbl, ...) opens a GUI to visualize the table 'tbl'.
%   Allows dynamic variable selection, grouping, and interactive point
%   selection (lasso a region, or drag a point and double-click to assign).
%
%   INPUTS:
%       tbl         - (table) The data table to visualize.
%
%   OPTIONAL KEY-VALUE PAIRS:
%       'xVar'/'yVar'/'szVar'/'grpVar' - (char) initial variables.
%       'clr'       - (m x 3) RGB color matrix for groups.
%       'alpha'     - (scalar) Marker transparency (0-1).
%       'varsExclude' - (cell) Variables to exclude from dropdowns.
%       'Parent'    - (handle) uifigure or uifigure container.
%       'SelectionCallback'/'GroupByCallback' - (func) host coordination.
%       'xScale'/'yScale'/'fitType' - (char) initial scale / fit.
%
%   Built on the shared graphics/+tblgui layer (uifigure + uigridlayout).
%
%   See also: LME_ANALYSE, TBLGUI_XY

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
addParameter(p, 'clr', [], @(x) isempty(x) || size(x, 2) == 3);
addParameter(p, 'Parent', [], @(x) isempty(x) || isgraphics(x));
addParameter(p, 'SelectionCallback', [], @(x) isempty(x) || isa(x, 'function_handle'));
addParameter(p, 'GroupByCallback', [], @(x) isempty(x) || isa(x, 'function_handle'));
addParameter(p, 'xScale', '', @(x) any(strcmpi(x, {'', 'linear', 'log'})));
addParameter(p, 'yScale', '', @(x) any(strcmpi(x, {'', 'linear', 'log', 'x-linked'})));
addParameter(p, 'fitType', '', @(x) any(strcmpi(x, {'', 'none', 'linear', 'ortho'})));
parse(p, tbl, varargin{:});

tbl         = p.Results.tbl;
varsExclude = p.Results.varsExclude;
hParent     = p.Results.Parent;
selCbk      = p.Results.SelectionCallback;
grpCbk      = p.Results.GroupByCallback;
defAlpha    = p.Results.alpha;
defClr      = p.Results.clr;
xVarIn      = p.Results.xVar;
yVarIn      = p.Results.yVar;
szVarIn     = p.Results.szVar;
grpVarIn    = p.Results.grpVar;

%% ========================================================================
%  INITIALIZATION
%  ========================================================================

% Scalar-per-row numerics for X/Y/Size (exclude id columns and matrices).
[numericVars, catVars] = tblgui.classifyVars(tbl, 'Shape', 'vector', 'Exclude', varsExclude);
catVars = setdiff(catVars, varsExclude, 'stable');

if isempty(numericVars)
    error('Input table must contain at least one numeric variable.');
end

% Defaults
curX = numericVars{1};
curY = numericVars{min(2, length(numericVars))};
curSize = 'None';
curGrp  = 'None';
if ~isempty(catVars), curGrp = catVars{1}; end

if ~isempty(xVarIn) && ismember(xVarIn, numericVars),   curX = char(xVarIn); end
if ~isempty(yVarIn) && ismember(yVarIn, numericVars),   curY = char(yVarIn); end
if ~isempty(szVarIn) && ismember(szVarIn, numericVars), curSize = char(szVarIn); end
if ~isempty(grpVarIn) && ismember(grpVarIn, catVars),   curGrp = char(grpVarIn); end

% Initial scale / fit dropdown values
xScale0 = mapChoice(p.Results.xScale, {'linear', 'Linear'; 'log', 'Log'}, 'Linear');
yScale0 = mapChoice(p.Results.yScale, {'linear', 'Linear'; 'log', 'Log'; 'x-linked', 'X-Linked'}, 'Linear');
fit0    = mapChoice(p.Results.fitType, {'none', 'None'; 'linear', 'Linear'; 'ortho', 'Ortho'}, 'None');

% Figure / parent
if isempty(hParent)
    hContainer = uifigure('Name', 'Table Visualizer', 'Position', [100, 100, 1110, 900]);
    hFig = hContainer;
else
    hContainer = hParent;
    hFig = ancestor(hContainer, 'figure');
end

guiData = struct();
guiData.tbl          = tbl;
guiData.numericVars  = numericVars;
guiData.catVars      = catVars;
guiData.defAlpha     = defAlpha;
guiData.defClr       = defClr;
guiData.selCbk       = selCbk;
guiData.grpCbk       = grpCbk;
guiData.highlightFcn   = @highlightPoints;
guiData.setGroupVarFcn = @setGroupVar;
guiData.setXYVarsFcn   = @setXYVars;
guiData.hHighlight   = [];
guiData.fitHandles   = [];
guiData.hEquality    = [];
guiData.chkGrp       = gobjects(0);
guiData.pointHandles = [];

%% ========================================================================
%  LAYOUT
%  ========================================================================

[~, gPlot, gCtrl, gActions] = tblgui.layout(hContainer, 'CtrlWidth', 230);

% Plot area: 2x2 grid -> top histogram (X), scatter, right histogram (Y).
gScat = uigridlayout(gPlot, [2, 2], 'RowHeight', {'1x', '3.5x'}, ...
    'ColumnWidth', {'4x', '1x'}, 'Padding', [2 2 2 2], 'RowSpacing', 2, 'ColumnSpacing', 2);
guiData.hAxHistX = uiaxes(gScat);   guiData.hAxHistX.Layout.Row = 1; guiData.hAxHistX.Layout.Column = 1;
guiData.hAxScatter = uiaxes(gScat); guiData.hAxScatter.Layout.Row = 2; guiData.hAxScatter.Layout.Column = 1;
guiData.hAxHistY = uiaxes(gScat);   guiData.hAxHistY.Layout.Row = 2; guiData.hAxHistY.Layout.Column = 2;

% Controls (stacked in the scrollable control column)
guiData.ddX = tblgui.labeledControl(gCtrl, 'dropdown', 'X Variable:', ...
    'Items', numericVars, 'Value', curX, 'ValueChangedFcn', @onUpdatePlot);
guiData.ddXScale = tblgui.labeledControl(gCtrl, 'dropdown', 'X Scale:', ...
    'Items', {'Linear', 'Log'}, 'Value', xScale0, 'ValueChangedFcn', @onUpdatePlot);
guiData.ddY = tblgui.labeledControl(gCtrl, 'dropdown', 'Y Variable:', ...
    'Items', numericVars, 'Value', curY, 'ValueChangedFcn', @onUpdatePlot);
guiData.ddYScale = tblgui.labeledControl(gCtrl, 'dropdown', 'Y Scale:', ...
    'Items', {'Linear', 'Log', 'X-Linked'}, 'Value', yScale0, 'ValueChangedFcn', @onUpdatePlot);
guiData.chkDisc = tblgui.labeledControl(gCtrl, 'checkbox', '', ...
    'Text', 'Discretize (binned)', 'ValueChangedFcn', @onUpdatePlot);
guiData.chkAdapt = tblgui.labeledControl(gCtrl, 'checkbox', '', ...
    'Text', 'Adaptive bins', 'ValueChangedFcn', @onUpdatePlot);
guiData.chkPrct = tblgui.labeledControl(gCtrl, 'checkbox', '', ...
    'Text', 'Bins: Mean +/- SEM', 'ValueChangedFcn', @onUpdatePlot);
guiData.ddFit = tblgui.labeledControl(gCtrl, 'dropdown', 'Fit:', ...
    'Items', {'None', 'Linear', 'Ortho'}, 'Value', fit0, 'ValueChangedFcn', @onUpdatePlot);
guiData.ddSize = tblgui.labeledControl(gCtrl, 'dropdown', 'Size:', ...
    'Items', [{'None'}, numericVars], 'Value', curSize, 'ValueChangedFcn', @onUpdatePlot);
guiData.ddGrp = tblgui.labeledControl(gCtrl, 'dropdown', 'Group:', ...
    'Items', [{'None'}, catVars], 'Value', curGrp, 'ValueChangedFcn', @onGrpChange);
tblgui.labeledControl(gCtrl, 'label', 'Filter:');
guiData.pnlGrp = tblgui.labeledControl(gCtrl, 'panel', '', 'RowHeight', '1x');

% Action buttons (pinned at the bottom; hosts may append here, e.g. utypes_gui)
guiData.btnSelect = tblgui.labeledControl(gActions, 'button', '', ...
    'Text', 'Select Group', 'ButtonPushedFcn', @onSelectRegion);
guiData.btnSelectDot = tblgui.labeledControl(gActions, 'button', '', ...
    'Text', 'Select Dot', 'ButtonPushedFcn', @onSelectDot);
guiData.btnSave = tblgui.labeledControl(gActions, 'button', '', ...
    'Text', 'Save Table to Workspace', 'ButtonPushedFcn', @onSaveTable);
guiData.gActions = gActions;

hContainer.UserData = guiData;

% Initial population & plot
onGrpChange(hContainer);
onUpdatePlot(hContainer, []);

%% ========================================================================
%  CALLBACKS
%  ========================================================================

    function onGrpChange(~, ~)
        data = hContainer.UserData;
        populateCheckboxes(data);
        data = hContainer.UserData;
        if ~isempty(data.grpCbk)
            [~, allCats] = tblgui.selectedCats(data.chkGrp);
            data.grpCbk(data.ddGrp.Value, allCats, hContainer);
        end
        onUpdatePlot(hContainer, []);
    end

    function onFilterChange(~, ~)
        data = hContainer.UserData;
        onUpdatePlot(hContainer, []);
        if ~isempty(data.grpCbk)
            activeCats = tblgui.selectedCats(data.chkGrp);
            data.grpCbk(data.ddGrp.Value, activeCats, hContainer);
        end
    end

    function populateCheckboxes(data)
        grpName = data.ddGrp.Value;
        if strcmp(grpName, 'None')
            delete(allchild(data.pnlGrp));
            data.chkGrp = gobjects(0);
        else
            cats = tblgui.catList(data.tbl.(grpName));
            data.chkGrp = tblgui.filterPanel(data.pnlGrp, cats, @onFilterChange);
        end
        hContainer.UserData = data;
    end

    function onUpdatePlot(~, ~)
        data = hContainer.UserData;
        tbl  = data.tbl;

        xName = data.ddX.Value;
        yName = data.ddY.Value;

        % --- Size logic (quadratic scaling) ---
        if ~strcmp(data.ddSize.Value, 'None')
            rawSz = tbl.(data.ddSize.Value);
            rawSz = rawSz + abs(min(rawSz));
            lims  = prctile(rawSz, [5 95]);
            minSz = lims(1); maxSz = lims(2);
            if maxSz > minSz
                clippedSz = max(min(rawSz, maxSz), minSz);
                normSz    = (clippedSz - minSz) / (maxSz - minSz);
                szData    = 20 + 80 * (normSz .^ 5);
            else
                szData = repmat(20, length(rawSz), 1);
            end
        else
            szData = repmat(20, height(tbl), 1);
        end

        % --- Group logic & colors ---
        grpName = data.ddGrp.Value;
        groups  = ones(height(tbl), 1);
        grpLabels = {'All'};
        fullCatList = {};
        idxOf = @(~) 1;
        baseColors = lines(1);

        isGrpActive = ~strcmp(grpName, 'None');
        if isGrpActive
            rawGrp = tbl.(grpName);
            if islogical(rawGrp) || ~iscategorical(rawGrp), rawGrp = categorical(rawGrp); end
            groups = rawGrp;

            fullCatList = tblgui.catList(rawGrp);            % present categories
            [baseColors, idxOf] = tblgui.groupColors(fullCatList, 'BaseColors', data.defClr);

            % Iterate over checkbox-selected categories (color stays stable)
            if ~isempty(data.chkGrp)
                active = tblgui.selectedCats(data.chkGrp);
                grpLabels = intersect(fullCatList, active, 'stable');
            else
                grpLabels = fullCatList;
            end
        end

        xData = tbl.(xName);
        yData = tbl.(yName);
        isDisc = data.chkDisc.Value;

        % --- Scales & limits ---
        scaleX = 'linear';
        if strcmp(data.ddXScale.Value, 'Log'), scaleX = 'log'; end
        isLinked = strcmp(data.ddYScale.Value, 'X-Linked');
        if isLinked
            scaleY = scaleX;
        else
            scaleY = 'linear';
            if strcmp(data.ddYScale.Value, 'Log'), scaleY = 'log'; end
        end

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
            xLim = [jointMin, jointMax]; yLim = [jointMin, jointMax];
        end

        % Global bin edges
        nBins = 30;
        if strcmp(scaleX, 'log')
            if xLim(1) <= 0, xLim(1) = min(xData(xData > 0)); if isempty(xLim(1)), xLim(1) = 0.1; end; end
            xEdges = logspace(log10(xLim(1)), log10(xLim(2)), nBins);
        else
            xEdges = linspace(xLim(1), xLim(2), nBins);
        end
        if strcmp(scaleY, 'log')
            if yLim(1) <= 0, yLim(1) = min(yData(yData > 0)); if isempty(yLim(1)), yLim(1) = 0.1; end; end
            yEdges = logspace(log10(yLim(1)), log10(yLim(2)), nBins);
        else
            yEdges = linspace(yLim(1), yLim(2), nBins);
        end

        set(data.hAxScatter, 'XScale', scaleX, 'YScale', scaleY);
        set(data.hAxHistX,   'XScale', scaleX, 'YScale', 'linear');
        set(data.hAxHistY,   'XScale', 'linear', 'YScale', scaleY);

        % --- Drawing ---
        if ~isempty(data.fitHandles), delete(data.fitHandles(isgraphics(data.fitHandles))); end
        data.fitHandles = [];
        if ~isempty(data.hEquality), delete(data.hEquality(isgraphics(data.hEquality))); end
        data.hEquality = [];

        delete(allchild(data.hAxScatter));
        delete(allchild(data.hAxHistX));
        delete(allchild(data.hAxHistY));

        hold(data.hAxScatter, 'on');
        hold(data.hAxHistX, 'on');
        hold(data.hAxHistY, 'on');

        statsData = {};

        % Pass 1: scatter & histograms per group
        for iG = 1:length(grpLabels)
            if isGrpActive
                currCat = grpLabels{iG};
                idx = (groups == currCat);
                cG = baseColors(idxOf(currCat), :);
            else
                idx = true(height(tbl), 1);
                currCat = 'All';
                cG = baseColors(1, :);
            end
            if ~any(idx), continue; end

            xG = xData(idx); yG = yData(idx); szG = szData(idx);
            nPoints = sum(idx);

            scatAlpha = data.defAlpha;
            if isDisc, scatAlpha = scatAlpha * 0.1; szG = szG * 0.3; end

            curFit = data.ddFit.Value;
            lbl = sprintf('%s (n=%d)', string(currCat), nPoints);

            [hS, hF] = plot_scat([], xG, yG, 'hAx', data.hAxScatter, ...
                'sz', szG, 'c', cG, 'alpha', scatAlpha, 'marker', 'o', ...
                'fitType', curFit, 'flgStats', true, 'dispName', lbl);

            set(hS, 'HitTest', 'off', 'PickableParts', 'none');
            if isgraphics(hF), data.fitHandles = [data.fitHandles; hF]; end

            % Binned stats (computed pass 1, drawn pass 2)
            if isDisc
                sStruct = calcBinnedStats(xG, yG, cG, data.chkAdapt.Value, ...
                    data.chkPrct.Value, xEdges, scaleX);
                if ~isempty(sStruct), statsData{end+1} = sStruct; end %#ok<AGROW>
            end

            % Marginal histograms
            if strcmp(scaleX, 'log'), normX = 'probability'; else, normX = 'pdf'; end
            plot_hist([], xG, 'hAx', data.hAxHistX, 'bins', xEdges, 'c', cG, ...
                'scale', scaleX, 'orient', 'vertical', 'norm', normX, 'flgKDE', true, 'flgStat', false);
            if strcmp(scaleY, 'log'), normY = 'probability'; else, normY = 'pdf'; end
            plot_hist([], yG, 'hAx', data.hAxHistY, 'bins', yEdges, 'c', cG, ...
                'scale', scaleY, 'orient', 'horizontal', 'norm', normY, 'flgKDE', true, 'flgStat', false);
        end

        % Pass 2: error bars on top
        if isDisc
            for k = 1:length(statsData)
                s = statsData{k};
                valid = ~isnan(s.meds);
                if ~any(valid), continue; end
                errColor = max(0, s.color * 0.7);
                hErr = errorbar(data.hAxScatter, s.centers(valid), s.meds(valid), ...
                    s.neg(valid), s.pos(valid), 'Color', errColor, 'LineWidth', 2, ...
                    'LineStyle', '-', 'Marker', 'o', 'MarkerSize', 6, ...
                    'MarkerFaceColor', errColor, 'CapSize', 0);
                hErr.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end
        end

        xlabel(data.hAxScatter, xName, 'Interpreter', 'none');
        ylabel(data.hAxScatter, yName, 'Interpreter', 'none');
        if isGrpActive
            legend(data.hAxScatter, 'Location', 'best', 'Interpreter', 'tex');
        end

        if isLinked
            data.hEquality = plot(data.hAxScatter, xLim, xLim, 'k--', 'LineWidth', 1, ...
                'HitTest', 'off', 'PickableParts', 'none', 'HandleVisibility', 'off');
        end

        applyPaddedLimits(data.hAxScatter, xLim, yLim, scaleX, scaleY);
        grid(data.hAxScatter, 'on');
        data.hAxHistX.XAxis.Visible = 'off'; data.hAxHistX.YAxis.Visible = 'off';
        data.hAxHistY.XAxis.Visible = 'off'; data.hAxHistY.YAxis.Visible = 'off';

        drawnow;
        linkaxes([data.hAxScatter, data.hAxHistX], 'x');

        hold(data.hAxScatter, 'off');
        hold(data.hAxHistX, 'off');
        hold(data.hAxHistY, 'off');

        data.btnSelect.Enable = matlab.lang.OnOffSwitchState(isGrpActive);
        data.btnSelectDot.Enable = matlab.lang.OnOffSwitchState(isGrpActive);

        hContainer.UserData = data;
    end

% --- Helper: limit calculation ---
    function lims = calcLimits(vec, scaleType)
        if strcmp(scaleType, 'log')
            v = vec(~isnan(vec) & ~isinf(vec) & vec > 0);
            if isempty(v), lims = [0.1 1]; else, lims = [min(v), max(v)]; end
        else
            v = vec(~isnan(vec) & ~isinf(vec));
            if isempty(v), lims = [0 1]; else, lims = [min(v), max(v)]; end
        end
    end

% --- Helper: apply padded limits ---
    function applyPaddedLimits(ax, xL, yL, sX, sY)
        if strcmp(sX, 'log')
            logMin = log10(xL(1)); logMax = log10(xL(2));
            span = logMax - logMin; if span == 0, span = 1; end
            xlim(ax, [10^(logMin - 0.05*span), 10^(logMax + 0.05*span)]);
        else
            span = diff(xL); if span == 0, span = 1; end
            xlim(ax, [xL(1) - 0.05*span, xL(2) + 0.05*span]);
        end
        if strcmp(sY, 'log')
            logMin = log10(yL(1)); logMax = log10(yL(2));
            span = logMax - logMin; if span == 0, span = 1; end
            ylim(ax, [10^(logMin - 0.05*span), 10^(logMax + 0.05*span)]);
        else
            span = diff(yL); if span == 0, span = 1; end
            ylim(ax, [yL(1) - 0.05*span, yL(2) + 0.05*span]);
        end
    end

% --- Helper: binned stats ---
    function sStruct = calcBinnedStats(xG, yG, cG, isAdapt, isPrct, xEdges, scaleX)
        sStruct = [];
        if isempty(xG) || isempty(yG), return; end
        valid = ~isnan(xG) & ~isnan(yG);
        xG = xG(valid); yG = yG(valid);
        if isempty(xG), return; end

        if isAdapt
            pcts  = linspace(0, 100, 16);
            edges = unique(prctile(xG, pcts));
            if length(edges) < 2, return; end
            binIdx     = discretize(xG, edges);
            uBins      = unique(binIdx(~isnan(binIdx)))';
            nBinSlots  = length(uBins);
            useDataCtr = true; preCenters = [];
        else
            edges  = xEdges;
            binIdx = discretize(xG, edges);
            if strcmp(scaleX, 'log')
                logEdges   = log10(edges);
                logCenters = (logEdges(1:end-1) + logEdges(2:end)) / 2;
                preCenters = 10 .^ logCenters;
            else
                preCenters = (edges(1:end-1) + edges(2:end)) / 2;
            end
            uBins      = 1:length(preCenters);
            nBinSlots  = length(uBins);
            useDataCtr = false;
        end

        bMeds = nan(1, nBinSlots); bLow = nan(1, nBinSlots);
        bHigh = nan(1, nBinSlots); bCtrs = nan(1, nBinSlots);
        hasData = false;

        for iB = 1:nBinSlots
            currBinIdx = uBins(iB);
            inBin = (binIdx == currBinIdx);
            if sum(inBin) >= 5
                hasData = true;
                ySub = yG(inBin);
                if isPrct
                    mu = mean(ySub, 'omitnan'); sd = std(ySub, 'omitnan');
                    sem = sd / sqrt(sum(~isnan(ySub)));
                    bLow(iB) = mu - sem; bMeds(iB) = mu; bHigh(iB) = mu + sem;
                else
                    pp = prctile(ySub, [10, 50, 90]);
                    bLow(iB) = pp(1); bMeds(iB) = pp(2); bHigh(iB) = pp(3);
                end
                if useDataCtr
                    bCtrs(iB) = median(xG(inBin), 'omitnan');
                else
                    bCtrs(iB) = preCenters(currBinIdx);
                end
            end
        end

        if hasData
            sStruct.centers = bCtrs;
            sStruct.meds    = bMeds;
            sStruct.neg     = bMeds - bLow;
            sStruct.pos     = bHigh - bMeds;
            sStruct.color   = cG;
        end
    end

%% ========================================================================
%  SELECTION
%  ========================================================================

    function onSelectRegion(~, ~)
        data = hContainer.UserData;
        roi = drawpolygon(data.hAxScatter);
        if isempty(roi.Position), delete(roi); return; end
        processSelection(hContainer, data, roi, false);
    end

    function onSelectDot(~, ~)
        data = hContainer.UserData;
        ax = data.hAxScatter;
        delete(findobj(ax, 'Type', 'images.roi.Point'));
        hPoint = drawpoint(ax, 'Label', 'Target');
        if isempty(hPoint.Position), delete(hPoint); return; end
        snapToData(hPoint, data);
        addlistener(hPoint, 'ROIClicked', @(roiSrc, evt) onDotDoubleClick(roiSrc, evt, hContainer));
        addlistener(hPoint, 'MovingROI', @(roiSrc, evt) snapToData(roiSrc, hContainer.UserData));
        title(ax, 'Double-click to assign group. Drag to browse traces.', 'Color', 'r');
    end

    function snapToData(hPoint, data)
        pos = hPoint.Position;
        xData = data.tbl.(data.ddX.Value);
        yData = data.tbl.(data.ddY.Value);
        ax = data.hAxScatter;
        xlims = ax.XLim; ylims = ax.YLim;
        xNorm = (xData - xlims(1)) / diff(xlims);
        yNorm = (yData - ylims(1)) / diff(ylims);
        pNorm = [(pos(1)-xlims(1))/diff(xlims), (pos(2)-ylims(1))/diff(ylims)];
        distSq = (xNorm - pNorm(1)).^2 + (yNorm - pNorm(2)).^2;
        [~, minIdx] = min(distSq);
        newPos = [xData(minIdx), yData(minIdx)];
        if sum((hPoint.Position - newPos).^2) > 1e-10
            hPoint.Position = newPos;
        end
        hPoint.UserData = minIdx;
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
            processSelectionIdx(hContainer, data, idx);
            delete(hPoint);
            title(data.hAxScatter, '');
        end
    end

    function processSelection(src, data, roi, isSinglePoint)
        xData = data.tbl.(data.ddX.Value);
        yData = data.tbl.(data.ddY.Value);
        inPoints = inpolygon(xData, yData, roi.Position(:,1), roi.Position(:,2));
        nSelected = sum(inPoints);
        if nSelected == 0
            tblgui.notify(hContainer, 'No points selected.', 'info');
            delete(roi); return;
        end
        if isSinglePoint && ~isempty(data.selCbk)
            data.selCbk(inPoints);
        end
        assignGroup(src, data, inPoints, nSelected);
        delete(roi);
    end

    function processSelectionIdx(src, data, idx)
        inPoints = false(height(data.tbl), 1);
        inPoints(idx) = true;
        if ~isempty(data.selCbk), data.selCbk(inPoints); end
        assignGroup(src, data, inPoints, 1);
    end

    function assignGroup(~, data, inPoints, nSelected)
        grpName = data.ddGrp.Value;
        if strcmp(grpName, 'None'), return; end
        currentGrpCol = data.tbl.(grpName);

        if iscategorical(currentGrpCol)
            cats = categories(currentGrpCol);
        elseif islogical(currentGrpCol)
            cats = {'false', 'true'};
        else
            cats = cellstr(unique(string(currentGrpCol)));
        end

        selectedCat = tblgui.chooseDialog(hContainer, ...
            sprintf('Assign %d points to:', nSelected), cats);
        if isempty(selectedCat), return; end

        if iscategorical(currentGrpCol)
            data.tbl.(grpName)(inPoints) = selectedCat;
        elseif islogical(currentGrpCol)
            data.tbl.(grpName)(inPoints) = strcmpi(selectedCat, 'true');
        else
            data.tbl.(grpName)(inPoints) = selectedCat;
        end

        hContainer.UserData = data;
        onUpdatePlot(hContainer, []);
        fprintf('Updated %d points.\n', nSelected);
    end

    function onSaveTable(~, ~)
        data = hContainer.UserData;
        assignin('base', 'fetTbl_mod', data.tbl);
        tblgui.notify(hContainer, 'Table saved to workspace as "fetTbl_mod".', 'success');
    end

%% ========================================================================
%  EXTERNAL SETTERS (host coordination)
%  ========================================================================

    function highlightPoints(indices)
        data = hContainer.UserData;
        if isnumeric(indices)
            tmp = false(height(data.tbl), 1); tmp(indices) = true; indices = tmp;
        end
        if ~isempty(data.hHighlight), delete(data.hHighlight(isgraphics(data.hHighlight))); end
        data.hHighlight = [];
        if any(indices)
            xData = data.tbl.(data.ddX.Value);
            yData = data.tbl.(data.ddY.Value);
            data.hHighlight = plot(data.hAxScatter, xData(indices), yData(indices), ...
                'o', 'MarkerSize', 12, 'LineWidth', 2, 'Color', [1 0 1], ...
                'HandleVisibility', 'off', 'HitTest', 'off', 'PickableParts', 'none');
        end
        hContainer.UserData = data;
    end

    function setGroupVar(varName, activeCats)
        data = hContainer.UserData;
        if ~ismember(varName, data.ddGrp.Items), return; end
        needsUpdate = false;
        if ~strcmp(data.ddGrp.Value, varName)
            data.ddGrp.Value = varName;
            populateCheckboxes(data);
            data = hContainer.UserData;
            needsUpdate = true;
        end
        if nargin > 1 && ~isempty(data.chkGrp)
            for i = 1:numel(data.chkGrp)
                val = ismember(data.chkGrp(i).Text, activeCats);
                if data.chkGrp(i).Value ~= val
                    data.chkGrp(i).Value = val;
                    needsUpdate = true;
                end
            end
        end
        if needsUpdate, onUpdatePlot(hContainer, []); end
    end

    function setXYVars(xName, yName)
        data = hContainer.UserData;
        needsUpdate = false;
        if ismember(xName, data.ddX.Items) && ~strcmp(data.ddX.Value, xName)
            data.ddX.Value = xName; needsUpdate = true;
        end
        if ismember(yName, data.ddY.Items) && ~strcmp(data.ddY.Value, yName)
            data.ddY.Value = yName; needsUpdate = true;
        end
        if needsUpdate, onUpdatePlot(hContainer, []); end
    end

end

function out = mapChoice(in, pairs, default)
% Map an input choice (case-insensitive) to a dropdown item; '' -> default.
out = default;
if isempty(in), return; end
for i = 1:size(pairs, 1)
    if strcmpi(in, pairs{i, 1}), out = pairs{i, 2}; return; end
end
end
