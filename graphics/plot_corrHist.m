function [hFig, hTyl, hGrid] = plot_corrHist(tbl, varargin)
% PLOT_CORRHIST Creates a matrix of scatter plots and histograms.
%
% SUMMARY:
%   This function creates a matrix of scatter plots and histograms to visualize
%   correlations between variables in a table. It offers flexibility for
%   customizing which variables are plotted on rows and columns, and for
%   grouping data points by a specified variable. It's an alternative to
%   MATLAB's built-in corrplot and scatterhistogram.
%
% SYNTAX:
%   [hAxAll, hTyl] = plot_corrHist(tbl)
%   [hAxAll, hTyl] = plot_corrHist(tbl, Name, Value, ...)
%
% INPUTS (Required):
%   tbl - A MATLAB table containing the data to plot.
%
% INPUTS (Name-Value Pairs):
%   'varsInc'    - Cell array of variable names to include in the plot.
%                  (Default: all numeric variables in tbl)
%   'grpIdx'     - Name or index of the grouping variable in tbl. This will
%                  be used to color data points.
%   'varsRow'    - Cell array of variable names to plot on the rows.
%   'varsCol'    - Cell array of variable names to plot on the columns.
%   'clrGrp'     - nGroups x 4 matrix specifying [R G B Alpha] for each group.
%   'thrOut'     - Scalar between 0 and 100 to exclude data outside this
%                  percentile range. (Default: 100, no exclusion).
%
% OUTPUTS:
%   hFig    - Handle to the figure.
%   hTyl    - Handle to the tiled layout.
%   hGrid   - Cell array of axes handles.
%
% DEPENDENCIES:
%   Statistics and Machine Learning Toolbox (for unique 'stable')
%
% HISTORY:
%   Oct 2024 - Original version created.
%
% See also: corrplot, scatterhistogram

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addRequired(p, 'tbl', @istable);

% Find numeric variables by default
defaultVars = tbl.Properties.VariableNames(vartype('numeric'));
addParameter(p, 'varsInc', defaultVars, @iscellstr);
addParameter(p, 'grpIdx', [], @(x) ischar(x) || isstring(x) || isnumeric(x));
addParameter(p, 'varsRow', {}, @iscellstr);
addParameter(p, 'varsCol', {}, @iscellstr);
addParameter(p, 'clrGrp', [], @(x) isnumeric(x) && size(x,2) == 4);
addParameter(p, 'thrOut', 100, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 100);

parse(p, tbl, varargin{:});

varsInc = p.Results.varsInc;
grpIdx = p.Results.grpIdx;
varsRow = p.Results.varsRow;
varsCol = p.Results.varsCol;
clrGrp = p.Results.clrGrp;
thrOut = p.Results.thrOut;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA PREPARATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filter to only include vars that are actually in the table
varsInc = varsInc(ismember(varsInc, tbl.Properties.VariableNames));

% Prepare group information
if ~isempty(grpIdx)
    grpData = tbl.(grpIdx);
    uGrps = unique(grpData, 'stable');
    nGrps = numel(uGrps);
    if isempty(clrGrp)
        clrGrp = [lines(nGrps), 0.5.*ones(nGrps,1)];
    end
else
    nGrps = 1;
    uGrps = 1; % Dummy group
    grpData = ones(height(tbl), 1);
    if isempty(clrGrp)
        clrGrp = [0 0 0 0.5]; 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT LAYOUT ORGANIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[plotMat, nRows, nCols] = org_tiles(varsInc, varsRow, varsCol);

hFig = plot_axSize('szOnly', false, 'axShape', 'square', 'axHeight', 150 * nRows);
hTyl = tiledlayout(nRows, nCols, 'TileSpacing', 'none', 'Padding', 'tight');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hGrid = cell(nRows, nCols);
hLgd = [];

for iRow = 1:nRows
    for iCol = 1:nCols
        plotInfo = plotMat{iRow,iCol};

        % Skip empty tiles (e.g., upper triangle in default mode)
        if isempty(plotInfo)
            if isempty(varsRow) && isempty(varsCol) && iCol > iRow
                hAx = nexttile;
                axis(hAx, 'off');
                hGrid{iRow, iCol} = hAx;
            elseif ~isempty(varsRow) || ~isempty(varsCol)
                if iRow == 1 && iCol == nCols
                    hAx = nexttile;
                    axis(hAx, 'off');
                    hGrid{iRow, iCol} = hAx;
                end
            end
            continue;
        end

        % Arrange axis
        hAx = nexttile;
        hGrid{iRow, iCol} = hAx;
        set_aesthetics(hAx);
        
        if strcmp(plotInfo.type, 'histogram')
            data_orig = tbl.(plotInfo.var);
            for iGrp = 1:nGrps
                grpMask = grpData == uGrps(iGrp);
                data = data_orig(grpMask);
                
                % Preprocess data
                [data, ~] = prep_data(data, thrOut);               
                plot_hist(hAx, data, clrGrp(iGrp,:), plotInfo.orientation);
            end
        else % scatter
            xData = tbl.(plotInfo.varX);
            yData = tbl.(plotInfo.varY);           
            hSct = gobjects(nGrps, 1);
            txtLgd = cell(nGrps, 1);

            for iGrp = 1:nGrps
                grpMask = grpData == uGrps(iGrp);
                xGrp = xData(grpMask);
                yGrp = yData(grpMask);
                
                % Preprocess x and y data separately and combine indices
                [~, idxGoodX] = prep_data(xGrp, thrOut);
                [~, idxGoodY] = prep_data(yGrp, thrOut);
                idxGood = idxGoodX & idxGoodY;
                xGrp = xGrp(idxGood);
                yGrp = yGrp(idxGood);
                
                [hSct(iGrp), txtLgd{iGrp}] = plot_sct(hAx, xGrp, yGrp, clrGrp(iGrp,:));
            end
            
            % Create legend if there are valid plots
            valid_plots = isgraphics(hSct);
            if any(valid_plots)
                if nGrps > 1
                    legend(hAx, hSct(valid_plots), txtLgd(valid_plots),...
                        'Interpreter', 'tex', 'Location', 'best');
                else
                    legend(hAx, hSct(valid_plots), txtLgd{valid_plots},...
                        'Interpreter', 'tex', 'Location', 'best');
                end
            end
            
            % Store handles for the main figure legend
            if ~isempty(hSct(valid_plots)) && isempty(hLgd)
                hLgd = hSct(valid_plots);
            end

            % Y-axis label only for the first column of scatter plots
            if iCol == 1
                ylabel(hAx, plotInfo.varY, 'Interpreter', 'none');
            else
                set(hAx, 'YTickLabel', []);
            end

            % X-axis label only for the bottom row of scatter plots
            if iRow == nRows
                xlabel(hAx, plotInfo.varX, 'Interpreter', 'none');
            else
                set(hAx, 'XTickLabel', []);
            end
        end
        hold(hAx, 'off');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AXES LINKING AND FINAL TOUCHES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note only one axis can be linked because linking must occur
% simultaneously on both axis.

% Create a logical mask to identify scatter plot axes from the plotMatrix
sctMask = cellfun(@(s) ~isempty(s) && strcmp(s.type, 'scatter'), plotMat);

% Link all X-axes for each column of scatter plots
for iCol = 1:nCols
    axes_in_col = hGrid(sctMask(:, iCol), iCol);
    axes_to_link = [axes_in_col{:}]; % Convert from cell to array
    if numel(axes_to_link) > 1
        linkaxes(axes_to_link, 'x');
    end
end

% Add a single, shared legend for the entire figure if grouping is used
if ~isempty(grpIdx)   
    txtLgd = cellstr(string(uGrps));   
    hAx = nexttile(tilenum(hTyl, 1, nCols));
    set_aesthetics(hAx);

    % Create dummy data points for legend
    hold(hAx, 'on');
    hDummy = gobjects(nGrps, 1);
    for iGrp = 1:nGrps
        hDummy(iGrp) = scatter(hAx, NaN, NaN, 10, 'filled', ...
            'MarkerFaceColor', clrGrp(iGrp, 1:3), 'MarkerFaceAlpha', clrGrp(iGrp, 4));
    end
    hold(hAx, 'off');
    
    % Create legend with dummy data
    hLgd = legend(hAx, hDummy, txtLgd, 'Interpreter', 'none');
    hLgd.Location = 'northwest';
    axis(hAx, 'off');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION: ORG_TILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [plotMatrix, nRows, nCols] = org_tiles(varsInc, varsRow, varsCol)
% ORG_TILES Organizes the layout of plots in the tiled figure.
    if isempty(varsRow) && isempty(varsCol)
        % Default case: Lower triangle of all numeric variables with histograms on the diagonal.
        nRows = length(varsInc);
        nCols = length(varsInc);
        plotMatrix = cell(nRows, nCols);

        for iRow = 1:nRows
            for iCol = 1:nCols
                if iRow == iCol
                    % Diagonal: Histogram
                    plotMatrix{iRow, iCol} = struct('type', 'histogram', 'var', varsInc{iRow}, 'orientation', 'vertical');
                elseif iCol < iRow
                    % Lower triangle: Scatter plot
                    plotMatrix{iRow, iCol} = struct('type', 'scatter', 'varX', varsInc{iCol}, 'varY', varsInc{iRow});
                end
                % Upper triangle remains empty
            end
        end
    else
        % User-defined rows and columns with histograms on top and right
        if isempty(varsRow)
            varsRow = varsInc;
        end
        if isempty(varsCol)
            varsCol = varsInc;
        end

        nScatterRows = length(varsRow);
        nScatterCols = length(varsCol);

        nRows = nScatterRows + 1; % +1 for top histogram row
        nCols = nScatterCols + 1; % +1 for right histogram col

        plotMatrix = cell(nRows, nCols);

        % Top row of histograms (for varsCol)
        for iCol = 1:nScatterCols
            plotMatrix{1, iCol} = struct('type', 'histogram', 'var', varsCol{iCol}, 'orientation', 'vertical');
        end

        % Right column of histograms (for varsRow)
        for iRow = 1:nScatterRows
             % The plot matrix row index is iRow+1 because of the top hist row
            plotMatrix{iRow + 1, nCols} = struct('type', 'histogram', 'var', varsRow{iRow}, 'orientation', 'horizontal');
        end

        % Scatter plots in the middle
        for iRow = 1:nScatterRows
            for iCol = 1:nScatterCols
                % The plot matrix indices are offset by 1
                plotMatrix{iRow + 1, iCol} = struct('type', 'scatter', 'varX', varsCol{iCol}, 'varY', varsRow{iRow});
            end
        end

        % Top-right corner remains empty
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION: PLOT_HIST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_hist(hAx, data, clr, orientation)
% PLOT_HIST Plots a histogram for a single group of data.
    if nargin < 4
        orientation = 'vertical';
    end
    
    % Plot histogram
    histogram(hAx, data, 'Normalization', 'pdf', ...
        'FaceColor', clr(1:3), 'FaceAlpha', clr(4), ...
        'Orientation', orientation, 'EdgeColor', 'none');
    
    % Add mean line
    hold(hAx, 'on');
    meanVal = mean(data);
    if strcmp(orientation, 'vertical')
        ylims = ylim(hAx);
        plot(hAx, [meanVal meanVal], ylims, '--', 'Color', clr(1:3), 'LineWidth', 2);
    else
        xlims = xlim(hAx);
        plot(hAx, xlims, [meanVal meanVal], '--', 'Color', clr(1:3), 'LineWidth', 2);
    end
    
    title(hAx, '');
    set(hAx, 'XTickLabel', [], 'YTickLabel', []);
    set(hAx, 'box', 'off', 'XColor', 'none', 'YColor', 'none');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION: PLOT_SCT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hSct, txtLgd] = plot_sct(hAx, xData, yData, clr)
% PLOT_SCT Plots a scatter plot for a single group of data.
    [rho, pval] = corr(xData, yData, 'Type', 'Spearman', 'Rows', 'complete');

    hSct = scatter(hAx, xData, yData, 10, 'filled', ...
        'MarkerFaceColor', clr(1:3), 'MarkerFaceAlpha', clr(4));
    
    % pFit = polyfit(xData, yData, 1);
    % xFit = xlim(hAx);
    % yFit = polyval(pFit, xFit);
    % plot(hAx, xFit, yFit, 'Color', clr(1:3), 'LineWidth', 1.5);

    if pval < 0.0001
        pStr = 'p < 0.0001';
    else
        pStr = sprintf('p = %.4f', pval);
    end
    txtLgd = sprintf('\\rho = %.2f, %s', rho, pStr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION: SET_PLOT_AESTHETICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function set_aesthetics(hAx)
% SET_PLOT_AESTHETICS Applies standard font and style settings to an axis.
    fntSize = 12;
    fntName = 'Arial';

    set(hAx, 'FontName', fntName, 'FontSize', fntSize);
    set(hAx, 'XTickLabelRotation', 0);
    box(hAx, 'OFF');
    hold(hAx, 'on');

    hTtl = get(hAx, 'Title');
    set(hTtl, 'FontName', fntName, 'FontSize', fntSize + 4);

    hLbl = get(hAx, 'XLabel');
    set(hLbl, 'FontName', fntName, 'FontSize', fntSize + 4);
    hLbl = get(hAx, 'YLabel');
    set(hLbl, 'FontName', fntName, 'FontSize', fntSize + 4);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION: PREP_DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data, idxGood] = prep_data(data, thrOut)
% PREP_DATA Preprocesses data by removing NaNs and outliers.

% Filter NaNs
idxGood = ~isnan(data);

% Filter by percentile if needed
if thrOut < 100
    limits = prctile(data, [100-thrOut, thrOut]);
    idxGood = data >= limits(1) & data <= limits(2) & ~isnan(data);
    data = data(idxGood);
end

end