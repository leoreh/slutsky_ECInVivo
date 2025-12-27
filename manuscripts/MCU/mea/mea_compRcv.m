function mea_compRcv(tbl, varsInc, winCalc)
% MEA_COMPRCV Compares means of variables in two time windows of a table.
%
% SUMMARY:
% This function calculates the mean for each row in tbl.varsInc during
% two specified time windows (winCalc) and uses tblGUI_scatHist to plot
% the relationship between the two windows.
%
% INPUT (Required):
%   tbl         - Input table containing vector columns.
%   varsInc     - Cell array or string of variable names to process.
%   winCalc     - 2x2 numeric matrix defining the two windows.
%                 Row 1: [Start, End] for Window 1.
%                 Row 2: [Start, End] for Window 2.
%
% OUTPUT:
%   None (opens a GUI).
%
% EXAMPLE:
%   win = [1, 600; 3001, 3600];
%   mea_compRcv(myTbl, {'FR', 'BurstRate'}, win);
%
% DEPENDENCIES:
%   tblGUI_scatHist
%
%   See also: TBLGUI_SCATHIST, TBL_TNORM

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addRequired(p, 'varsInc', @(x) iscell(x) || ischar(x) || isstring(x));
addRequired(p, 'winCalc', @(x) isnumeric(x) && all(size(x) == [2 2]));

parse(p, tbl, varsInc, winCalc);

varsInc = p.Results.varsInc;
if ischar(varsInc) || isstring(varsInc)
    varsInc = cellstr(varsInc);
end

%% ========================================================================
%  CALCULATION
%  ========================================================================

tblOut = tbl;
plotVars = {};

% Windows

for iVar = 1:length(varsInc)
    varName = varsInc{iVar};

    if ~ismember(varName, tbl.Properties.VariableNames)
        warning('Variable "%s" not found in table. Skipping.', varName);
        continue;
    end

    data = tbl.(varName);

    if ~isnumeric(data)
        warning('Variable "%s" is not numeric. Skipping.', varName);
        continue;
    end

    % Check dimensions
    nC = size(data, 2);

    for iWin = 1:2
        % --- Window iWin ---
        wS = winCalc(iWin, 1);
        wE = winCalc(iWin, 2);

        currS = max(1, wS);
        currE = min(nC, wE);

        winData = data(:, currS:currE);
        mu = mean(winData, 2, 'omitnan');
        sd = std(winData, [], 2, 'omitnan');

        name = sprintf('%s_W%d', varName, iWin);
        tblOut.(name) = sd ./ mu;

        if iWin == 1
            name1 = name;
        else
            name2 = name;
        end
    end

    plotVars = [plotVars; {name1, name2}]; %#ok<AGROW>
end

%% ========================================================================
%  PLOTTING
%  ========================================================================

if isempty(plotVars)
    error('No valid variables processed.');
end

% Set defaults to the first pair
xDefault = plotVars{1, 1};
yDefault = plotVars{1, 2};

% Launch GUI
hFig = tblGUI_scatHist(tblOut, ...
    'xVar', xDefault, ...
    'yVar', yDefault, ...
    'grpVar', 'Group');

% Access scatter axes from GUI data
if isvalid(hFig)
    data = hFig.UserData;
    if isfield(data, 'hAxScatter')
        ax = data.hAxScatter;

        % Get current limits
        xl = xlim(ax);
        yl = ylim(ax);

        % Calculate global min/max to square the plot
        globalMin = min(xl(1), yl(1));
        globalMax = max(xl(2), yl(2));

        newLim = [globalMin, globalMax];

        % Apply to both axes
        xlim(ax, newLim);
        ylim(ax, newLim);

        % Add identity line
        hold(ax, 'on');
        plot(ax, newLim, newLim, 'k--', 'DisplayName', 'Identity', 'LineWidth', 1.5);
        hold(ax, 'off');
    end
end

end
