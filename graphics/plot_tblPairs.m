function hndFig = plot_tblPairs(tbl, varargin)
% PLOT_TBLPAIRS Creates a tiled layout with scatter plots for pairs of table variables.
%
% SUMMARY:
% This function creates a figure with a tiled layout where each tile contains
% a scatter plot of a pair of variables from the input table. It uses
% plot_scatterCorr for each individual scatter plot.
%
% INPUT (Required):
%   tbl           - Table containing variables to plot
%
% INPUT (Optional Key-Value Pairs):
%   pairIdx       - nx2 matrix of indices to table variables to plot {[]}
%                   If empty, plots all possible pairs
%   flgOtl        - Logical flag to remove outliers {false}
%   tileLayout    - 2-element vector [rows, cols] for tile layout {[]}
%                   If empty, automatically calculated
%   varsExc       - Cell array of variable names to exclude from analysis {[]}
%   varsInc       - Cell array of variable names to include in analysis {[]}
%                   If provided, only these variables will be used (overrides varsExc)
%
% OUTPUT:
%   hndFig        - Handle to the figure containing the plots
%
% DEPENDENCIES:
%   plot_scatterCorr.m
%
% HISTORY:
%   Aug 2024 (AI Assisted) - Initial version

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING & INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define input parameters and their defaults/validation functions
p = inputParser;
addRequired(p, 'tbl', @istable);
addParameter(p, 'pairIdx', [], @(x) isempty(x) || (isnumeric(x) && size(x, 2) == 2));
addParameter(p, 'flgOtl', false, @islogical);
addParameter(p, 'tileLayout', [], @(x) isempty(x) || (isnumeric(x) && length(x) == 2));
addParameter(p, 'varsExc', [], @(x) iscell(x) || isempty(x));
addParameter(p, 'varsInc', [], @(x) iscell(x) || isempty(x));

% Parse input arguments
parse(p, tbl, varargin{:});
tbl = p.Results.tbl;
pairIdx = p.Results.pairIdx;
flgOtl = p.Results.flgOtl;
tileLayout = p.Results.tileLayout;
varsExc = p.Results.varsExc;
varsInc = p.Results.varsInc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA VALIDATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get table variable names (excluding non-numeric columns)
varNames = tbl.Properties.VariableNames;
numericVars = varfun(@isnumeric, tbl, 'OutputFormat', 'uniform');
numericVarNames = varNames(numericVars);

% Handle variable inclusion/exclusion
if ~isempty(varsInc)
    % Include only specified variables (takes precedence over varsExc)
    includeIdx = ismember(numericVarNames, varsInc);
    numericVarNames = numericVarNames(includeIdx);
    
    % Warn if any specified variables were not found
    notFound = setdiff(varsInc, varNames);
    if ~isempty(notFound)
        warning('The following variables were not found in the table: %s', ...
            strjoin(notFound, ', '));
    end
    
    % Warn if any specified variables are not numeric
    nonNumeric = setdiff(varsInc, numericVarNames);
    if ~isempty(nonNumeric)
        warning('The following variables are not numeric and will be ignored: %s', ...
            strjoin(nonNumeric, ', '));
    end
elseif ~isempty(varsExc)
    % Exclude specified variables
    excludeIdx = ismember(numericVarNames, varsExc);
    numericVarNames = numericVarNames(~excludeIdx);
    
    % Warn if any specified variables were not found
    notFound = setdiff(varsExc, varNames);
    if ~isempty(notFound)
        warning('The following variables were not found in the table: %s', ...
            strjoin(notFound, ', '));
    end
end

nVars = length(numericVarNames);

if nVars < 2
    error('Table must contain at least 2 numeric variables for pairwise plotting.');
end

% Generate pair indices if not provided
if isempty(pairIdx)
    % Create all possible pairs
    [i, j] = find(triu(ones(nVars), 1));
    pairIdx = [i, j];
end

% Validate pair indices
if ~isempty(pairIdx)
    maxIdx = max(pairIdx(:));
    if maxIdx > nVars
        error('Pair indices exceed the number of numeric variables in the table.');
    end
    if min(pairIdx(:)) < 1
        error('Pair indices must be positive integers.');
    end
end

nPairs = size(pairIdx, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TILE LAYOUT CALCULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate optimal tile layout if not provided
if isempty(tileLayout)
    % Calculate optimal grid
    nCols = ceil(sqrt(nPairs));
    nRows = ceil(nPairs / nCols);
    tileLayout = [nRows, nCols];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE CREATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create figure with tiled layout
hndFig = figure('Name', 'Table Variable Pairs', 'Position', [100, 100, 1200, 800]);
th = tiledlayout(tileLayout(1), tileLayout(2), 'TileSpacing', 'compact', 'Padding', 'compact');

% Set figure title
title(th, 'Pairwise Variable Relationships', 'FontSize', 16, 'FontWeight', 'bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT CREATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create scatter plots for each pair
for iPair = 1:nPairs
    % Get variable indices for this pair
    idx1 = pairIdx(iPair, 1);
    idx2 = pairIdx(iPair, 2);
    
    % Get variable names
    var1Name = numericVarNames{idx1};
    var2Name = numericVarNames{idx2};
    
    % Get data for this pair
    x = tbl.(var1Name);
    y = tbl.(var2Name);
    
    % Create tile
    axh = nexttile(th);
    
    % Create scatter plot using plot_scatterCorr
    plot_scatterCorr(x, y, ...
        'xLbl', var1Name, ...
        'yLbl', var2Name, ...
        'flgOtl', flgOtl, ...
        'hndAx', axh);

end

% Adjust figure appearance
set(hndFig, 'DefaultAxesFontSize', 10);

end 