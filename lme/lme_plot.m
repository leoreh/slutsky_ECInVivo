function fh = lme_plot(lmeData, lmeMdl, varargin)

% plots linear mixed effects model results according to fixed effects.
% first fixed effect determines different lines/groups, second fixed effect
% determines x axis values. For single fixed effect, plots as x value
%
% INPUT
%   lmeData     table with data used in lme analysis
%   lmeMdl      fitted linear mixed effects model
%   ptype       string specifying plot type {'line', 'box', 'bar'}
%   axh         axis handle
%   figShape    string specifying figure shape {'square', 'tall', 'wide'}
%
% OUTPUT
%   fh         handle to figure
%

% 10 Jan 24

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% handle input
p = inputParser;
addOptional(p, 'clr', []);
addOptional(p, 'axh', []);
addOptional(p, 'ptype', 'line');
addOptional(p, 'figShape', 'square');

parse(p, varargin{:})
clr                 = p.Results.clr;
axh                 = p.Results.axh;
ptype               = p.Results.ptype;
figShape            = p.Results.figShape;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% organize colors
if isempty(clr)
    clr = lines;
    clr(1, :) = [0.3 0.3 0.3];
    clr(2, :) = [0.784 0.667 0.392];
end
clrAlpha = 0.3;

% initialize figure
if isempty(axh)
    fh = figure;
    set(fh, 'Color', 'w');                      
    fhUnits = get(fh, 'Units'); 
    set(fh, 'Units', 'pixels');
    fhPos = get(fh, 'Position');
    switch figShape
        case 'square'
            fhPos(3) = fhPos(4);
        case 'tall'
            fhPos(3) = fhPos(4) * 0.62;
        case 'wide'
            fhPos(3) = fhPos(4) * 1.62;
    end
    set(fh, 'Position', fhPos);
    set(fh, 'Units', fhUnits);
    tlayout = [1, 1];
    th = tiledlayout(tlayout(1), tlayout(2));
    th.TileSpacing = 'none';
    th.Padding = 'tight';
    axh = nexttile(th, 1, [1, 1]); cla; hold on
    % axis(axh, 'square');        % Set axes proportions to a perfect square
    fntSize = 16;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get variables from formula
frml = char(lmeMdl.Formula);
[varsFxd, varRsp] = get_vars(frml);

% organize data for plotting based on number of variables
[dataGrp, xVals, varLbls] = get_data(lmeData, varsFxd, varRsp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ptype
    case 'line'
        % plot data based on number of fixed effects
        if length(varsFxd) >= 2
            % multiple lines case (grouped by first variable)
            for igrp = 1 : length(dataGrp)
                ph(igrp) = plot_stdShade('dataMat', dataGrp{igrp},...
                    'xVal', [1:length(xVals{igrp})], ...
                    'axh', axh, 'clr', clr(igrp, :), 'alpha', clrAlpha);
            end
            xlabel(varsFxd{2})
            xticks([1:length(xVals{igrp})])
            xlim([1-0.2, length(xVals{igrp})+0.2])
            xticklabels(xVals{igrp})
            legend(ph, varLbls, 'Location', 'northwest')
        else
            % single variable case
            ph = plot_stdShade('dataMat', dataGrp, 'xVal', [1:length(xVals)], ...
                'axh', axh, 'clr', clr(1,:), 'alpha', clrAlpha);
            xlabel(varsFxd{1})
            xticklabels(xVals)
            xlim([1-0.2, length(xVals)+0.2])
        end
        

    case {'box', 'bar', 'allPnts'}
        if length(varsFxd) >= 2
            % pass cell array of matrices to plot_boxMean
            ph = plot_boxMean('dataMat', dataGrp, 'xVal', 1 : length(xVals{1}), ...
                'clr', clr(1 : length(dataGrp), :), 'alphaIdx', 1, ...
                'plotType', ptype, 'axh', axh);
            xlabel(varsFxd{2})
            xticklabels(xVals{1})
            legend(ph, varLbls, 'Location', 'northwest')
        else
            % single variable case
            plot_boxMean('dataMat', dataGrp, 'xVal', 1 : size(dataGrp, 1),...
                'clr', clr(1, :), 'alphaIdx', clrAlpha, ...
                'plotType', ptype, 'axh', axh);
            xlabel(varsFxd{1})
            xticklabels(xVals)
        end

        % update graphics to the case of only two groups (assumes Control
        % and MCU-KO)
        hndlBar = findobj(axh, 'Type', 'Bar');
        nBars = numel(hndlBar.YData);
        if nBars == 2 && length(hndlBar) == 1
            hndlBar.FaceColor = 'flat';
            hndlBar.CData(1, :) = clr(1, :);
            hndlBar.CData(2, :) = clr(2, :);
            hndlBar.FaceAlpha = 1;
        end
end

% add response variable label and title
ylabel(varRsp)
title(axh, lme_frml2char(frml, 'rmRnd', false))

% Set fonts and font sizes for axes elements
set(axh, 'FontName', 'Arial', 'FontSize', fntSize); % Affects axis labels, ticks

% Set font and font size for the title
hndlTtl = get(axh, 'Title');
set(hndlTtl, 'FontName', 'Arial', 'FontSize', fntSize + 4);

% Set font and font size for the legend
hndlLgnd = findobj(fh, 'Type', 'Legend');
set(hndlLgnd, 'FontName', 'Arial', 'FontSize', fntSize);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [varsFxd, varRsp] = get_vars(frml)
% extract variable names from formula

% get y variable (response)
varRsp = regexp(frml, '(\w+)\s*~', 'tokens');
varRsp = varRsp{1}{1};

% get x variables (fixed effects) excluding random effects and intercept
xStr = regexp(frml, '~\s*(.*?)\s*(\(|$)', 'tokens'); % get everything between ~ and (, or end
xStr = xStr{1}{1}; % extract matched string

% Initialize varNames as empty cell array
varsFxd = {};

% first check for interaction terms
interactVars = regexp(xStr, '(\w+)\s*[*]\s*(\w+)', 'tokens');
if ~isempty(interactVars)
    % Process all interaction pairs
    for i = 1:length(interactVars)
        varsFxd = [varsFxd, interactVars{i}];
    end
    % Remove duplicates
    varsFxd = unique(varsFxd, 'stable');
else
    % if no interactions, process normally
    xStr = regexprep(xStr, '\s*1\s*\+?\s*', '');      % remove intercept term
    xStr = regexprep(xStr, '\([^)]*\)', '');          % remove random effects
    xVars = strtrim(strsplit(xStr, '+'));             % split by + and trim whitespace
    varsFxd = xVars(~cellfun(@isempty, xVars));       % remove any empty cells
end

end


function [dataGrp, xVals, varLbls] = get_data(tbl, varsFxd, varRsp)
% organize data for plotting based on number of fixed effects
% Returns raw data matrices grouped by variables

if length(varsFxd) >= 2
    % case: two variables - group by first var, x-axis is second var
    var1 = varsFxd{1};
    var2 = varsFxd{2};

    % get unique values maintaining order
    var1Vals = unique(tbl.(var1), 'stable');

    % initialize outputs
    nGrps = length(var1Vals);
    dataGrp = cell(nGrps, 1);
    xVals = cell(nGrps, 1);
    varLbls = strings(nGrps, 1);

    % get data for each group
    for igrp = 1:nGrps
        % get data for current group
        idx = tbl.(var1) == var1Vals(igrp);
        grpTbl = tbl(idx, :);

        % get unique x values using categorical order if possible
        if iscategorical(grpTbl.(var2))
            xUnique = categories(grpTbl.(var2));
        else
            xUnique = unique(grpTbl.(var2), 'stable');
        end
        nX = length(xUnique);

        % organize data matrix
        dataMat = nan(nX, sum(idx));
        for ix = 1:nX
            xIdx = grpTbl.(var2) == xUnique{ix};
            yVals = grpTbl.(varRsp)(xIdx);
            dataMat(ix, 1:length(yVals)) = yVals';
        end

        % store results
        dataGrp{igrp} = dataMat;
        xVals{igrp} = string(xUnique);
        varLbls{igrp} = char(var1Vals(igrp));
    end

else
    % case: single variable
    var1 = varsFxd{1};

    % get unique values using categorical order if possible
    if iscategorical(tbl.(var1))
        xUnique = categories(tbl.(var1));
    else
        xUnique = unique(tbl.(var1), 'stable');
    end
    nX = length(xUnique);

    % initialize data matrix
    maxPts = max(histcounts(categorical(tbl.(var1))));
    dataMat = nan(nX, maxPts);

    % fill data matrix
    for ix = 1:nX
        idx = tbl.(var1) == xUnique{ix};
        yVals = tbl.(varRsp)(idx);
        dataMat(ix, 1:length(yVals)) = yVals';
    end

    % prepare outputs
    dataGrp = dataMat;
    xVals = char(xUnique);
    varLbls = [];
end

end