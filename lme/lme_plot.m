function hFig = lme_plot(lmeData, lmeMdl, varargin)

% plots linear mixed effects model results according to fixed effects.
% first fixed effect determines different lines/groups, second fixed effect
% determines x axis values. For single fixed effect, plots as x value
%
% INPUT
%   lmeData     table with data used in lme analysis
%   lmeMdl      fitted linear mixed effects model
%   ptype       string specifying plot type {'line', 'box', 'bar'}
%   hAx         axis handle
%   axShape     string specifying figure shape {'square', 'tall', 'wide'}
%   idxRow      row index for significance lines
%   lmeStats    statistics for significance lines
%
% OUTPUT
%   hFig         handle to figure

% 10 Jan 24

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% handle input
p = inputParser;
addOptional(p, 'clr', []);
addOptional(p, 'hAx', []);
addOptional(p, 'ptype', 'line');
addOptional(p, 'axShape', 'square');
addOptional(p, 'idxRow', []);
addOptional(p, 'grpVar', 'Group');
addOptional(p, 'lmeStats', []);

parse(p, varargin{:})
clr                 = p.Results.clr;
hAx                 = p.Results.hAx;
ptype               = p.Results.ptype;
axShape             = p.Results.axShape;
idxRow              = p.Results.idxRow;
grpVar              = p.Results.grpVar;
lmeStats            = p.Results.lmeStats;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% organize colors
if isempty(clr)
    clr = mcu_clr();
    clr = clr.grp;
end
clrAlpha = 0.3;

% initialize figure
if isempty(hAx)
    [hFig, hAx] = plot_axSize('axShape', axShape, 'szOnly', false,...
        'flgPos', true);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% organize data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get variables from formula
frml = char(lmeMdl.Formula);
[varsFxd, varRsp] = lme_frml2vars(frml);

% organize data for plotting based on number of variables
[dataGrp, xVals, varLbls] = lme_tbl2cell(lmeData, varsFxd, varRsp, 'flgCats', true, 'grpVar', grpVar);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine x-axis variable name
if length(varsFxd) >= 2
    if ~isempty(grpVar) && ismember(grpVar, varsFxd)
        xVarName = setdiff(varsFxd, grpVar, 'stable');
        if length(xVarName) > 1
            xVarName = xVarName{1}; % Use first non-group variable
        end
    else
        xVarName = varsFxd{2}; % Default: second variable
    end
else
    xVarName = varsFxd{1};
end

switch ptype
    case 'line'
        % plot data based on number of fixed effects
        if length(varsFxd) >= 2
            % multiple lines case (grouped by first variable)
            for igrp = 1 : length(dataGrp)
                ph(igrp) = plot_stdShade('dataMat', dataGrp{igrp},...
                    'xVal', [1:length(xVals{igrp})], ...
                    'hAx', hAx, 'clr', clr(igrp, :), 'alpha', clrAlpha);
            end
            xlabel(xVarName)
            xticks([1:length(xVals{igrp})])
            xlim([1-0.2, length(xVals{igrp})+0.2])
            xticklabels(xVals{igrp})
            legend(ph, varLbls, 'Location', 'northwest')
        else
            % single variable case
            ph = plot_stdShade('dataMat', dataGrp, 'xVal', [1:length(xVals)], ...
                'hAx', hAx, 'clr', clr(1,:), 'alpha', clrAlpha);
            xlabel(xVarName)
            xticklabels(xVals)
            xlim([1-0.2, length(xVals)+0.2])
        end
        

    case {'box', 'bar', 'allPnts'}
        if length(varsFxd) >= 2
            % pass cell array of matrices to plot_boxMean
            ph = plot_boxMean('dataMat', dataGrp, 'xVal', 1 : length(xVals{1}), ...
                'clr', clr(1 : length(dataGrp), :), 'alphaIdx', 1, ...
                'plotType', ptype, 'hAx', hAx);
            xlabel(xVarName)
            xticklabels(xVals{1})
            legend(ph, varLbls, 'Location', 'northwest')
        else
            % single variable case
            plot_boxMean('dataMat', dataGrp, 'xVal', 1 : size(dataGrp, 1),...
                'clr', clr(1, :), 'alphaIdx', clrAlpha, ...
                'plotType', ptype, 'hAx', hAx);
            xlabel(xVarName)
            xticklabels(xVals)
        end

        % update graphics to the case of only two groups (assumes Control
        % and MCU-KO)
        hBar = findobj(hAx, 'Type', 'Bar');
        nBars = numel(hBar.YData);
        if nBars == 2 && length(hBar) == 1
            hBar.FaceColor = 'flat';
            hBar.CData(1, :) = clr(1, :);
            hBar.CData(2, :) = clr(2, :);
            hBar.FaceAlpha = 1;
        end
end

% add response variable label and title
ylabel(varRsp)
title(hAx, lme_frml2char(frml, 'rmRnd', false))

% add significance lines if idxRow is provided
if ~isempty(idxRow)

    % Generate significance lines
    [barIdx, barLbl] = lme_sigLines(lmeStats, lmeMdl, lmeData, 'idxRow', idxRow, 'grpVar', grpVar);
    
    % Add significance lines to the plot
    if ~isempty(barIdx)
        plot_sigLines(hAx, barIdx, barLbl, 'flgNS', false);
    end
end

% adjust figure size and aesthetics
drawnow;
plot_axSize('hFig', hFig, 'axShape', axShape, 'szOnly', true);

end

