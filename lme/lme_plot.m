function fh = lme_plot(lme_tbl, lme_mdl, varargin)

% plots linear mixed effects model results according to fixed effects.
% first fixed effect determines different lines/groups, second fixed effect
% determines x axis values. For single fixed effect, plots as x value
%
% INPUT
%   lme_tbl     table with data used in lme analysis
%   lme_mdl     fitted linear mixed effects model
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
clr_alpha = 0.3;

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
frml = char(lme_mdl.Formula);
[varsFxd, varRsp] = get_vars(frml);

% organize data for plotting based on number of variables
[data_grp, x_vals, var_lbls] = get_data(lme_tbl, varsFxd, varRsp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ptype
    case 'line'
        % plot data based on number of fixed effects
        if length(varsFxd) >= 2
            % multiple lines case (grouped by first variable)
            for igrp = 1 : length(data_grp)
                ph(igrp) = plot_stdShade('dataMat', data_grp{igrp},...
                    'xVal', [1:length(x_vals{igrp})], ...
                    'axh', axh, 'clr', clr(igrp, :), 'alpha', clr_alpha);
            end
            xlabel(varsFxd{2})
            xticks([1:length(x_vals{igrp})])
            xlim([1-0.2, length(x_vals{igrp})+0.2])
            xticklabels(x_vals{igrp})
            legend(ph, var_lbls, 'Location', 'northwest')
        else
            % single variable case
            ph = plot_stdShade('dataMat', data_grp, 'xVal', [1:length(x_vals)], ...
                'axh', axh, 'clr', clr(1,:), 'alpha', clr_alpha);
            xlabel(varsFxd{1})
            xticklabels(x_vals)
            xlim([1-0.2, length(x_vals)+0.2])
        end
        

    case {'box', 'bar', 'allPnts'}
        if length(varsFxd) >= 2
            % pass cell array of matrices to plot_boxMean
            ph = plot_boxMean('dataMat', data_grp, 'xVal', 1 : length(x_vals{1}), ...
                'clr', clr(1 : length(data_grp), :), 'alphaIdx', 1, ...
                'plotType', ptype, 'axh', axh);
            xlabel(varsFxd{2})
            xticklabels(x_vals{1})
            legend(ph, var_lbls, 'Location', 'northwest')
        else
            % single variable case
            plot_boxMean('dataMat', data_grp, 'xVal', 1 : size(data_grp, 1),...
                'clr', clr(1, :), 'alphaIdx', clr_alpha, ...
                'plotType', ptype, 'axh', axh);
            xlabel(varsFxd{1})
            xticklabels(x_vals)
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
title(axh, lme_frml2char(frml, 'rm_rnd', false))

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
x_str = regexp(frml, '~\s*(.*?)\s*(\(|$)', 'tokens'); % get everything between ~ and (, or end
x_str = x_str{1}{1}; % extract matched string

% Initialize varNames as empty cell array
varsFxd = {};

% first check for interaction terms
interact_vars = regexp(x_str, '(\w+)\s*[*]\s*(\w+)', 'tokens');
if ~isempty(interact_vars)
    % Process all interaction pairs
    for i = 1:length(interact_vars)
        varsFxd = [varsFxd, interact_vars{i}];
    end
    % Remove duplicates
    varsFxd = unique(varsFxd, 'stable');
else
    % if no interactions, process normally
    x_str = regexprep(x_str, '\s*1\s*\+?\s*', '');      % remove intercept term
    x_str = regexprep(x_str, '\([^)]*\)', '');          % remove random effects
    x_vars = strtrim(strsplit(x_str, '+'));             % split by + and trim whitespace
    varsFxd = x_vars(~cellfun(@isempty, x_vars));       % remove any empty cells
end

end


function [data_grp, x_vals, var_lbls] = get_data(tbl, varsFxd, varRsp)
% organize data for plotting based on number of fixed effects
% Returns raw data matrices grouped by variables

if length(varsFxd) >= 2
    % case: two variables - group by first var, x-axis is second var
    var1 = varsFxd{1};
    var2 = varsFxd{2};

    % get unique values maintaining order
    var1_vals = unique(tbl.(var1), 'stable');

    % initialize outputs
    n_grps = length(var1_vals);
    data_grp = cell(n_grps, 1);
    x_vals = cell(n_grps, 1);
    var_lbls = strings(n_grps, 1);

    % get data for each group
    for igrp = 1:n_grps
        % get data for current group
        idx = tbl.(var1) == var1_vals(igrp);
        grp_tbl = tbl(idx, :);

        % get unique x values using categorical order if possible
        if iscategorical(grp_tbl.(var2))
            x_unique = categories(grp_tbl.(var2));
        else
            x_unique = unique(grp_tbl.(var2), 'stable');
        end
        n_x = length(x_unique);

        % organize data matrix
        data_mat = nan(n_x, sum(idx));
        for ix = 1:n_x
            x_idx = grp_tbl.(var2) == x_unique{ix};
            y_vals = grp_tbl.(varRsp)(x_idx);
            data_mat(ix, 1:length(y_vals)) = y_vals';
        end

        % store results
        data_grp{igrp} = data_mat;
        x_vals{igrp} = string(x_unique);
        var_lbls{igrp} = char(var1_vals(igrp));
    end

else
    % case: single variable
    var1 = varsFxd{1};

    % get unique values using categorical order if possible
    if iscategorical(tbl.(var1))
        x_unique = categories(tbl.(var1));
    else
        x_unique = unique(tbl.(var1), 'stable');
    end
    n_x = length(x_unique);

    % initialize data matrix
    max_pts = max(histcounts(categorical(tbl.(var1))));
    data_mat = nan(n_x, max_pts);

    % fill data matrix
    for ix = 1:n_x
        idx = tbl.(var1) == x_unique{ix};
        y_vals = tbl.(varRsp)(idx);
        data_mat(ix, 1:length(y_vals)) = y_vals';
    end

    % prepare outputs
    data_grp = data_mat;
    x_vals = char(x_unique);
    var_lbls = [];
end

end