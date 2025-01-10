function fh = mcu_lmePlot(lme_tbl, lme_mdl, varargin)

% plots linear mixed effects model results according to fixed effects.
% first fixed effect determines different lines, second fixed effect
% determines x axis values. For single fixed effect, plots as x value
%
% INPUT
%   lme_tbl     table with data used in lme analysis
%   lme_mdl     fitted linear mixed effects model
%   ptype'      string specifying plot type {'line'}
%
% OUTPUT
%   fh         handle to figure
%

% 07 Jan 24

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% handle input
p = inputParser;
addOptional(p, 'clr', []);
addOptional(p, 'ptype', 'line');

parse(p, varargin{:})
clr                 = p.Results.clr;
ptype               = p.Results.ptype;

% get variables from formula
frml = char(lme_mdl.Formula);
[varsFxd, varRsp] = get_vars(frml);

% organize data for plotting
[y_data, grp_data, x_data] = get_data(lme_tbl, varsFxd, var2, varRsp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize figure
fh = figure;
hold on

% set colormap
colormap(p.clr);
clr = colormap;
n_lines = height(grp_data);
clr_idx = round(linspace(1, size(clr, 1), n_lines));

% plot according to ptype
switch p.ptype
    case 'line'
        plot_lines(y_data, x_data, clr(clr_idx, :))
end

% add legend if multiple groups
if ~isempty(varsFxd)
    legend(grp_data.(varsFxd))
end

% labels
xlabel(var2)
ylabel(varRsp)

end

% EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [varsFxd, varRsp] = get_vars(frml)
% extract variable names from formula
% Returns:
% varNames - cell array containing all x variables
% y_name - response variable name

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
    varsFxd = unique(varsFxd);
else
    % if no interactions, process normally
    x_str = regexprep(x_str, '\s*1\s*\+?\s*', '');      % remove intercept term
    x_str = regexprep(x_str, '\([^)]*\)', '');          % remove random effects
    x_vars = strtrim(strsplit(x_str, '+'));             % split by + and trim whitespace
    varsFxd = x_vars(~cellfun(@isempty, x_vars));       % remove any empty cells
end

end


function [y_data, grp_data, x_data] = get_data(tbl, var1, var2, y_name)
% organize data for plotting by grouping according to variables

if ~isempty(var1)
    % get unique groups for first variable
    grp_data = unique(tbl(:, {var1}));
    n_grps = height(grp_data);
    
    % initialize cell arrays for data
    y_means = cell(n_grps, 1);
    y_sems = cell(n_grps, 1);
    x_vals = cell(n_grps, 1);
    
    % calculate stats for each group
    for igrp = 1 : n_grps
        % get data for current group
        grp_idx = tbl.(var1) == grp_data.(var1)(igrp);
        grp_tbl = tbl(grp_idx, :);
        
        % calculate stats per subgroup
        [g_idx, g_lbls] = findgroups(grp_tbl(:, {var2}));
        y_stats = splitapply(@(x) [mean(x) std(x)/sqrt(length(x))], grp_tbl.(y_name), g_idx);
        
        % store results
        y_means{igrp} = y_stats(:, 1);
        y_sems{igrp} = y_stats(:, 2);
        x_vals{igrp} = g_lbls;
    end
    
    % organize output
    y_data = table(y_means, y_sems, 'VariableNames', {'mean', 'sem'});
    x_data = x_vals;
    
else
    % single fixed effect
    grp_data = table();
    [g_idx, g_lbls] = findgroups(tbl(:, {var2}));
    y_stats = splitapply(@(x) [mean(x) std(x)/sqrt(length(x))], tbl.(y_name), g_idx);
    
    % organize as cell for consistency
    y_data = table({y_stats(:, 1)}, {y_stats(:, 2)}, 'VariableNames', {'mean', 'sem'});
    x_data = {g_lbls};
end
end


function plot_lines(y_data, x_data, clr)
% plot data as lines with error bars

n_grps = height(y_data);

% plot each group
for igrp = 1 : n_grps
    y_mean = y_data.mean{igrp};
    y_sem = y_data.sem{igrp};
    x_vals = 1 : length(y_mean);
    
    errorbar(x_vals, y_mean, y_sem, '-o', 'Color', clr(igrp, :),...
        'LineWidth', 2, 'MarkerFaceColor', clr(igrp, :))
end

% adjust x axis based on last group
set(gca, 'XTick', x_vals, 'XTickLabel', x_data{end}.(1))

end


