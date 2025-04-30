function mcu_lmeStats(fh, lme_mdl)
% MCU_LMESTATS Plot statistics from linear mixed effects model on figure
%
% INPUTS:
%   fh          handle to figure containing the plot
%   lme_mdl     fitted linear mixed effects model
%
% Plots p-values for all fixed effects (except intercept) on the figure.
% Uses dashed lines for interactions and solid lines for main effects.
% Shows exact p-values up to 3 decimal places regardless of significance.
%
% Example:
%   mcu_lmeStats(gcf, lme_mdl)
%
% 10 Jan 2024

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse model and extract statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get axis handle from figure
axh = findobj(fh, 'Type', 'axes');

% Extract model information
coef_tbl = lme_mdl.Coefficients;
coef_names = coef_tbl.Name;
p_values = coef_tbl.pValue;

% Skip intercept
coef_idx = ~contains(coef_names, 'Intercept');
coef_names = coef_names(coef_idx);
p_values = p_values(coef_idx);

% reverse order of estimates
coef_names = flipud(coef_names);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold(axh, 'on');

% Track maximum y value to update limits at end
max_y_used = ylim(axh);
max_y_used = max_y_used(2);

% Initialize y_offset and % of y-range for each step
yLimit = ylim(axh);
y_offset = yLimit(2) - 0.3;
y_offset_step = range(ylim(axh)) * 0.1;

% Plot each coefficient's statistics
for i = 1:length(coef_names)
    % Get coordinates for current estimate
    [x_coords, y_coords, is_interaction] = get_estimate_coords(char(coef_names{i}), lme_mdl, axh);

    % Skip if no valid coordinates found
    if isempty(x_coords) || isempty(y_coords)
        continue;
    end

    % Adjust y coordinates with offset (except for legend vertical line)
    if diff(x_coords) ~= 0  % For horizontal lines
        y_coords = repmat(y_offset + i * y_offset_step, [1, 2]);
    end

    % Update maximum y value used
    max_y_used = max(max_y_used, max(y_coords(:)));

    % Set line properties
    if is_interaction
        line_style = '--';
    else
        line_style = '-';
    end

    % Format p-value string
    p_str = sprintf('(%.3f)', p_values(i));

    % Single line
    plot(axh, x_coords, y_coords, 'k-', ...
        'LineStyle', line_style, 'LineWidth', 1.5, ...
        'HandleVisibility', 'off');

    % Position text based on line orientation
    if diff(x_coords) == 0  % Vertical line
        text_x = x_coords(1) + range(xlim(axh)) * 0.02;
        text_y = mean(y_coords);
    else  % Horizontal line
        text_x = mean(x_coords);
        text_y = y_coords(1) + range(ylim(axh)) * 0.04;
    end

    % Add p-value text
    text(text_x, text_y, p_str, ...
        'HorizontalAlignment', 'center', ...
        'HandleVisibility', 'off', ...
        'Parent', axh);

    % Update y-axis limits to accommodate all lines and text
    yLimit = ylim(axh);
    ylim(axh, [yLimit(1), max_y_used + range(ylim(axh)) * 0.1]);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x_coords, y_coords, is_interaction] = get_estimate_coords(estimate_name, lme_mdl, axh)
% GET_ESTIMATE_COORDS Returns coordinates for plotting significance lines
%
% Maps a given model estimate to coordinates on the plot where significance
% lines should be drawn, handling both main effects and interactions
%
% INPUTS:
%   estimate_name    name of the estimate from lme coefficients (e.g., 'Group_MCU-KO')
%   lme_mdl         fitted linear mixed effects model
%   axh             axis handle to the plot
%
% OUTPUTS:
%   x_coords        vector of x coordinates to connect [x1 x2]
%   y_coords        vector of y coordinates to connect [y1 y2]
%   is_interaction  boolean indicating if this is an interaction term

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse formula and get plot information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get current axis information
x_ticks = xticks(axh);
x_labels = xticklabels(axh);
y_lim = ylim(axh);

% Get variables from formula
frml = char(lme_mdl.Formula);
[varsFxd, varRsp] = get_vars(frml);

% Check if interaction model
has_interaction = contains(frml, '*');
is_interaction = contains(estimate_name, ':');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse estimate and determine coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse estimate name to get variable and level
terms = split(estimate_name, '_');
var_name = terms{1};
level_name = terms{2};

% Initialize coordinates
x_coords = [];
y_coords = [];

% Check if this is first variable estimate in additive model
is_first_var = strcmp(var_name, varsFxd{1}) && ~has_interaction && length(varsFxd) > 1;

if is_first_var
    % Place vertical line near legend
    leg = findobj(axh.Parent, 'Type', 'legend');
    if ~isempty(leg)
        leg_pos = leg.Position;
        x_base = leg_pos(1) * 0.9; % Slightly left of legend
        y_base = leg_pos(2); % Same height as legend
        x_coords = [x_base x_base];
        y_coords = [y_base y_base + 0.1];
    end
else
    % Horizontal line above relevant conditions
    y_coords = [y_lim(2) * 1.05 y_lim(2) * 1.05];  

    if is_interaction
        % Parse both parts of interaction term
        terms = split(estimate_name, ':');

        % Second term has the state/condition
        level_terms = split(terms{2}, '_');
        comp_level = level_terms{2};  % e.g., 'NREM'

        % Find bar positions at this state/condition
        comp_idx = find(strcmp(x_labels, comp_level));

        if ~isempty(comp_idx)
            % Create horizontal line spanning the two bars at this state
            bar_width = 0.15;  % Approximate bar width
            x_coords = [x_ticks(comp_idx)-bar_width x_ticks(comp_idx)+bar_width];
        end
    else

        % Get x-axis label to determine variable on x-axis
        x_var = get(axh.XLabel, 'String');

        % Find positions for reference and comparison levels
        if strcmp(var_name, x_var)
            % Variable is on x-axis
            ref_idx = 1; % First level is reference
            comp_idx = find(strcmp(x_labels, level_name));

            if ~isempty(comp_idx)
                x_coords = [x_ticks(ref_idx) x_ticks(comp_idx)];
            end
        else
            % Variable is not on x-axis - span all x values
            x_coords = [x_ticks(1) x_ticks(end)];
        end
    end
end
end


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
    varsFxd = unique(varsFxd);
else
    % if no interactions, process normally
    x_str = regexprep(x_str, '\s*1\s*\+?\s*', '');      % remove intercept term
    x_str = regexprep(x_str, '\([^)]*\)', '');          % remove random effects
    x_vars = strtrim(strsplit(x_str, '+'));             % split by + and trim whitespace
    varsFxd = x_vars(~cellfun(@isempty, x_vars));       % remove any empty cells
end

end