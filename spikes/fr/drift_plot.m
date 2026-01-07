function drift_plot(drft, ax)
% DRIFT_PLOT Plots population drift statistics.
%
%   DRIFT_PLOT(DRFT, AX) plots the correlation decay over time lags,
%   overlaid with the linear fit.
%
%   INPUTS:
%       drft        - (struct) Output from drift_calc.m
%       ax          - (axes) Optional axes handle.
%
%   See also: DRIFT_CALC, FR_NETWORK

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

if nargin < 2 || isempty(ax)
    figure('Color', 'w', 'Name', 'Drift Analysis');
    ax = gca;
end

%% ========================================================================
%  PLOT
%  ========================================================================

% Prepare Data
m_corr = drft.m_corr;
dt_corr = drft.dt_corr;
lin_coef = drft.lin_coef;
xVal = 1:length(m_corr);

hold(ax, 'on');

% 1. Plot Individual Correlations (Grey Dots)
% -------------------------------------------
% Unwrap cell array for plotting
for i = 1:length(dt_corr)
    vals = dt_corr{i};
    if ~isempty(vals)
        plot(ax, repmat(i, size(vals)), vals, '.', ...
            'Color', [0.7 0.7 0.7], 'MarkerSize', 8);
    end
end

% 2. Plot Mean Correlation (Black Dots)
% -------------------------------------
plot(ax, xVal, m_corr, 'o', ...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w', 'MarkerSize', 8);

% 3. Plot Linear Fit (Red Line)
% -----------------------------
if ~isempty(lin_coef) && ~any(isnan(lin_coef))
    yFit = polyval(lin_coef, xVal);
    plot(ax, xVal, yFit, '-', 'Color', '#D95319', 'LineWidth', 2);
end

% Aesthetics
% ----------
grid(ax, 'on');
set(ax, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 12);

xlabel(ax, 'Time Lag (Windows)');
ylabel(ax, 'PV Correlation');

titleStr = 'Population Drift';
if ~isempty(drft.drate)
    titleStr = sprintf('Drift Rate: %.4f / win', drft.drate);
end
title(ax, titleStr, 'FontWeight', 'bold');

ylim(ax, [0 1.05]);
xlim(ax, [0.5, length(m_corr) + 0.5]);


end     % EOF
