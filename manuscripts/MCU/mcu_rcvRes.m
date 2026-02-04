function hFig = mcu_rcvRes(tbl, varargin)
% MCU_RCVRES Residual Analysis of Recovery Strategy vs Contribution.
%
%   HFIG = MCU_RCVRES(TBL) performs a residual analysis to disentangle the
%   relationship between Burst and Single spike recovery, independent of
%   baseline firing rate and burstiness.
%
%   It calculates residuals from an LME model: 
%       Metric ~ (pBspk + fr) * Group + (1|Name)
%
%   It compares two distinct perspectives in a 2-row figure:
%       Relative Log-Fold Change
%       Absolute Hz Change
%
%   INPUTS:
%       tbl         - (table) Data table containing:
%                       'Group', 'Name', 'fr', 'pBspk'
%                       'frBspk', 'ss_frBspk' (Burst Rates)
%                       'frSspk', 'ss_frSspk' (Single Rates)
%                       
%   OPTIONAL (Key-Value):
%       'flgPlot'   - (logical) Whether to generate the figure. {Default: true}
%
%   OUTPUTS:
%       hFig        - (handle) Figure handle.
%

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
addParameter(p, 'flgPlot', true, @islogical);
parse(p, tbl, varargin{:});

flgPlot = p.Results.flgPlot;
hFig = [];

% Ensure required columns exist
reqVars = {'Group', 'Name', 'fr', 'pBspk', ...
           'frBspk', 'ss_frBspk', 'frSspk', 'ss_frSspk'};
if ~all(ismember(reqVars, tbl.Properties.VariableNames))
    error('[MCU_RCVRES] Missing required variables in table.');
end


%% ========================================================================
%  CALCULATE METRICS
%  ========================================================================

% Pseudocount for Log Ratios (1 spike / hour)
c = 1 / 3600;

% 1. Relative Growth -> Log Fold Change
% ------------------------------------------------
% How much did it grow relative to its own baseline?
tbl.rel_Brst = log((tbl.ss_frBspk + c) ./ (tbl.frBspk + c));
tbl.rel_Sngl = log((tbl.ss_frSspk + c) ./ (tbl.frSspk + c));

% 2. Absolute Growth (Volume) -> Delta Hz
% ------------------------------------------------
% How many actual spikes were added/subtracted?
tbl.abs_Brst = tbl.ss_frBspk - tbl.frBspk;
tbl.abs_Sngl = tbl.ss_frSspk - tbl.frSspk;


%% ========================================================================
%  RESIDUAL ANALYSIS
%  ========================================================================
%  Regress out Baseline Firing Rate (fr) and Burstiness (pBspk) to see
%  the underlying correlation between Burst and Single changes.

fprintf('[MCU_RCVRES] Calculating Residuals...\n');

% Define Base Formula
frmlBase = ' ~ (pBspk + fr) * Group + (1|Name)';

% --- A. RELATIVE RESIDUALS ---
% Burst
[lmeMdl, ~, ~, ~] = lme_analyse(tbl, ['rel_Brst' frmlBase], ...
    'dist', 'normal', 'verbose', false);
tbl.R_rel_Brst = residuals(lmeMdl, 'ResidualType', 'Pearson');

% Single
[lmeMdl, ~, ~, ~] = lme_analyse(tbl, ['rel_Sngl' frmlBase], ...
    'dist', 'normal', 'verbose', false);
tbl.R_rel_Sngl = residuals(lmeMdl, 'ResidualType', 'Pearson');


% --- B. ABSOLUTE RESIDUALS ---
% Transform: Inverse Hyperbolic Sine (asinh)
% Handles negative values while compressing heavy tails (log-like behavior).
tbl.ahs_Brst = asinh(tbl.abs_Brst);
tbl.ahs_Sngl = asinh(tbl.abs_Sngl);

% Burst
[lmeMdl, ~, ~, ~] = lme_analyse(tbl, ['ahs_Brst' frmlBase], ...
    'dist', 'normal', 'verbose', false);
tbl.R_abs_Brst = residuals(lmeMdl, 'ResidualType', 'Pearson');

% Single
[lmeMdl, ~, ~, ~] = lme_analyse(tbl, ['ahs_Sngl' frmlBase], ...
    'dist', 'normal', 'verbose', false);
tbl.R_abs_Sngl = residuals(lmeMdl, 'ResidualType', 'Pearson');


%% ========================================================================
%  PLOTTING
%  ========================================================================

if flgPlot
    hFig = figure('Color', 'w', ...
          'Name', 'Recovery');
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Shared Params
    groups = unique(tbl.Group);
    cols = lines(length(groups));
    alphaVal = 0.5;
    
    % --- SUBPLOT 1: RELATIVE ---
    nexttile;
    hold on;
    title({'Relative', '(Log Fold Change Residuals)'});
    plot_scatter_residuals(tbl, groups, cols, 'R_rel_Brst', 'R_rel_Sngl', alphaVal);
    xlabel('Burst Residuals (Log)');
    ylabel('Single Residuals (Log)');
    legend('Location', 'best');
    
    % --- SUBPLOT 2: ABSOLUTE ---
    nexttile;
    hold on;
    title({'Absolute', '(asinh \Delta Hz Residuals)'});
    plot_scatter_residuals(tbl, groups, cols, 'R_abs_Brst', 'R_abs_Sngl', alphaVal);
    xlabel('Burst Residuals (asinh)');
    ylabel('Single Residuals (asinh)');
end

end

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function plot_scatter_residuals(tbl, groups, cols, xVar, yVar, alphaVal)
    % PLOT_SCATTER_RESIDUALS Helper to plot scatter and orthogonal regression
    
    % 1. Equality Line (y=x)
    lims = max(abs([tbl.(xVar); tbl.(yVar)])) * 1.1;
    plot([-lims, lims], [-lims, lims], '--k', 'DisplayName', 'Equality (y=x)');
    
    % 2. Zero Axes
    xline(0, '-', 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
    yline(0, '-', 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
    
    axis square;
    grid on;
    xlim([-lims, lims]);
    ylim([-lims, lims]);
    
    % 3. Group Plotting
    for iGrp = 1:length(groups)
        idx = tbl.Group == groups(iGrp);
        xData = tbl.(xVar)(idx);
        yData = tbl.(yVar)(idx);
        
        % Scatter
        scatter(xData, yData, 30, cols(iGrp, :), 'filled', ...
            'MarkerFaceAlpha', alphaVal, 'MarkerEdgeColor', 'none', ...
            'HandleVisibility', 'off');
        
        % Orthogonal Regression (PCA)
        if length(xData) > 2
            % Center Data
            dataXY = [xData, yData];
            [coeff, ~, ~] = pca(dataXY);
            v1 = coeff(:, 1); % First Principal Component
            mu = mean(dataXY);
            
            % Slope & Intercept
            slope = v1(2) / v1(1);
            yInt  = mu(2) - slope * mu(1);
            
            % Plot Fit (Clipped to axes)
            fLine = @(x) slope * x + yInt;
            fplot(fLine, [-lims, lims], 'Color', cols(iGrp, :), ...
                'LineWidth', 2, 'DisplayName', char(groups(iGrp)));
            
            % Annotate Angle
            regAngle = atan2d(v1(2), v1(1));
            % Ensure angle is in intuitive range (-90 to 90 or 0 to 180 depending on preference)
            % Here we just keep it simple.
            
            % Offset text slightly based on group index to avoid overlap
            yTxt = lims * (0.8 - (iGrp-1)*0.1); 
            xTxt = -lims * 0.9;
            
            text(xTxt, yTxt, sprintf('%s: %.2f (%.0f\\circ)', ...
                char(groups(iGrp)), slope, regAngle), ...
                'Color', cols(iGrp, :), 'FontSize', 8, 'FontWeight', 'bold');
        end
    end
end
