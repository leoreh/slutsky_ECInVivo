function hFig = mcu_rcvRes(tbl, varargin)
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
% Transform: SymLog
% Handles negative values while compressing heavy tails (log-like behavior).
tbl.sl_Brst = symlog(tbl.abs_Brst);
tbl.sl_Sngl = symlog(tbl.abs_Sngl);

% Burst
[lmeMdl, ~, ~, ~] = lme_analyse(tbl, ['sl_Brst' frmlBase], ...
    'dist', 'normal', 'verbose', false);
tbl.R_sl_Brst = residuals(lmeMdl, 'ResidualType', 'Pearson');

% Single
[lmeMdl, ~, ~, ~] = lme_analyse(tbl, ['sl_Sngl' frmlBase], ...
    'dist', 'normal', 'verbose', false);
tbl.R_sl_Sngl = residuals(lmeMdl, 'ResidualType', 'Pearson');


%% ========================================================================
%  PLOTTING
%  ========================================================================

if flgPlot
    hFig = figure('Color', 'w', ...
          'Name', 'Recovery');
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Shared Params
    groups = unique(tbl.Group);
    cfg = mcu_cfg;
    clr = cfg.clr.grp;
    alphaVal = 0.4;
    
    % --- SUBPLOT 1: RELATIVE ---
    nexttile;
    hold on;
    title({'Relative', '(Log Fold Change Residuals)'});
    plot_scatRes(tbl, groups, clr, 'R_rel_Brst', 'R_rel_Sngl', alphaVal);
    xlabel('Burst Residuals (Log)');
    ylabel('Single Residuals (Log)');
    
    % --- SUBPLOT 2: ABSOLUTE ---
    nexttile;
    hold on;
    title({'Absolute', '(SymLog \Delta Hz Residuals)'});
    plot_scatRes(tbl, groups, clr, 'R_sl_Brst', 'R_sl_Sngl', alphaVal);
    xlabel('Burst Residuals (SymLog)');
    ylabel('Single Residuals (SymLog)');
end

end

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function plot_scatRes(tbl, groups, cols, xVar, yVar, alphaVal)
    % PLOT_SCATTER_RESIDUALS Helper to plot scatter and orthogonal regression
    
    % 1. Equality Line (y=x)
    lims = max(abs([tbl.(xVar); tbl.(yVar)])) * 1.1;
    plot([-lims, lims], [-lims, lims], '--k', 'HandleVisibility', 'off');
    
    % 2. Zero Axes
    xline(0, '-', 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
    yline(0, '-', 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
    
    axis square;
    grid off;
    xlim([-lims, lims]);
    ylim([-lims, lims]);
    
    % 3. Group Plotting
    for iGrp = 1:length(groups)
        idx = tbl.Group == groups(iGrp);
        xData = tbl.(xVar)(idx);
        yData = tbl.(yVar)(idx);
        
        % Scatter
        scatter(xData, yData, 20, cols(iGrp, :), 'filled', ...
            'MarkerFaceAlpha', alphaVal, 'MarkerEdgeColor', 'none', ...
            'HandleVisibility', 'off');
        
        % Orthogonal Regression via plot_linReg (No embedded text, manual placement)
        Stats = plot_linReg(xData, yData, 'hAx', gca, 'type', 'ortho', ...
            'clr', cols(iGrp, :), 'flgTxt', false);
        
        if ~isnan(Stats.slope)
            % Annotate Angle with custom offset
            yTxt = lims * (0.8 - (iGrp-1)*0.1);
            xTxt = -lims * 0.9;
            
            text(xTxt, yTxt, sprintf('%s: %.2f (%.0f\\circ)', ...
                char(groups(iGrp)), Stats.slope, Stats.angle), ...
                'Color', cols(iGrp, :), 'FontSize', 8, 'FontWeight', 'bold');
        end
    end
end
