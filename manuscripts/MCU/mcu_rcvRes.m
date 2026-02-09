function hFig = mcu_rcvRes(tbl, varargin)
%   HFIG = MCU_RCVRES(TBL) performs a residual analysis to disentangle the
%   relationship between Burst and Single spike recovery, independent of
%   baseline firing rate and burstiness.
%
%   It compares two distinct perspectives in a 2x3 figure:
%       Rows:    Relative (Log-Fold), Absolute (SymLog Hz)
%       Columns: Allocation, Burst Prediction, Single Prediction
%
%   INPUTS:
%       tbl         - (table) Data table containing:
%                       'Group', 'Name', 'fr', 'pBspk'
%                       'frBspk', 'ss_frBspk', 'frSspk', 'ss_frSspk'
%                       'pBspk_trans' (Logit transformed burstiness)
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
%  ALLOCATION ANALYSIS 
%  ========================================================================
%  Regress out Baseline Firing Rate (fr) and Burstiness (pBspk) to see
%  the underlying correlation between Burst and Single changes.

fprintf('[MCU_RCVRES] Calculating Allocation Residuals...\n');

frmlBase = ' ~ (pBspk + fr) * Group + (1|Name)';

% Relative
mdl = lme_analyse(tbl, ['dBrst_rel' frmlBase], 'dist', 'normal', 'verbose', false);
tbl.R_dBrst_rel = residuals(mdl, 'ResidualType', 'Pearson');

mdl = lme_analyse(tbl, ['dSngl_rel' frmlBase], 'dist', 'normal', 'verbose', false);
tbl.R_dSngl_rel = residuals(mdl, 'ResidualType', 'Pearson');

% Absolute
mdl = lme_analyse(tbl, ['dBrst_abs' frmlBase], 'dist', 'normal', 'verbose', false);
tbl.R_dBrst_abs = residuals(mdl, 'ResidualType', 'Pearson');

mdl = lme_analyse(tbl, ['dSngl_abs' frmlBase], 'dist', 'normal', 'verbose', false);
tbl.R_dSngl_abs = residuals(mdl, 'ResidualType', 'Pearson');


%% ========================================================================
%  PREDICTION ANALYSIS
%  ========================================================================
%  Does Baseline Burstiness (pBspk) predict the deviation from the 
%  standard plasticity rule? Regress out the "Other" component.

fprintf('[MCU_RCVRES] Calculating Prediction Residuals...\n');

% --- Relative ---
% Burst (Control for Single)
frml_Brst = ' ~ (fr + dSngl_rel) * Group + (1|Name)';
mdl = lme_analyse(tbl, ['dBrst_rel' frml_Brst], 'dist', 'normal', 'verbose', false);
tbl.P_dBrst_rel = residuals(mdl, 'ResidualType', 'Pearson');

% Single (Control for Burst)
frml_Sngl = ' ~ (fr + dBrst_rel) * Group + (1|Name)';
mdl = lme_analyse(tbl, ['dSngl_rel' frml_Sngl], 'dist', 'normal', 'verbose', false);
tbl.P_dSngl_rel = residuals(mdl, 'ResidualType', 'Pearson');

% --- Absolute ---
% Burst (Control for Single)
frml_Brst = ' ~ (fr + dSngl_abs) * Group + (1|Name)';
mdl = lme_analyse(tbl, ['dBrst_abs' frml_Brst], 'dist', 'normal', 'verbose', false);
tbl.P_dBrst_abs = residuals(mdl, 'ResidualType', 'Pearson');

% Single (Control for Burst)
frml_Sngl = ' ~ (fr + dBrst_abs) * Group + (1|Name)';
mdl = lme_analyse(tbl, ['dSngl_abs' frml_Sngl], 'dist', 'normal', 'verbose', false);
tbl.P_dSngl_abs = residuals(mdl, 'ResidualType', 'Pearson');


%% ========================================================================
%  PLOTTING
%  ========================================================================

if ~flgPlot
    return;
end

hFig = figure('Color', 'w', 'Name', 'Recovery: Allocation & Prediction');
hTile = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Definitions for Loop (Values correspond to suffixes of calculated table vars)
rowDefs = {
    'Relative', 'dBrst_rel',    'dSngl_rel',    'Single Gain',     'Burst Gain';
    'Absolute', 'dBrst_abs',    'dSngl_abs',    '\Delta Single',   '\Delta Burst'
};

for iRow = 1:2
    typeStr = rowDefs{iRow, 1};
    sufB = rowDefs{iRow, 2}; % Suffix Burst
    sufS = rowDefs{iRow, 3}; % Suffix Single
    lblS = rowDefs{iRow, 4}; % Label Single
    lblB = rowDefs{iRow, 5}; % Label Burst
    
    % --- 1. Allocation (Burst vs Single) ---
    nexttile; hold on;
    title([typeStr ' Allocation']);
    
    xVar = ['R_' sufB];
    yVar = ['R_' sufS];
    
    % Reference Lines
    lims = max(abs([tbl.(xVar); tbl.(yVar)])) * 1.1;
    
    if strcmpi(typeStr, 'Absolute')
        lims = symlog(lims);
    end
    
    isSymlog = strcmpi(typeStr, 'Absolute');
    plot_grpScatter(tbl, xVar, yVar, lims, isSymlog);

    xlim([-lims, lims]); ylim([-lims, lims]);
    xlabel(lblB, 'Interpreter', 'tex'); 
    ylabel(lblS, 'Interpreter', 'tex');
     
    plot([-lims, lims], [-lims, lims], '--k', 'HandleVisibility', 'off'); % Identity
    xline(0, '-', 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
    yline(0, '-', 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
    
    % --- 2. Burst Prediction ---
    nexttile; hold on;
    title([typeStr ' Burst Recovery']);
    
    yVar = ['P_' sufB];
    xVar = 'pBspk_trans';
    
    yline(0, '-', 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
    
    isSymlog = strcmpi(typeStr, 'Absolute');
    plot_grpScatter(tbl, xVar, yVar, [], isSymlog);
    
    xlabel('Baseline Burstiness'); 
    ylabel(lblB, 'Interpreter', 'tex');

    % --- 3. Single Prediction ---
    nexttile; hold on;
    title([typeStr ' Single Recovery']);
    
    yVar = ['P_' sufS];
    xVar = 'pBspk_trans';
    
    yline(0, '-', 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
    
    isSymlog = strcmpi(typeStr, 'Absolute');
    plot_grpScatter(tbl, xVar, yVar, [], isSymlog);
    
    xlabel('Baseline Burstiness'); 
    ylabel(lblS, 'Interpreter', 'tex');
end

end

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function plot_grpScatter(tbl, xVar, yVar, lims, isSymlog)
    % PLOT_GRPSCATTER Helper for common grouping scatter/regression logic
    
    if nargin < 5; isSymlog = false; end
    
    regType = 'ortho';
    
    cfg = mcu_cfg;
    clr = cfg.clr.grp;
    groups = unique(tbl.Group);
    alphaVal = 0.4;
    
    axis square;
    
    % 1. Plot Scatter (Raw)
    for iGrp = 1:length(groups)
        idx = tbl.Group == groups(iGrp);
        x = tbl.(xVar)(idx); 
        y = tbl.(yVar)(idx);
        c = clr(iGrp, :);
        
        scatter(x, y, 20, c, 'filled', 'MarkerFaceAlpha', alphaVal, ...
            'MarkerEdgeColor', 'none', 'HandleVisibility', 'off');
    end

    % 2. Transform Axes
    if isSymlog
        if contains(xVar, 'pBspk')
             symlog(gca, 'y');
        else
             symlog(gca, 'xy');
        end
    end

    % 3. Plot Regression (on potentially transformed data)
    for iGrp = 1:length(groups)
        idx = tbl.Group == groups(iGrp);
        x = tbl.(xVar)(idx); 
        y = tbl.(yVar)(idx);
        c = clr(iGrp, :);

        if isSymlog
            if contains(xVar, 'pBspk')
                 % Only Y is symlog
                 yReg = symlog(y);
                 xReg = x;
            else
                 % Both
                 xReg = symlog(x);
                 yReg = symlog(y);
            end
        else
            xReg = x;
            yReg = y;
        end
        
        % Regression
        Stats = plot_linReg(xReg, yReg, 'hAx', gca, 'type', regType, ...
            'clr', c, 'flgTxt', false);
        
        % Annotation
        if ~isnan(Stats.slope)
            % Orthogonal: Show Angle/Slope using Normalized Units (Robust to limits/symlog)
            % Position: Top-Left
            txtX = 0.02; 
            txtY = 0.98 - (iGrp-1) * 0.06;
            
            text(txtX, txtY, sprintf('%s: %.2f (%.0f\\circ)', ...
                char(groups(iGrp)), Stats.slope, Stats.angle), ...
                'Units', 'normalized', ...
                'Color', c, 'FontSize', 8, 'FontWeight', 'bold');
        end
    end
end
