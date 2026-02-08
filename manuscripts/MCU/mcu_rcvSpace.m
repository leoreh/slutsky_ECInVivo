function [hFigScat, hFigPol] = mcu_rcvSpace(tbl, varargin)
% MCU_RCVSPACE Plots Recovery State Space Trajectories & Strategy Angles.
%
%   [hFigScat, hFigPol] = MCU_RCVSPACE(tbl) plots the State Space Trajectories
%   (Recovery Vectors) comparing Burst Spike Gain vs Single Spike Gain.
%
%   It produces two figures:
%       1. Scatter Plots (hFigScat): 2x2 Layout
%          Row 1: Log Fold Change (Relative)
%          Row 2: Absolute Difference (SynLog transformed)
%
%       2. Polar Histograms (hFigPol): 1x2 Layout
%          showing the distribution of recovery angles based on Log Fold Change.
%
%   INPUTS:
%       tbl     - (table) Data table containing recovery metrics.
%                 Must contain: 'Group', 'Name', 'fr', 'pBspk'.
%                 Must contain for calcs: 'frBspk', 'ss_frBspk', 'frSspk', 'ss_frSspk'
%
%   OUTPUTS:
%       hFigScat - (Handle) Figure handle for Scatter plots.
%       hFigPol  - (Handle) Figure handle for Polar plots.
%
%   See also: MEA_RCVVEC, MCU_FRQQ, LME_ANALYSE

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
parse(p, tbl, varargin{:});

%% ========================================================================
%  PREP
%  ========================================================================

grps = {'Control', 'MCU-KO'};

% Marker Size (Fixed)
tbl.scatSz = repmat(25, height(tbl), 1);

% Color Data
% Row 1: Burstiness
cData_Brst = tbl.dpBspk;
cLabel_Brst = '\Delta Burstiness (log-odds)';
clim_Brst = [min(cData_Brst, [], 'all'), max(cData_Brst, [], 'all')];
if diff(clim_Brst) == 0; clim_Brst = clim_Brst + [-0.1 0.1]; end

% Row 2: dFr
cData_dFr = tbl.dFr;
cLabel_dFr = '\Delta FR (Fold Change)';
clim_dFr = [min(cData_dFr, [], 'all'), max(cData_dFr, [], 'all')];
if diff(clim_dFr) == 0; clim_dFr = clim_dFr + [-0.1 0.1]; end



%% ========================================================================
%  FIGURE: SCATTER PLOTS (State Space)
%  ========================================================================
%   Row 1: Log Fold Change
%   Row 2: Absolute Difference (SynLog)

hFigScat = figure('Position', [50 50 800 800], 'Color', 'w', ...
    'Name', 'Recovery State Space: Scatter');
hTile = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- ROW 1: LOG FOLD CHANGE ---
for iGrp = 1:2
    grp = grps{iGrp};
    nexttile;
    
    plot_stateSpace(tbl, grp, 'dBrst_rel', 'dSngl_rel', cData_Brst, cLabel_Brst, false);
    
    title(sprintf('%s', grp));
    xlabel('Burst Gain (log-fold)');
    ylabel('Single Gain (log-fold)');
    clim(clim_Brst);
end
c = colorbar; 
c.Label.String = cLabel_Brst;
c.Label.FontWeight = 'bold';


% --- ROW 2: ABSOLUTE DIFFERENCE ---
for iGrp = 1:2
    grp = grps{iGrp};
    nexttile;
    
    plot_stateSpace(tbl, grp, 'dBrst_abs', 'dSngl_abs', cData_dFr, cLabel_dFr, true);
    
    xlabel('\Delta Burst (Hz)', 'Interpreter', 'tex');
    ylabel('\Delta Single (Hz)', 'Interpreter', 'tex');
    clim(clim_dFr);
end
c = colorbar; 
c.Label.String = cLabel_dFr;
c.Label.FontWeight = 'bold';


% Global Colormap
colormap(hFigScat, turbo);


end


%% ========================================================================
%  HELPER: PLOT STATE SPACE
%  ========================================================================

function plot_stateSpace(tbl, grp, xVar, yVar, cData, ~, isSymlog)
    
    if nargin < 7; isSymlog = false; end
    
    subTbl = tbl(strcmpi(string(tbl.Group), grp), :);
    
    hold on;
    axis square;
    grid off;
    
    if isempty(subTbl)
        text(0,0,'n=0', 'HorizontalAlignment','center');
        return; 
    end

    % Determine Axis Limits (Symmertic)
    % We want 0,0 in center.
    allVals = [tbl.(xVar); tbl.(yVar)];
    maxAbs = max(abs(allVals));
    if maxAbs == 0; maxAbs = 1; end
    
    if isSymlog
        % Pre-calc maxVal in symlog space for limits
        % Note: we apply padding AFTER transformation to ensure visual margin
        limVal = symlog(maxAbs) * 1.1; 
    else
        limVal = maxAbs * 1.1;
    end
       
    % Filter NaNs
    xData = subTbl.(xVar);
    yData = subTbl.(yVar);
    mk = ~isnan(xData) & ~isnan(yData);
    xData = xData(mk);
    yData = yData(mk);
    curCData = cData(strcmpi(string(tbl.Group), grp));
    curCData = curCData(mk);    
    
    % Scatter Plot
    scatter(xData, yData, subTbl.scatSz(mk), curCData, ...
        'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');
  
    set(gca, 'FontSize', 12);

    if isSymlog
        % 1. Transform Axes (Scatter points become SymLog)
        symlog(gca, 'xy');
        
        % 2. Transform Data for Regression (to match visual)
        xReg = symlog(xData);
        yReg = symlog(yData);
    else
        xReg = xData;
        yReg = yData;
    end
    
    xlim([-limVal limVal]);
    ylim([-limVal limVal]);

    % Orthogonal Regression
    % If symlog, this plots a straight line in the transformed space (visually correct)
    % The Slope/Angle calculated will be for the transformed data.
    plot_linReg(xReg, yReg, 'hAx', gca, 'type', 'ortho', ...
        'clr', 'k', 'flgTxt', true);

    % Reference Lines (Visual Space)
    % Identity
    plot([-limVal limVal], [-limVal limVal], 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    % Quadrant Separators
    xline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    yline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

end




