function [hFigScat, hFigPol] = mcu_rcvSpace(tbl, varargin)
% MCU_RCVSPACE Plots Recovery State Space Trajectories & Strategy Angles.
%
%   [hFigScat, hFigPol] = MCU_RCVSPACE(tbl) plots the State Space Trajectories
%   (Recovery Vectors) comparing Burst Spike Gain vs Single Spike Gain.
%
%   It produces two figures:
%       1. Scatter Plots (hFigScat): 2x2 Layout
%          Row 1: Log Fold Change (Relative)
%          Row 2: Absolute Difference (asinh transformed)
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
%  PRE-PROCESS: CALCULATE METRICS
%  ========================================================================

% Pseudocount for Log Ratios (1 spike / hour)
c = 1 / 3600;

% Check for required columns
reqVars = {'frBspk', 'ss_frBspk', 'frSspk', 'ss_frSspk'};
if ~all(ismember(reqVars, tbl.Properties.VariableNames))
    error('[MCU_RCVSPACE] Table is missing required columns: %s', strjoin(reqVars, ', '));
end

% 1. Log Fold Change (Relative) -> "dBrst", "dSngl"
% ------------------------------------------------
% log(Post / Pre)
tbl.log_Brst = log((tbl.ss_frBspk + c) ./ (tbl.frBspk + c));
tbl.log_Sngl = log((tbl.ss_frSspk + c) ./ (tbl.frSspk + c));

% 2. Absolute Difference (Absolute) -> "ahs_Brst", "ahs_Sngl"
% ------------------------------------------------
% asinh(Post - Pre)
tbl.abs_Brst = tbl.ss_frBspk - tbl.frBspk;
tbl.abs_Sngl = tbl.ss_frSspk - tbl.frSspk;

tbl.ahs_Brst = asinh(tbl.abs_Brst);
tbl.ahs_Sngl = asinh(tbl.abs_Sngl);

% 3. Cartesian to Polar (Based on Log Fold Change)
% ------------------------------------------------
% NOTE: Inverted Axis -> X = log_Brst, Y = log_Sngl
[tbl.vecTheta, tbl.vecR] = cart2pol(tbl.log_Brst, tbl.log_Sngl);
tbl.vecDeg = rad2deg(tbl.vecTheta);


%% ========================================================================
%  PREP: GROUPS & PALETTE
%  ========================================================================

cMap = containers.Map;
cMap('Control') = [0.2 0.2 0.2]; % Dark Gray
cMap('MCU-KO')  = [0.8 0.0 0.0]; % Red

grps = {'Control', 'MCU-KO'};
cols = {cMap('Control'), cMap('MCU-KO')};

% Marker Size Scaling (Global)
if ismember('fr', tbl.Properties.VariableNames)
    frLog = log10(tbl.fr);
    lims  = prctile(frLog, [5 95]); 
    if diff(lims) == 0; lims = [lims(1)-0.1, lims(2)+0.1]; end
    scatNorm = (frLog - lims(1)) / diff(lims);
    scatNorm(scatNorm < 0) = 0; scatNorm(scatNorm > 1) = 1; 
    tbl.scatSz = scatNorm * 40 + 10; 
else
    tbl.scatSz = repmat(20, height(tbl), 1);
end

% Color by pBspk (for scatter points)
cData = tbl.pBspk; 
% Ensure pBspk is transformed if needed, but existing code used 'pBspk_trans'. 
% If 'pBspk_trans' exists, use it, otherwise use pBspk.
if ismember('pBspk_trans', tbl.Properties.VariableNames)
    cData = tbl.pBspk_trans;
    cLabel = 'Burstiness (logit)';
else
    cLabel = 'Burstiness';
end


%% ========================================================================
%  FIGURE 1: SCATTER PLOTS (State Space)
%  ========================================================================
%   Row 1: Log Fold Change
%   Row 2: Absolute Difference (asinh)

hFigScat = figure('Position', [50 50 800 800], 'Color', 'w', ...
    'Name', 'Recovery State Space: Scatter');
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- ROW 1: LOG FOLD CHANGE ---
for iGrp = 1:2
    grp = grps{iGrp};
    nexttile;
    
    plot_state_space(tbl, grp, 'log_Brst', 'log_Sngl', cData, cLabel);
    
    title(sprintf('%s (Log Fold)', grp));
    xlabel('Burst Log-Fold');
    ylabel('Single Log-Fold');
end

% --- ROW 2: ABSOLUTE DIFFERENCE ---
for iGrp = 1:2
    grp = grps{iGrp};
    nexttile;
    
    plot_state_space(tbl, grp, 'ahs_Brst', 'ahs_Sngl', cData, cLabel);
    
    title(sprintf('%s (Abs Diff)', grp));
    xlabel('Burst \DeltaHz (asinh)');
    ylabel('Single \DeltaHz (asinh)');
end

% Add Colorbar to the last tile of each row or just once?
% Doing it per row helps if ranges differ, but here cData is same (Burstiness).
% Let's add it to the last tile plotted.
c = colorbar;
c.Label.String = cLabel;
colormap(hFigScat, turbo);


%% ========================================================================
%  FIGURE 2: POLAR HISTOGRAMS
%  ========================================================================

hFigPol = figure('Position', [900 100 800 400], 'Color', 'w', ...
    'Name', 'Recovery Strategies: Polar');
tP = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for iGrp = 1:2
    nexttile;
    grp = grps{iGrp};
    subTbl = tbl(strcmpi(string(tbl.Group), grp), :);
    
    if isempty(subTbl); continue; end
    
    % Polar Histogram
    polarhistogram(subTbl.vecTheta, 20, 'FaceColor', cols{iGrp}, ...
        'FaceAlpha', 0.6, 'EdgeColor', 'w');
    
    title(sprintf('%s Directions', grp));
    set(gca, 'FontSize', 12);
    
    % Mean Direction Vector
    hold on;
    % Compute vector sum (Resultant Vector)
    z = sum(exp(1i * subTbl.vecTheta));
    mAngle = angle(z);
    
    % Scale vector for visibility (to the limit of the plot)
    rLim = rlim;
    polarplot([0, mAngle], [0, rLim(2)], 'Color', cols{iGrp}, 'LineWidth', 3);
end

sgtitle(hFigPol, 'Recovery Strategies (Log Fold Angles)');


end


%% ========================================================================
%  HELPER: PLOT STATE SPACE
%  ========================================================================

function plot_state_space(tbl, grp, xVar, yVar, cData, ~)
    
    subTbl = tbl(strcmpi(string(tbl.Group), grp), :);
    
    hold on;
    axis square;
    grid on;
    
    if isempty(subTbl)
        text(0,0,'n=0', 'HorizontalAlignment','center');
        return; 
    end

    % 1. Determine Axis Limits (Symmertic)
    % ------------------------------------------------
    % We want 0,0 in center.
    allVals = [tbl.(xVar); tbl.(yVar)];
    maxVal = max(abs(allVals)) * 1.1;
    if maxVal == 0; maxVal = 1; end
    
    % 2. Orthogonal Regression
    % ------------------------------------------------
    dataXY = [subTbl.(xVar), subTbl.(yVar)];
    
    % Exclude NaNs if any
    mk = all(~isnan(dataXY), 2);
    dataXY = dataXY(mk, :);
    subTbl = subTbl(mk, :);
    curCData = cData(strcmpi(string(tbl.Group), grp));
    curCData = curCData(mk);
    
    if size(dataXY, 1) > 2
        [coeff, ~, ~] = pca(dataXY);
        v1 = coeff(:, 1); 
        mu = mean(dataXY);
        
        slope = v1(2) / v1(1);
        yInt  = mu(2) - slope * mu(1);
        
        % Plot Regression Line
        fLine = @(x) slope * x + yInt;
        fplot(fLine, [-maxVal maxVal], 'k-', 'LineWidth', 2);
        
        % Annotate Slope
        regAngle = atan2d(v1(2), v1(1));
        text(-maxVal*0.9, maxVal*0.9, sprintf('Slope: %.2f (%.1f\\circ)', slope, regAngle), ...
            'FontSize', 10, 'FontWeight', 'bold');
    end
    
    % 3. Reference Lines
    % ------------------------------------------------
    % Identity
    plot([-maxVal maxVal], [-maxVal maxVal], 'k:', 'LineWidth', 1.5);
    % Quadrant Separators
    xline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    yline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    
    % 4. Scatter Plot
    % ------------------------------------------------
    scatter(subTbl.(xVar), subTbl.(yVar), subTbl.scatSz, curCData, ...
        'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');
    
    xlim([-maxVal maxVal]);
    ylim([-maxVal maxVal]);
    set(gca, 'FontSize', 12);
end
