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
%  PRE-PROCESS: CALCULATE METRICS
%  ========================================================================

% Pseudocount for Log Ratios (1 spike / hour)
c = 1 / 3600;

% Check for required columns
reqVars = {'frBspk', 'ss_frBspk', 'frSspk', 'ss_frSspk'};
if ~all(ismember(reqVars, tbl.Properties.VariableNames))
    error('[MCU_RCVSPACE] Table is missing required columns: %s', strjoin(reqVars, ', '));
end

% Log Fold Change (Relative) -> "dBrst", "dSngl"
% log(Post / Pre)
tbl.log_Brst = log((tbl.ss_frBspk + c) ./ (tbl.frBspk + c));
tbl.log_Sngl = log((tbl.ss_frSspk + c) ./ (tbl.frSspk + c));

% Absolute Difference (Absolute) -> "sl_Brst", "sl_Sngl"
% synLog(Post - Pre)
tbl.abs_Brst = tbl.ss_frBspk - tbl.frBspk;
tbl.abs_Sngl = tbl.ss_frSspk - tbl.frSspk;

tbl.sl_Brst = symlog(tbl.abs_Brst);
tbl.sl_Sngl = symlog(tbl.abs_Sngl);

%% ========================================================================
%  PREP
%  ========================================================================

grps = {'Control', 'MCU-KO'};

% Marker Size (Fixed)
tbl.scatSz = repmat(25, height(tbl), 1);

% Color Data
% Row 1: Burstiness
cData_Brst = tbl.dpBspk;
cLabel_Brst = 'Burstiness (logit)';
clim_Brst = [min(cData_Brst, [], 'all'), max(cData_Brst, [], 'all')];
if diff(clim_Brst) == 0; clim_Brst = clim_Brst + [-0.1 0.1]; end

% Row 2: dFr
cData_dFr = tbl.dFr;
cLabel_dFr = 'dFr (Fold Change)';
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
    
    plot_stateSpace(tbl, grp, 'log_Brst', 'log_Sngl', cData_Brst, cLabel_Brst);
    
    title(sprintf('%s', grp));
    xlabel('Burst (Log-Fold)');
    ylabel('Single (Log-Fold)');
    clim(clim_Brst);
end
c = colorbar; 
c.Label.String = cLabel_Brst;
c.Label.FontWeight = 'bold';


% --- ROW 2: ABSOLUTE DIFFERENCE ---
for iGrp = 1:2
    grp = grps{iGrp};
    nexttile;
    
    plot_stateSpace(tbl, grp, 'sl_Brst', 'sl_Sngl', cData_dFr, cLabel_dFr);
    
    xlabel('Burst (\DeltaHz SymLog)');
    ylabel('Single (\DeltaHz SymLog)');
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

function plot_stateSpace(tbl, grp, xVar, yVar, cData, ~)
    
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
    maxVal = max(abs(allVals)) * 1.1;
    if maxVal == 0; maxVal = 1; end
       
    % Filter NaNs
    xData = subTbl.(xVar);
    yData = subTbl.(yVar);
    mk = ~isnan(xData) & ~isnan(yData);
    xData = xData(mk);
    yData = yData(mk);
    curCData = cData(strcmpi(string(tbl.Group), grp));
    curCData = curCData(mk);
    
    % Reference Lines
    % Identity
    plot([-maxVal maxVal], [-maxVal maxVal], 'k:', 'LineWidth', 1.5);
    % Quadrant Separators
    xline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    yline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    
    % Scatter Plot
    scatter(xData, yData, subTbl.scatSz(mk), curCData, ...
        'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');

    xlim([-maxVal maxVal]);
    ylim([-maxVal maxVal]);
    set(gca, 'FontSize', 12);

    % Orthogonal Regression
    plot_linReg(xData, yData, 'hAx', gca, 'type', 'ortho', ...
        'clr', 'k', 'flgTxt', true);
end




