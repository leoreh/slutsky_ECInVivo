function hFigScat = mcu_rcvSpace(tbl, varargin)
% MCU_RCVSPACE Plots Recovery State Space Trajectories & Strategy Angles.
%
%   hFigScat = MCU_RCVSPACE(tbl) plots the State Space Trajectories
%   (Recovery Vectors) comparing Burst Spike Gain vs Single Spike Gain.
%
%   It produces a standard scatter plot figure (hFigScat): 2x2 Layout
%      Row 1: Log Fold Change (Relative)
%      Row 2: Absolute Difference
%
%   INPUTS:
%       tbl     - (table) Data table containing recovery metrics.
%                 Must contain: 'Group', 'Name', 'fr', 'pBspk'.
%                 Must contain for calcs: 'frBspk', 'ss_frBspk', 'frSspk', 'ss_frSspk'
%
%   OUTPUTS:
%       hFigScat - (Handle) Figure handle for Scatter plots.
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
%   Row 2: Absolute Difference

hFigScat = figure('Position', [50 50 800 800], 'Color', 'w', ...
    'Name', 'Recovery State Space: Scatter');
hTile = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- ROW 1: LOG FOLD CHANGE ---
for iGrp = 1:2
    grp = grps{iGrp};
    nexttile;
    
    % Filter Data
    mask = strcmpi(string(tbl.Group), grp);
    subTbl = tbl(mask, :);
    subCData = cData_Brst(mask);
    
    % Plot Scatter (via plot_scat)
    plot_scat(subTbl, 'dBrst_rel', 'dSngl_rel', ...
        'c', subCData, ...
        'hAx', gca, ...
        'fitType', 'ortho', ...
        'flgStats', true, ...
        'alpha', 0.6);

    % Plot Equality & Zero Lines (and square axis)
    plot_lineEq('hAx', gca, 'flgSqr', true);
    
    legend('Location', 'northwest')
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
    
    % Filter Data
    mask = strcmpi(string(tbl.Group), grp);
    subTbl = tbl(mask, :);
    subCData = cData_dFr(mask);

    % Plot Scatter (via plot_scat)
    plot_scat(subTbl, 'dBrst_abs', 'dSngl_abs', ...
        'c', subCData, ...
        'hAx', gca, ...
        'fitType', 'ortho', ...
        'flgStats', true, ...
        'alpha', 0.6);

    % Plot Equality & Zero Lines (and square axis)
    plot_lineEq('hAx', gca, 'flgSqr', true);
    
    legend('Location', 'northwest')
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







