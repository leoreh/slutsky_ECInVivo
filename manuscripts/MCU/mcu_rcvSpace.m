function hFig = mcu_rcvSpace(tbl, varargin)
% MCU_RCVSPACE Plots Recovery State Space Trajectories & Strategy Angles.
%
%   hFig = MCU_RCVSPACE(tbl) plots the State Space Trajectories (Recovery Vectors)
%   comparing Burst Spike Gain vs Single Spike Gain.
%
%   INPUTS:
%       tbl     - (table) Data table containing recovery metrics.
%                 Must contain: 'Group', 'Name', 'fr', 'pBspk'.
%                 Should contain: 'dBrst', 'dSngl'. 
%
%       varargin - (param/value) Optional parameters (currently none used, 
%                  but kept for extensibility/consistency).
%
%   OUTPUTS:
%       hFig    - (Handle) Figure handle.
%
%   See also: MEA_RCVVEC, MCU_FRQQ, LME_ANALYSE

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'tbl', @istable);
parse(p, tbl, varargin{:});


%% ========================================================================
%  PRE-PROCESS: CALCULATE VECTORS
%  ========================================================================

% Cartesian to Polar
% NOTE: Inverted Axis -> X = dBrst, Y = dSngl
[tbl.vecTheta, tbl.vecR] = cart2pol(tbl.dBrst, tbl.dSngl);
tbl.vecDeg = rad2deg(tbl.vecTheta);

%% ========================================================================
%  PLOTTING: STATE SPACE & ANGLES
%  ========================================================================

hFig = figure('Position', [50 100 1000 800], 'Color', 'w', 'Name', 'Recovery State Space');
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Color Palette
cMap = containers.Map;
cMap('Control') = [0.2 0.2 0.2]; % Dark Gray
cMap('MCU-KO')  = [0.8 0.0 0.0]; % Red

grps = {'Control', 'MCU-KO'};
cols = {cMap('Control'), cMap('MCU-KO')};

% Max axis range for symmetry
maxVal = max(abs([tbl.dSngl; tbl.dBrst])) * 1.1;
if maxVal == 0; maxVal = 1; end % Safety

% Global Scaling for Plotting (Marker Size)
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

% Color by pBspk if available
cData = tbl.dpBspk;
cLabel = 'Burstiness (logit)';


%% ROW 1: SCATTER PLOTS (State Space)
for iGrp = 1:2
    nexttile;
    hold on;
    axis square;
    grid on;
    
    grp = grps{iGrp};
    subTbl = tbl(strcmpi(string(tbl.Group), grp), :);
    
    if isempty(subTbl)
        title(sprintf('%s (n=0)', grp));
        continue;
    end

    % --- Orthogonal Regression (PCA) ---
    % X = dBrst, Y = dSngl
    dataXY = [subTbl.dBrst, subTbl.dSngl];
    [coeff, ~, ~] = pca(dataXY);
    v1 = coeff(:, 1); % Direction of max variance
    mu = mean(dataXY);
    
    % Slope & Intercept
    % y - muY = m(x - muX)
    slope = v1(2) / v1(1);
    yInt  = mu(2) - slope * mu(1);
    
    % Plot Regression Line
    fLine = @(x) slope * x + yInt;
    fplot(fLine, [-maxVal maxVal], 'k-', 'LineWidth', 2);
    
    % Identity Line (y=x)
    plot([-maxVal maxVal], [-maxVal maxVal], 'k:', 'LineWidth', 1.5, 'DisplayName', 'Identity');
    
    % Quadrants
    xline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    yline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

    % Display Slope/Angle
    regAngle = atan2d(v1(2), v1(1));
    text(-maxVal*0.9, maxVal*0.9, sprintf('Slope: %.2f (%.1f\\circ)', slope, regAngle), ...
        'FontSize', 10, 'FontWeight', 'bold');
      
    % Scatter
    % X = dBrst, Y = dSngl
    scatter(subTbl.dBrst, subTbl.dSngl, subTbl.scatSz, cData(strcmpi(string(tbl.Group), grp)), ...
        'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');
    
    % Formatting
    xlim([-maxVal maxVal]);
    ylim([-maxVal maxVal]);
    xlabel('Burst FR Gain');
    ylabel('Single FR Gain');
    title(sprintf('%s (n=%d)', grp, height(subTbl)));
    
    colormap(gca, turbo);
    if iGrp == 2
        c = colorbar;
        c.Label.String = cLabel;
    end
    
    set(gca, 'FontSize', 12);
end


%% ROW 2: POLAR HISTOGRAMS (Strategies)
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
    
    % Mean Direction
    hold on;
    % Compute vector sum
    z = sum(exp(1i * subTbl.vecTheta));
    mAngle = angle(z);
    
    % For visualization, scaling the vector to the max count helps visibility
    % or just to the limit
    rLim = rlim;
    polarplot([0, mAngle], [0, rLim(2)], 'Color', cols{iGrp}, 'LineWidth', 3);
end

sgtitle(hFig, 'Recovery State Space & Strategies');

end
