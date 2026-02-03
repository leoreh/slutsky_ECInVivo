%% ========================================================================
%  MEA RECOVERY VECTOR ANALYSIS (State Space Trajectory)
%  ========================================================================
%  Analyzes the recovery "strategy" of neurons by decomposing firing rate
%  recovery into Single-Spike and Burst-Spike components.
%
%  Calculates a vector from Baseline (0,0) to Steady State (dSingle, dBurst)
%  and analyzes the Magnitude (Strength of change) and Angle (Strategy).

%% ========================================================================
%  LOAD DATA
%  ========================================================================

% Load with steady state variables
presets = {'steadyState'};
% [tbl, xVec, basepaths, v] = mcu_tblMea('presets', presets, 'flgOtl', true);

%% ========================================================================
%  PRE-PROCESS: CALCULATE VECTORS
%  ========================================================================

% 1. Calculate Deltas (Steady State - Baseline)
% ---------------------------------------------
% frBspk = Baseline Burst Rate
% ss_frBspk = Steady State Burst Rate
tbl.dBrst = tbl.ss_frBspk - tbl.frBspk;
tbl.dSngl = tbl.ss_frSspk - tbl.frSspk;

% Pseudocount
c = 1 / 3600;          

tbl.dBrst = log((tbl.ss_frBspk + c) ./ (tbl.frBspk + c));
tbl.dSngl = log((tbl.ss_frSspk + c) ./ (tbl.frSspk + c));

% 2. Calculate Vector Properties (Polar Coordinates)
% --------------------------------------------------
% Convert Cartesian (dSngl, dBrst) to Polar (Theta, Radius)
[tbl.vecTheta, tbl.vecR] = cart2pol(tbl.dSngl, tbl.dBrst);

% Convert Theta to Degrees for easier interpretation
tbl.vecDeg = rad2deg(tbl.vecTheta);

% 3. Calculate "Strategy" Categories via Quadrants
% ------------------------------------------------
% Q1: Indiscriminate Increase (+S, +B)
% Q2: Burst Compensation (-S, +B) -> Expected for Control
% Q3: Global Depression (-S, -B)
% Q4: Single Compensation (+S, -B)
tbl.Strategy = categorical(zeros(height(tbl), 1));

idxQ1 = tbl.dSngl > 0 & tbl.dBrst > 0;
idxQ2 = tbl.dSngl < 0 & tbl.dBrst > 0;
idxQ3 = tbl.dSngl < 0 & tbl.dBrst < 0;
idxQ4 = tbl.dSngl > 0 & tbl.dBrst < 0;

tbl.Strategy(idxQ1) = 'Global Inc';
tbl.Strategy(idxQ2) = 'Burst Comp';
tbl.Strategy(idxQ3) = 'Global Dec';
tbl.Strategy(idxQ4) = 'Sngl Comp';
tbl.Strategy = removecats(tbl.Strategy);

%% ========================================================================
%  STATISTICS: RECOVERY MAGNITUDE & DIRECTION
%  ========================================================================

fprintf('\n================================================================\n');
fprintf(' STATISTICS: RECOVERY VECTOR ANALYSIS\n');
fprintf('================================================================\n');

% 1. Vector Magnitude (LME)
% -------------------------
% Does the magnitude of change differ between groups?
% fprintf('\n[ANALYSIS] Vector Magnitude (Radius)\n');
% [lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, 'vecR ~ Group + (1|Name)');

% 2. Vector Angle (LME - Linear approximation)
% --------------------------------------------
% Note: Strictly speaking, angles are circular. However, if vectors are 
% directional (e.g., mostly pointing up/left), linear LME can be a proxy. 
% Check distribution first.
% fprintf('\n[ANALYSIS] Vector Angle (Degrees)\n');
% lme_analyse(tbl, 'vecDeg ~ Group + (1|Name)');

% 3. Strategy Distribution (Chi-Square)
% -------------------------------------
fprintf('\n[ANALYSIS] Strategy Distribution (Contingency Table)\n');
[statTbl, chi2, pVal] = crosstab(tbl.Group, tbl.Strategy);
disp('Contingency Table (Rows=Group, Cols=Strategy):');
disp(statTbl);
fprintf('Chi2: %.2f, p-value: %.4g\n', chi2, pVal);


%% ========================================================================
%  PLOTTING: STATE SPACE TRAJECTORIES
%  ========================================================================

fig = figure('Position', [100 100 1200 500], 'Color', 'w');
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Color Palette
cMap = containers.Map;
cMap('Control') = [0.2 0.2 0.2]; % Dark Gray
cMap('MCU-KO')  = [0.8 0.0 0.0]; % Red

colWt = cMap('Control');
colKo = cMap('MCU-KO');

% Groups
grps = {'Control', 'MCU-KO'};
cols = {colWt, colKo};

% Max axis range for symmetry
maxVal = max(abs([tbl.dSngl; tbl.dBrst])) * 1.1;

for iGrp = 1:2
    nexttile;
    hold on;
    axis square;
    grid on;
    
    grp = grps{iGrp};
    subTbl = tbl(tbl.Group == grp, :);
    
    % Draw Axes
    xline(0, 'k--', 'LineWidth', 1);
    yline(0, 'k--', 'LineWidth', 1);
    
    % Draw Quivers (Arrows)
    % Origin (0,0) to (dSngl, dBrst)
    % We use quiver(X, Y, U, V) where X,Y are origins (0,0) and U,V are components
    nU = height(subTbl);
    q = quiver(zeros(nU, 1), zeros(nU, 1), subTbl.dSngl, subTbl.dBrst, ...
        'Color', [cols{iGrp}, 0.3], 'AutoScale', 'off', 'LineWidth', 1);
    
    % Mean Vector
    mSngl = mean(subTbl.dSngl);
    mBrst = mean(subTbl.dBrst);
    qM = quiver(0, 0, mSngl, mBrst, 'Color', cols{iGrp}, ...
        'AutoScale', 'off', 'LineWidth', 3, 'MaxHeadSize', 0.5);
    
    % Formatting
    xlim([-maxVal maxVal]);
    ylim([-maxVal maxVal]);
    xlabel('Log Fold Change (Single)');
    ylabel('Log Fold Change (Burst)');
    title(sprintf('%s (n=%d)', grp, nU));
    
    % Add Quadrant Labels (Optional)
    text(-maxVal*0.8, maxVal*0.8, 'Burst Comp', 'Color', [.5 .5 .5]);
    text(maxVal*0.8, -maxVal*0.8, 'Sngl Comp', 'Color', [.5 .5 .5]);
    
    set(gca, 'FontSize', 12, 'LineWidth', 1.2);
end

sgtitle('Recovery State Space Trajectories');

%% ========================================================================
%  PLOTTING: POLAR HISTOGRAM (ANGLES)
%  ========================================================================

for iGrp = 1:2
    nexttile;
    grp = grps{iGrp};
    subTbl = tbl(tbl.Group == grp, :);
    
    % Polar Histogram
    polarhistogram(subTbl.vecTheta, 20, 'FaceColor', cols{iGrp}, ...
        'FaceAlpha', 0.6, 'EdgeColor', 'w');
    
    title(sprintf('%s Directions', grp));
    set(gca, 'FontSize', 12);
    
    % Add Mean Direction line
    hold on;
    mAngle = circ_mean(subTbl.vecTheta); % Requires CircStats or defined below
    rLim = rlim;
    polarplot([0, mAngle], [0, rLim(2)], 'Color', cols{iGrp}, 'LineWidth', 3);
end


%% ========================================================================
%  HELPER: CIRCULAR MEAN (Simple approximation if toolbox missing)
%  ========================================================================
function mu = circ_mean(alpha)
    % Computes mean direction for circular data
    % Input: alpha (radians)
    sa = sum(sin(alpha));
    ca = sum(cos(alpha));
    mu = atan2(sa, ca);
end
