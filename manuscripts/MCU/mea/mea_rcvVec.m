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
% Load with steady state variables
presets = {'steadyState'};
% [tbl, xVec, basepaths, v] = mcu_tblMea('presets', presets, 'flgOtl', true);

% Add logit pBspk
tblTrans = tbl_trans(tbl, 'varsInc', {'pBspk'}, 'logBase', 'logit');
tbl.pBspk_trans = tblTrans.pBspk;

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
%  PLOTTING: FIGURE 1 - RAW STATE SPACE & ANGLES
%  ========================================================================

fig1 = figure('Position', [50 100 1200 900], 'Color', 'w', 'Name', 'Raw Trajectories');
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Color Palette
cMap = containers.Map;
cMap('Control') = [0.2 0.2 0.2]; % Dark Gray
cMap('MCU-KO')  = [0.8 0.0 0.0]; % Red

colWt = cMap('Control');
colKo = cMap('MCU-KO');
grps = {'Control', 'MCU-KO'};
cols = {colWt, colKo};

% Max axis range for symmetry
maxVal = max(abs([tbl.dSngl; tbl.dBrst])) * 1.1;

% Global Scaling for Plotting
frLog = log10(tbl.fr);
lims  = prctile(frLog, [5 95]); 
scatNorm = (frLog - lims(1)) / diff(lims);
scatNorm(scatNorm < 0) = 0; scatNorm(scatNorm > 1) = 1; 
tbl.scatSz = scatNorm * 40 + 10; 

% --- ROW 1: SCATTER PLOTS ---
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
    
    % Equality Line (y=x)
    plot([-maxVal maxVal], [-maxVal maxVal], 'k:', 'LineWidth', 1.5);
    
    % Scatter (Color = Original Burstiness)
    % We use the transformed pBspk for clearer gradient
    scatter(subTbl.dSngl, subTbl.dBrst, subTbl.scatSz, subTbl.pBspk_trans, ...
        'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');
    
    % Mean Global Vector
    mSngl = mean(subTbl.dSngl);
    mBrst = mean(subTbl.dBrst);
    quiver(0, 0, mSngl, mBrst, 'Color', 'k', ...
        'AutoScale', 'off', 'LineWidth', 3, 'MaxHeadSize', 0.5);
    
    % Formatting
    xlim([-maxVal maxVal]);
    ylim([-maxVal maxVal]);
    xlabel('Log Fold Change (Single)');
    ylabel('Log Fold Change (Burst)');
    title(sprintf('%s (n=%d)', grp, height(subTbl)));
    
    colormap(gca, turbo);
    if iGrp == 2
        c = colorbar;
        c.Label.String = 'Baseline Burstiness (logit)';
    end
    
    set(gca, 'FontSize', 12);
end

% --- ROW 2: POLAR HISTOGRAMS ---
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
    mAngle = angle(sum(exp(1i * subTbl.vecTheta)));
    rLim = rlim;
    polarplot([0, mAngle], [0, rLim(2)], 'Color', cols{iGrp}, 'LineWidth', 3);
end

sgtitle(fig1, 'Figure 1: Recovery State Space & Strategies');



