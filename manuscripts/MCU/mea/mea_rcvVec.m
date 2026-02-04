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
[tbl, xVec, basepaths, v] = mcu_tblMea('presets', presets, 'flgOtl', true);

% Add logit pBspk
tblTrans = tbl_trans(tbl, 'varsInc', {'pBspk'}, 'logBase', 'logit');
tbl.pBspk_trans = tblTrans.pBspk;


flgPseudo = false; % Validate "Pseudo-Tracking" (Quantile Matching) on MEA data
flgQq = false;

if flgQq
    hFig = mcu_frQq(tbl, 'dataSet', 'mea');
end

%% ========================================================================
%  PRE-PROCESS: CALCULATE VECTORS
%  ========================================================================

% 1. Calculate Deltas (Steady State - Baseline)
% ---------------------------------------------

if flgPseudo
    % --- CONVERT TO LONG FORMAT ---
    % match_qntl expects: [Group, Name, Day, vars...]
    
    % Baseline Table
    varsBsl = {'fr', 'frBspk', 'frSspk', 'pBspk'};
    tBsl = tbl(:, [{'Group', 'Name'}, varsBsl]);
    tBsl.Day = repmat({'BSL'}, height(tBsl), 1);
    tBsl.Day = categorical(tBsl.Day);

    % Steady State Table (Rename ss_ vars)
    varsSs = {'ss_fr', 'ss_frBspk', 'ss_frSspk', 'pBspk'}; % pBspk roughly same, or use ss_pBspk if available
    % Note: mcu_tblMea provides ss_pBspk? Let's assume we map standard vars.
    % Actually, let's just grab the ss columns and rename.
    tSs = tbl(:, {'Group', 'Name', 'frSs', 'ss_frBspk', 'ss_frSspk', 'ss_pBspk'});
    tSs = renamevars(tSs, {'frSs', 'ss_frBspk', 'ss_frSspk', 'ss_pBspk'}, ...
                           {'fr',    'frBspk',    'frSspk',    'pBspk'});
    tSs.Day = repmat({'BAC3'}, height(tSs), 1); % match_qntl hardcodes 'BAC3' as the 2nd day
    tSs.Day = categorical(tSs.Day);

    % Combine
    tblLong = [tBsl; tSs];
    
    % --- RUN MATCHING ---
    nBins = 80;
    fprintf('VALIDATION: Running match_qntl on MEA data (nBins=%d)...\n', nBins);
    tblSynth = match_qntl(tblLong, nBins, 'flgPool', true, ...
        'vars', {'fr', 'frBspk', 'frSspk', 'pBspk'});
    
    % --- MAP BACK TO MEA FORMAT ---
    % MEA Script Expects: frBspk (BSL), ss_frBspk (SteadyState)
    
    tbl = tblSynth;
    
    % Baseline Mappings
    tbl.fr        = tbl.fr_BSL;      % For Size Scaling
    tbl.pBspk     = tbl.pBspk_BSL;   % For Color
    tbl.frBspk    = tbl.frBspk_BSL;
    tbl.frSspk    = tbl.frSspk_BSL;
    
    % Steady State Mappings
    tbl.ss_frBspk = tbl.frBspk_BAC3;
    tbl.ss_frSspk = tbl.frSspk_BAC3;
    
    fprintf('VALIDATION: Replaced real units with %d Synthetic Units.\n', height(tbl));
end

% frBspk = Baseline Burst Rate
% ss_frBspk = Steady State Burst Rate
tbl.dBrst = tbl.ss_frBspk - tbl.frBspk;
tbl.dSngl = tbl.ss_frSspk - tbl.frSspk;

% Pseudocount
c = 1 / 3600;          

tbl.dBrst = log((tbl.ss_frBspk + c) ./ (tbl.frBspk + c));
tbl.dSngl = log((tbl.ss_frSspk + c) ./ (tbl.frSspk + c));
tbl.dFr = log((tbl.frSs) ./ (tbl.fr));

% 2. Calculate Vector Properties (Polar Coordinates)
% --------------------------------------------------
% Convert Cartesian (dSngl, dBrst) to Polar (Theta, Radius)
[tbl.vecTheta, tbl.vecR] = cart2pol(tbl.dSngl, tbl.dBrst);

% Convert Theta to Degrees for easier interpretation
tbl.vecDeg = rad2deg(tbl.vecTheta);


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
    
    % Orthogonal Regression (PCA)
    % ---------------------------
    % The first principal component is the line of total least squares
    dataXY = [subTbl.dSngl, subTbl.dBrst];
    [coeff, ~, latent] = pca(dataXY);
    v1 = coeff(:, 1); % Direction of max variance
    mu = mean(dataXY);
    
    % Calculate Slope & Intercept for plotting
    % y - muY = m(x - muX)  =>  y = m*x + (muY - m*muX)
    slope = v1(2) / v1(1);
    yInt  = mu(2) - slope * mu(1);
    
    % Plot Regression Line
    fLine = @(x) slope * x + yInt;
    fplot(fLine, [-maxVal maxVal], 'k-', 'LineWidth', 2);
    
    % Equality Line (y=x)
    plot([-maxVal maxVal], [-maxVal maxVal], 'k:', 'LineWidth', 1.5);
    
    % Quadrant Lines
    xline(0, 'k--', 'LineWidth', 1);
    yline(0, 'k--', 'LineWidth', 1);

    % Display Slope/Angle
    regAngle = atan2d(v1(2), v1(1));
    text(-maxVal*0.9, maxVal*0.9, sprintf('Slope: %.2f (%.1f\\circ)', slope, regAngle), ...
        'FontSize', 10, 'FontWeight', 'bold');
      
    % Scatter (Color = Original Burstiness)
    % We use the transformed pBspk for clearer gradient
    scatter(subTbl.dSngl, subTbl.dBrst, subTbl.scatSz, subTbl.pBspk, ...
        'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');
    
    % Mean Global Vector
    % mSngl = mean(subTbl.dSngl);
    % mBrst = mean(subTbl.dBrst);
    % quiver(0, 0, mSngl, mBrst, 'Color', 'k', ...
    %     'AutoScale', 'off', 'LineWidth', 3, 'MaxHeadSize', 0.5);
    
    % Formatting
    xlim([-maxVal maxVal]);
    ylim([-maxVal maxVal]);
    xlabel('Single FR Gain');
    ylabel('Busrt FR Gain');
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


%% ========================================================================
%  ANALYSIS: ROTATED COORDINATES (D-METRIC)
%  ========================================================================
%  D = dBurst - dSingle
%  Deviations from the identity line (y=x).
%  D > 0: Burst Biased Strategy
%  D = 0: Balanced Recovery
%  D < 0: Single Biased Strategy

fprintf('\n================================================================\n');
fprintf(' STATISTICS: D-METRIC STRATEGY ANALYSIS\n');
fprintf('================================================================\n');

% 1. Calculate D
tbl.D = tbl.dBrst - tbl.dSngl;

% 2. LME Model
% Model: D ~ Group + (1|Name)
% Intercept represents the Control deviation from 0.
% Group coefficient represents the shift in MCU-KO.
frml = 'D ~ Group * (pBspk + fr) + (1|Name)';
frml = 'D ~ Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, 'dist', 'normal');

frml = 'dBrst ~ dSngl * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');

frml = 'dSngl ~ dBrst * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');

frml = 'dSngl ~ (fr + pBspk + dBrst) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');

frml = 'dBrst ~ (fr + pBspk) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');

frml = 'dFr ~ (dBrst + dSngl) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');

frml = 'frSs ~ (fr + pBspk) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'gamma');

% tblGUI_scatHist(tbl, 'xVar', 'pBspk_trans', 'yVar', 'dBrst', 'grpVar', 'Group');
% tblGUI_bar(tbl, 'yVar', 'pBspk', 'xVar', 'Group');
