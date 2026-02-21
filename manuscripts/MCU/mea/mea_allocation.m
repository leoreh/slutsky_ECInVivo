%% ========================================================================
%  PROPORTIONAL SPIKE ALLOCATION
%  ========================================================================
% Demonstrates the mathematical derivation of absolute firing rate recovery
% from the relative power-law relationship between single and burst spikes.
%
% This script empirically validates the mathematical relationship where
% a proportional (log-fold) recovery strategy in relative space dictates
% the absolute recovery trajectories based on initial baseline proportions.
%
% Derivation:
%   S1 = S0 * C * (B1 / B0)^beta
%   dS = S0 * [ C * (1 + dB / B0)^beta - 1 ]
%
% Initial Slope (Tangent at dB = 0):
%   d(dS)/d(dB) = C * beta * (S0 / B0)
%
% See "Proportional Spike Allocation.md" for full derivation details.
%
% See also: MCU_REGRESSION, MCU_RCVSPACE

%% ========================================================================
%  LOAD AND PRE-PROCESS DATA
%  ========================================================================

% Load with steady state variables
presets = {'steadyState'};
[tbl, ~, ~, ~] = mcu_tblMea('presets', presets, 'flgOtl', true);

% Relative (Log-Fold Change)
tbl.dBrst_rel = log((tbl.ss_frBspk) ./ (tbl.frBspk));
tbl.dSngl_rel = log((tbl.ss_frSspk) ./ (tbl.frSspk));

% Absolute Difference (Hz)
tbl.dBrst_abs = (tbl.ss_frBspk - tbl.frBspk);
tbl.dSngl_abs = (tbl.ss_frSspk - tbl.frSspk);


%% ========================================================================
%  1. RELATIVE SPACE (POWER LAW FIT)
%  ========================================================================
%  Extract the scaling exponent (beta) and constant (C) by fitting an 
%  orthogonal regression to the relative space log-fold changes.

fprintf('\n================================================================\n');
fprintf(' 1. RELATIVE SPACE FIT (POWER LAW)\n');
fprintf('================================================================\n');

grps = {'Control', 'MCU-KO'};
nGrps = length(grps);
allocModel = struct();

for iGrp = 1:nGrps
    grp = grps{iGrp};
    grpFld = strrep(grp, '-', '_');
    
    idx = tbl.Group == grp;
    
    xRel = tbl.dBrst_rel(idx);
    yRel = tbl.dSngl_rel(idx);
    
    % Remove invalid points
    mk = ~isnan(xRel) & ~isnan(yRel) & ~isinf(xRel) & ~isinf(yRel);
    xRel = xRel(mk);
    yRel = yRel(mk);
    
    % Fit Orthogonal Regression (PCA)
    % Log(S1/S0) = beta * Log(B1/B0) + Log(C)
    dataXY = [xRel, yRel];
    [coeff, ~, ~] = pca(dataXY);
    v1 = coeff(:, 1);
    mu = mean(dataXY);
    
    % Retrieve slope and intercept from 1st principal component
    beta = v1(2) / v1(1);
    logC = mu(2) - beta * mu(1);
    C = exp(logC);
    
    % Store
    allocModel.(grpFld).idx = idx;
    allocModel.(grpFld).beta = beta;
    allocModel.(grpFld).C = C;
    
    fprintf('--- %s ---\n', grp);
    fprintf('  Scaling Exponent (beta) : %6.3f\n', beta);
    fprintf('  Vertical Shift (C)      : %6.3f\n', C);
end


%% ========================================================================
%  2. ABSOLUTE SPACE (THEORETICAL PREDICTIONS VS EMPIRICAL)
%  ========================================================================
%  Demonstrate how the theoretical initial tangent calculated from baseline
%  parameters (S0, B0) predicts the absolute trajectory slopes.

fprintf('\n================================================================\n');
fprintf(' 2. ABSOLUTE RECOVERY SLOPES\n');
fprintf('================================================================\n');

for iGrp = 1:nGrps
    grp = grps{iGrp};
    grpFld = strrep(grp, '-', '_');
    
    idx = allocModel.(grpFld).idx;
    
    % Empirical Data for Absolute Fit
    xAbs = tbl.dBrst_abs(idx);
    yAbs = tbl.dSngl_abs(idx);
    
    S0 = tbl.frSspk(idx);
    B0 = tbl.frBspk(idx);
    
    % Remove invalid points
    mk = ~isnan(xAbs) & ~isnan(yAbs) & ~isinf(xAbs) & ~isinf(yAbs) & ...
         ~isnan(S0) & ~isnan(B0) & (B0 > 0);
    xAbs = xAbs(mk);
    yAbs = yAbs(mk);
    S0 = S0(mk);
    B0 = B0(mk);
    
    % Empirical Orthogonal Fit
    dataXY = [xAbs, yAbs];
    [coeff, ~, ~] = pca(dataXY);
    v1 = coeff(:, 1);
    empiricalSlope = v1(2) / v1(1);
    
    % Theoretical Tangents (evaluated at origin)
    C = allocModel.(grpFld).C;
    beta = allocModel.(grpFld).beta;
    
    tangentDist = C * beta * (S0 ./ B0);
    meanTangent = mean(tangentDist);
    tangentOfMeanBsl = C * beta * (mean(S0) / mean(B0));
    
    % Store
    allocModel.(grpFld).empiricalSlope = empiricalSlope;
    allocModel.(grpFld).tangentOfMeanBsl = tangentOfMeanBsl;
    
    fprintf('--- %s ---\n', grp);
    fprintf('  Empirical Fit (Orthogonal) : %6.3f\n', empiricalSlope);
    fprintf('  Predicted Initial Tangent  : %6.3f\n', tangentOfMeanBsl);
    %fprintf('  Mean of Per-Cell Tangents  : %6.3f\n', meanTangent);
end


%% ========================================================================
%  3. VISUALIZATION
%  ========================================================================

hFig = figure('Name', 'Proportional Spike Allocation', ...
    'Position', [50 50 1400 450], 'Color', 'w');
tTile = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Panel A: Relative Space (Power-law Conservation) ---
nexttile;
plot_scat(tbl, 'dBrst_rel', 'dSngl_rel', 'g', 'Group', ...
    'fitType', 'ortho', 'flgStats', true, 'alpha', 0.6);
plot_lineEq('hAx', gca, 'flgSqr', true);
xlabel('Burst Gain (log-fold)');
ylabel('Single Gain (log-fold)');
title('A. Relative Space (Power Law)');
legend('Location', 'northwest');


% --- Panel B: Theoretical Model (Absolute Prediction) ---
nexttile;
hold on;
clrs = lines(2); % Default colors matching plot_scat Auto mode

% Plot Empirical References as Background Shadows
for iGrp = 1:nGrps
    grp = grps{iGrp};
    grpFld = strrep(grp, '-', '_');
    clr = clrs(iGrp, :);
    idx = allocModel.(grpFld).idx;
    
    scatter(tbl.dBrst_abs(idx), tbl.dSngl_abs(idx), 20, clr, 'filled', ...
        'MarkerFaceAlpha', 0.15, 'HandleVisibility', 'off');
end

% Plot Trajectories
for iGrp = 1:nGrps
    grp = grps{iGrp};
    grpFld = strrep(grp, '-', '_');
    clr = clrs(iGrp, :);
    
    idx = allocModel.(grpFld).idx;
    
    S0_mean = mean(tbl.frSspk(idx), 'omitnan');
    B0_mean = mean(tbl.frBspk(idx), 'omitnan');
    
    C = allocModel.(grpFld).C;
    beta = allocModel.(grpFld).beta;
    
    % Extrapolate delta burst range based on limits
    minDB = min(tbl.dBrst_abs(idx), [], 'omitnan');
    maxDB = max(tbl.dBrst_abs(idx), [], 'omitnan');
    
    % Expand slightly for visualization, but avoid getting dB < -B0
    dbRng = maxDB - minDB;
    minEval = max(minDB - 0.1*dbRng, -B0_mean + 1e-4);
    dB = linspace(minEval, maxDB + 0.1*dbRng, 100); 
    
    % Equation: dS = S0 * [C * (1 + dB/B0)^beta - 1]
    dS_curve = S0_mean * ( C * (1 + dB ./ B0_mean).^beta - 1 );
    
    % Equation: dS = T * dB   where T = C*beta*(S0/B0)
    T = allocModel.(grpFld).tangentOfMeanBsl;
    dS_tang = T .* dB;
    
    % Plot Theoretical Curve
    plot(dB, dS_curve, '-', 'Color', clr, 'LineWidth', 2.5, ...
        'DisplayName', sprintf('%s (Predicted \\DeltaS)', grp));
    
    % Plot Initial Tangent at dB = 0
    plot(dB, dS_tang, '--', 'Color', clr, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('%s (Initial Tangent)', grp));
end

plot_lineEq('hAx', gca, 'flgSqr', true);
xlabel('\Delta Burst (Hz)', 'Interpreter', 'tex');
ylabel('\Delta Single (Hz)', 'Interpreter', 'tex');
title('B. Allocation Model Prediction');
legend('Location', 'northwest');


% --- Panel C: Empirical Absolute Space ---
nexttile;
plot_scat(tbl, 'dBrst_abs', 'dSngl_abs', 'g', 'Group', ...
    'fitType', 'ortho', 'flgStats', true, 'alpha', 0.6);
plot_lineEq('hAx', gca, 'flgSqr', true);
xlabel('\Delta Burst (Hz)', 'Interpreter', 'tex');
ylabel('\Delta Single (Hz)', 'Interpreter', 'tex');
title('C. Empirical Trajectories (Absolute)');
legend('Location', 'northwest');

