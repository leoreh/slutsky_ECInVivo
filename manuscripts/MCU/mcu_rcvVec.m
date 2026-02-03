%% ========================================================================
%  MCU RECOVERY VECTOR ANALYSIS (In Vivo Pseudo-Tracking)
%  ========================================================================
%  Analyzes the recovery "strategy" of neurons in vivo by creating
%  "synthetic units" via quantile matching.
%
%  Since we cannot track the same units over days in vivo, we rank units
%  by firing rate in Baseline and BAC3 (Steady State) and pair them
%  based on their rank (Quantile Matching).
%
%  Then, we perform the same vector analysis as in MEA:
%  - Calculate a vector from Baseline (0,0) to Steady State (dSingle, dBurst)
%  - Analyze Magnitude (Strength) and Angle (Strategy)
%
%  ========================================================================



%% ========================================================================
%  LOAD DATA
%  ========================================================================

% Load table
% basepaths = [mcu_basepaths('wt'), mcu_basepaths('mcu')];
% presets = {'frNet', 'brst'};
% tbl = mcu_tblVivo('basepaths', basepaths, 'flgClean', true, ...
%     'presets', presets);

% Filter: RS units only
tblLme = tbl(tbl.unitType == 'RS', :);

% Filter: Only BSL and BAC3 days
idxDay = ismember(tblLme.Day, {'BSL', 'BAC3'});
tblLme = tblLme(idxDay, :);

% Assert no zero values for logs. This is instead of a pseudocount.
tblTrans = tbl_trans(tblLme, 'flg0', true, 'verbose', true);


%% ========================================================================
%  ASSERT EQUAL DISTRIBUTIONS & FILTERING
%  ========================================================================
%  Check if the data is consistent and follows expected distributions (Log-Normal).
%  Filter out low firing rate units that deviate from the distribution.

flgQq = false;
cutoff_Z = -1.5; % User defined cut-off (Visual Inspection)

if flgQq
    idxWt = tblLme.Group == 'Control';
    idxMcu = tblLme.Group == 'MCU-KO';
    
    frBsl = tblLme.fr(tblLme.Day == 'BSL' & idxWt);
    frBac = tblLme.fr(tblLme.Day == 'BAC3' & idxWt);
    
    frBslMcu = tblLme.fr(tblLme.Day == 'BSL' & idxMcu);
    
    hFig = figure('Position', [300 300 900 900], 'Color', 'w', 'Name', 'Data Distribution Check');
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Log-Normality Check (The most important for cutoff)
    nexttile;
    hQ = qqplot(log(frBsl)); 
    xlabel('Theoretical Normal Quantiles');
    ylabel('Sample Log(FR) Quantiles');
    title('Log-Normality Check (BSL Control)');
    grid on; axis square;   
    hQ(1).LineWidth = 1;

    % Group Check (Control vs MCU)
    nexttile;
    hQ = qqplot(log(frBsl), log(frBslMcu));
    xlabel('Control BSL Log(FR)');
    ylabel('MCU BSL Log(FR)');
    title('Group: Control vs MCU (BSL)');
    grid on; axis square;
    refline(1,0);
    hQ(1).LineWidth = 1;
    
    % Stability Check (BSL vs BAC3)
    nexttile;  hold on;
    hQ = qqplot(log(frBsl), log(frBac));
    xlabel('BSL Log(FR) Quantiles');
    ylabel('BAC3 Log(FR) Quantiles');
    title('Stability: BSL vs BAC3 (Control)');
    grid on; axis square;
    refline(1,0); 
    hQ(1).LineWidth = 1;

    % Calculate cutoff value in log space for plotting
    logCutoff = prctile(log(frBsl), normcdf(cutoff_Z)*100);
    yline(logCutoff, 'r--', 'Cutoff', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');

    % Stability Check (MCU BSL vs BAC3)
    % Just to be thorough
    frBacMcu = tblLme.fr(tblLme.Day == 'BAC3' & idxMcu);
    nexttile;
    hQ = qqplot(log(frBslMcu), log(frBacMcu));
    xlabel('MCU BSL Log(FR)');
    ylabel('MCU BAC3 Log(FR)');
    title('Stability: BSL vs BAC3 (MCU)');
    grid on; axis square;
    refline(1,0); 
    hQ(1).LineWidth = 1;
end

% --- FILTERING ---
cutoff_p = normcdf(cutoff_Z);
% Define threshold based on Control Baseline (canonical)
frBslRef = tblLme.fr(tblLme.Day == 'BSL' & tblLme.Group == 'Control'); 
cutoff_Hz = prctile(frBslRef, cutoff_p * 100);

fprintf('\n================================================================\n');
fprintf(' FILTERING LOW FR UNITS\n');
fprintf('================================================================\n');
fprintf('Cut-off Z-score: %.2f (%.1f%%)\n', cutoff_Z, cutoff_p*100);
fprintf('Cut-off FR     : %.4f Hz\n', cutoff_Hz);

% Apply Filter
nBefore = height(tblLme);

% Identify Units to keep (Must have BSL FR > cutoff)
goodIdx = tblLme.fr >= cutoff_Hz;

% Filter original table
tblLme = tblLme(goodIdx, :);

fprintf('Removed %d units (%.1f%%). Remaining: %d units\n', ...
    sum(~goodIdx), sum(~goodIdx)/length(goodIdx)*100, sum(goodIdx));
fprintf('================================================================\n\n');


%% ========================================================================
%  PSEUDO-TRACKING (QUANTILE MATCHING)
%  ========================================================================
%  Create 'Synthetic Units' by matching BSL and BAC3 units by rank.

nBins = 7; % Deciles

% Initialize container for synthetic table
tblSynth = match_qntl(tblLme, nBins, 'flgPool', false);

% tblGUI_scatHist(tblSynth, 'xVar', 'pBspk_BSL', 'yVar', 'fr_BAC3', 'grpVar', 'Group');


%% ========================================================================
%  VECTOR ANALYSIS
%  ========================================================================

% 1. Calculate Deltas
% -------------------
% dBrst = log(frBspk_BAC3 / frBspk_BSL)
% dSngl = log(frSspk_BAC3 / frSspk_BSL)

tblSynth.dBrst = log((tblSynth.frBspk_BAC3) ./ (tblSynth.frBspk_BSL));
tblSynth.dSngl = log((tblSynth.frSspk_BAC3) ./ (tblSynth.frSspk_BSL));

% 2. Calculate Vector Properties
% ------------------------------
[tblSynth.vecTheta, tblSynth.vecR] = cart2pol(tblSynth.dSngl, tblSynth.dBrst);
tblSynth.vecDeg = rad2deg(tblSynth.vecTheta);

% 3. Calculate Strategy Index (D)
% -------------------------------
tblSynth.D = tblSynth.dBrst - tblSynth.dSngl;

%% ========================================================================
%  PLOTTING: VECTOR STRATEGY
%  ========================================================================

fig1 = figure('Position', [100 100 1000 500], 'Color', 'w', 'Name', 'Vivo Recovery Vectors');
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Color Palette
cMap = containers.Map;
cMap('Control') = [0.2 0.2 0.2]; % Dark Gray
cMap('MCU-KO')  = [0.8 0.0 0.0]; % Red

% Axis Limits
maxVal = max(abs([tblSynth.dSngl; tblSynth.dBrst])) * 1.1;

for iGrp = 1:length(grps)
    nexttile;
    hold on; axis square; grid on;

    grp = grps{iGrp};
    subTbl = tblSynth(tblSynth.Group == grp, :);
    col = cMap(grp);

    % Orthogonal Regression (PCA)
    dataXY = [subTbl.dSngl, subTbl.dBrst];
    [coeff, ~, ~] = pca(dataXY);
    v1 = coeff(:, 1);
    mu = mean(dataXY);

    % Slope & Intercept
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

    % Scatter Points (Size by Baseline FR, Color by Baseline Burstiness)
    % Scaling for size
    scatSz = (log10(subTbl.fr_BSL) + 2) * 20;
    scatSz(scatSz < 5) = 5;

    scatter(subTbl.dSngl, subTbl.dBrst, scatSz, log(subTbl.pBspk_BSL), ...
        'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', 'k');

    % Text Stats (Slope)
    regAngle = atan2d(v1(2), v1(1));
    text(-maxVal*0.9, maxVal*0.9, sprintf('Slope: %.2f (%.1f\\circ)', slope, regAngle), ...
        'FontSize', 10, 'FontWeight', 'bold');

    % Formatting
    xlim([-maxVal maxVal]);
    ylim([-maxVal maxVal]);
    xlabel('Log Fold Change (Single)');
    ylabel('Log Fold Change (Burst)');
    title(sprintf('%s (Synthetic n=%d)', grp, height(subTbl)));

    colormap(gca, turbo);
    clim([-3 -0.5]); % Adjust based on pBspk log range
    if iGrp == 2
        c = colorbar;
        c.Label.String = 'log(Baseline pBspk)';
    end
    set(gca, 'FontSize', 12);
end

sgtitle('In Vivo Recovery Strategy (Pseudo-Tracking)');

%% ========================================================================
%  STATISTICS: STRATEGY INDEX (D)
%  ========================================================================

fprintf('\n================================================================\n');
fprintf(' STATISTICS: STRATEGY INDEX (D = dBrst - dSngl)\n');
fprintf('================================================================\n');

% Simple ANOVA on Synthetic Units (Since n is small, 10 per group)
% D ~ Group
[pVal, tblAnova, stats] = anova1(tblSynth.D, tblSynth.Group, 'off');

fprintf('ANOVA (D ~ Group):\n');
fprintf('F(%d, %d) = %.2f, p = %.4f\n', ...
    tblAnova{2,2}, tblAnova{2,3}, tblAnova{2,5}, tblAnova{2,6});

% Means
meanCtl = mean(tblSynth.D(tblSynth.Group == 'Control'));
meanMcu = mean(tblSynth.D(tblSynth.Group == 'MCU-KO'));
fprintf('Mean D (Control): %.3f\n', meanCtl);
fprintf('Mean D (MCU-KO) : %.3f\n', meanMcu);

% Plot Strategy Index
fig2 = figure('Position', [600 100 400 500], 'Color', 'w', 'Name', 'Strategy Index');
boxchart(tblSynth.Group, tblSynth.D, 'GroupByColor', tblSynth.Group);
hold on;
yline(0, 'k--', 'Balanced', 'LineWidth', 1.5);
ylabel('Strategy Index (D)');
title('Recovery Strategy Bias');
set(gca, 'FontSize', 12);
colororder([cMap('Control'); cMap('MCU-KO')]);

%% ========================================================================
%  STATISTICS: INTERACTION (dBrst ~ dSngl * Group)
%  ========================================================================

% fprintf('\n================================================================\n');
% fprintf(' STATISTICS: INTERACTION (dBrst ~ dSngl * Group)\n');
% fprintf('================================================================\n');
% 
% % Fit Linear Model
% lm = fitlm(tblSynth, 'dBrst ~ dSngl * Group');
% disp(lm);
% 
% % Extract Interaction p-value
% pInt = lm.Coefficients.pValue('dSngl:Group_MCU-KO');
% fprintf('Interaction p-value: %.4f\n', pInt);




