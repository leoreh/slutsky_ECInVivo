%% ========================================================================
%  STATE SPACE PCA
%  ========================================================================
%  Hypothesis: Deviation from the conserved "Recovery Rule" (PC1) is 
%  predicted by Baseline Burstiness.
%
%  1. Calculate Allocation Residuals (controlling for BSL activity)
%  2. Perform PCA on the Residuals (PC1 = Rule, PC2 = Deviation)
%  3. Correlate PC2 with Baseline Burstiness

fprintf('\n================================================================\n');
fprintf(' STATE SPACE PCA\n');
fprintf('================================================================\n');

% --- 1. Residuals (Allocation) ---
% Regress out Baseline FR and Burstiness to isolate intrinsic plasticity
resType  = 'Raw';
varBrst = 'pBspk_trans';
regType = 'linear';

% Burst
frml = 'dBrst_rel ~ (fr + dSngl_rel) * Group + (1|Name)';
mdl = lme_analyse(tbl, frml, 'dist', 'normal', 'verbose', false);
tbl.R_dBrst_rel = residuals(mdl, 'ResidualType', resType);

% Single
frml = 'dSngl_rel ~ (fr + dBrst_rel) * Group + (1|Name)';
mdl = lme_analyse(tbl, ['dSngl_rel' frmlBase], 'dist', 'normal', 'verbose', false);
tbl.R_dSngl_rel = residuals(mdl, 'ResidualType', resType);


% --- 2. PCA ---
varsPCA = {'dBrst_rel', 'dSngl_rel'};
X = table2array(tbl(:, varsPCA));

% PCA (Rows=Observations, Cols=Variables)
[coeff, score, latent, tsquared, explained] = pca(X);

% Store Scores (Rotate the data)
tbl.PC1 = score(:,1); 
tbl.PC2 = score(:,2);

fprintf('PCA Variance Explained: PC1=%.1f%%, PC2=%.1f%%\n', explained(1), explained(2));
disp(array2table(coeff, 'RowNames', varsPCA, 'VariableNames', {'PC1', 'PC2'}));


% --- 3. Statistical Test (PC2 vs Burstiness) ---
fprintf('\n--- Predicting Deviation (PC2) ---\n');
frml = sprintf('PC2 ~ %s * Group + (1|Name)', varBrst);
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'verbose', true);


% --- 4. Visualization ---
hFig = figure('Color', 'w', 'Name', 'State Space PCA');
tlo = tiledlayout(hFig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
cfg = mcu_cfg;
clr = cfg.clr.grp;
groups = unique(tbl.Group);


% Panel A: Residual Space + PCA Vectors
% -------------------------------------
nexttile; hold on;
title('Residual State Space');

xVar = varsPCA{1};
yVar = varsPCA{2};
lims = max(abs([tbl.(xVar); tbl.(yVar)])) * 1.1;

% Scatter
for i = 1:length(groups)
    idx = tbl.Group == groups(i);
    scatter(tbl.(xVar)(idx), tbl.(yVar)(idx), 15, clr(i,:), 'filled', ...
        'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor', 'none');
end

% PCA Vectors
scale = lims * 0.8;
quiver(0, 0, coeff(1,1)*scale, coeff(2,1)*scale, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.5, 'DisplayName', 'PC1');
quiver(0, 0, coeff(1,2)*scale, coeff(2,2)*scale, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5, 'DisplayName', 'PC2');

xlabel('Burst Residuals (\DeltaBrst)'); 
ylabel('Single Residuals (\DeltaSngl)');
xlim([-lims lims]); ylim([-lims lims]);
xline(0, '--k', 'HandleVisibility', 'off');
yline(0, '--k', 'HandleVisibility', 'off');
axis square;


% Panel B: Rotated Space (PC1 vs PC2)
% -----------------------------------
nexttile; hold on;
title('Rotated Space (PCA)');

for i = 1:length(groups)
    idx = tbl.Group == groups(i);
    scatter(tbl.PC1(idx), tbl.PC2(idx), 15, clr(i,:), 'filled', ...
        'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor', 'none');
end

xlabel('PC1 (Synergy/Rule)');
ylabel('PC2 (Deviation)');
xline(0, '--k'); yline(0, '--k');
axis square;


% Panel C: PC2 Prediction
% -----------------------------------
nexttile; hold on;
title('Predicting Deviation');

xVar = varBrst;
yVar = 'PC2';

for i = 1:length(groups)
    idx = tbl.Group == groups(i);
    x = tbl.(xVar)(idx);
    y = tbl.(yVar)(idx);
    c = clr(i,:);
    
    scatter(x, y, 15, c, 'filled', ...
        'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor', 'none');
        
    % Regression Line
    plot_linReg(x, y, 'hAx', gca, 'type', regType, 'clr', c, 'flgTxt', true);
end

xlabel('Baseline Burstiness (Logit)');
ylabel('PC2 Score');
axis square;
