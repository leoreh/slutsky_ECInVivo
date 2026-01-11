% Test LME Ablation Stats
% Verifies that lme_ablation returns res.stats with correct LRT results.

root = 'd:\Code\slutsky_ECInVivo';
addpath(fullfile(root, 'lme'));
addpath(fullfile(root, 'io'));

fprintf('---------------------------------------------------\n');
fprintf('TEST: LME Ablation Statistics (LRT, AIC, p-values)\n');
fprintf('---------------------------------------------------\n');

%% 1. Create Dummy Data
% Scenario: Y = 2*X1 + 0*X2 + RandomEffect
n = 100;
rng(42); % Reproducibility

nameList = {'A', 'B', 'C', 'D', 'E'};
tbl = table();
tbl.Name = repmat(nameList(:), 20, 1);


% Predictors
tbl.X1 = randn(n, 1);          % Significant
tbl.X2 = randn(n, 1);          % Noise

% Fix Group generation (ensure 100x1)
tbl.Group = repmat({'Ctrl'; 'Tx'}, 50, 1);

% Response (Normal)
% Fixed effect: 2*X1
% Random effect intercept per Name
% Simple manual mapping for random effects validation
uNames = nameList;
nEffVals = randn(length(uNames), 1) * 2;
mapEff = containers.Map(uNames, nEffVals);
nameEff = values(mapEff, tbl.Name); % Returns cell array
nameEff = cell2mat(nameEff'); % Convert to column vector

tbl.Y = 2 * tbl.X1 + nameEff + randn(n, 1) * 0.5;

frml = 'Y ~ X1 + X2 + (1|Name)';

%% 2. Run Ablation (Normal / ML)
fprintf('\n--- Running Normal (LME) Ablation ---\n');
try
    res = lme_ablation(tbl, frml, 'dist', 'Normal', 'nReps', 2, 'flgPlot', false);

    if ~isfield(res, 'stats')
        error('Result struct missing .stats field');
    end

    disp('Result Stats Table:');
    disp(res.stats);

    % Validation
    fullAIC = res.stats.AIC(strcmp(res.stats.Predictor, 'None'));

    % X1 (Significant) -> Reduced model should be worse -> Higher AIC than Full -> DeltaAIC > 0
    % X1 Removed:
    rowX1 = res.stats(strcmp(res.stats.Predictor, 'X1'), :);
    pValX1 = rowX1.pValue;
    dAICX1 = rowX1.DeltaAIC;

    if pValX1 < 0.05
        fprintf('[PASS] X1 is significant (p=%.4g)\n', pValX1);
    else
        fprintf('[FAIL] X1 should be significant! (p=%.4g)\n', pValX1);
    end

    if dAICX1 > 0
        fprintf('[PASS] Removing X1 increases AIC (Delta=%.4f)\n', dAICX1);
    else
        fprintf('[FAIL] Removing X1 should increase AIC! (Delta=%.4f)\n', dAICX1);
    end

    % X2 (Noise) -> Reduced model should be similar/better -> DeltaAIC small or negative
    rowX2 = res.stats(strcmp(res.stats.Predictor, 'X2'), :);
    pValX2 = rowX2.pValue;

    if pValX2 > 0.05
        fprintf('[PASS] X2 is not significant (p=%.4g)\n', pValX2);
    else
        fprintf('[WARN] X2 found significant (p=%.4g) - might be chance\n', pValX2);
    end

catch ME
    fprintf('FAILURE: %s\n', ME.message);
    disp(ME.stack(1));
end

%% 3. Run Ablation (Binomial / Laplace)
fprintf('\n--- Running Binomial (GLME) Ablation ---\n');
try
    % Create Binary Response correlated with X1
    prob = 1 ./ (1 + exp(-(2*tbl.X1))); % Sigmoid
    tbl.YBin = double(rand(n, 1) < prob);

    frmlBin = 'YBin ~ X1 + X2 + (1|Name)';

    resBin = lme_ablation(tbl, frmlBin, 'dist', 'Binomial', 'nReps', 2, 'flgPlot', false);

    disp('Result Stats Table (Binomial):');
    disp(resBin.stats);

    % Check p-value for X1
    pValBinX1 = resBin.stats.pValue(strcmp(resBin.stats.Predictor, 'X1'));
    if pValBinX1 < 0.05
        fprintf('[PASS] X1 is significant in Binomial (p=%.4g)\n', pValBinX1);
    else
        fprintf('[FAIL] X1 should be significant in Binomial! (p=%.4g)\n', pValBinX1);
    end

catch ME
    fprintf('FAILURE: %s\n', ME.message);
    disp(ME.stack(1));
end
