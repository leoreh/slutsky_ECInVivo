function tests = test_vif
tests = functiontests(localfunctions);
end

function testLmeVif_Direct(testCase)
% Create correlated data
rng(42);
x1 = randn(100, 1);
x2 = x1 * 0.9 + randn(100, 1) * 0.1; % Highly correlated
x3 = randn(100, 1); % Independent

tbl = table(x1, x2, x3);

% Call lme_vif
vifTbl = lme_vif(tbl, 'y ~ x1 + x2 + x3');

% Check Structure
verifyEqual(testCase, vifTbl.Properties.VariableNames, {'Predictor', 'VIF'});
verifyEqual(testCase, height(vifTbl), 3);

% Check Values
% VIF(x1) should be high
rowX1 = vifTbl(strcmp(vifTbl.Predictor, 'x1'), :);
verifyTrue(testCase, rowX1.VIF > 5, 'VIF for x1 should be high');

% VIF(x3) should be near 1
rowX3 = vifTbl(strcmp(vifTbl.Predictor, 'x3'), :);
verifyTrue(testCase, rowX3.VIF < 2, 'VIF for x3 should be low');
end

function testLmeAnalyse_Integration(testCase)
% Create data with adequate predictors (>=2 numeric)
tbl = table(randn(20,1), randn(20,1), randn(20,1), 'VariableNames', {'y', 'x1', 'x2'});
% Create simple lme formula (no random effects for simplicity of test, but lme_analyse might expect one?)
% lme_analyse calls lme_fit. lme_fit usually handles GLME/LME. If I use naive formula y~x1+x2, it might fail if lme_fit expects RE.
% Let's assume lme_analyse works for simple regression if no RE, or just add a dummy group.
tbl.Name = repmat({'A'; 'B'}, 10, 1);
frml = 'y ~ x1 + x2 + (1|Name)';

% Run analyse
[~, ~, lmeInfo] = lme_analyse(tbl, frml, 'verbose', false);

% Check VIF in lmeInfo
verifyTrue(testCase, isfield(lmeInfo, 'vif'), 'lmeInfo should have vif field');
verifyFalse(testCase, isempty(lmeInfo.vif), 'VIF table should not be empty');
verifyEqual(testCase, height(lmeInfo.vif), 2);
end

function testLmeAbation_Integration(testCase)
% Create data
tbl = table(randn(20,1), randn(20,1), randn(20,1), 'VariableNames', {'y', 'x1', 'x2'});
tbl.Name = repmat({'A'; 'B'}, 10, 1);
frml = 'y ~ x1 + x2 + (1|Name)';

% Run ablation (minimal reps for speed)
res = lme_ablation(tbl, frml, 'nReps', 1, 'flgPlot', false);

% Check VIF in res
verifyTrue(testCase, isfield(res, 'vif'), 'res should have vif field');
verifyFalse(testCase, isempty(res.vif), 'VIF table should not be empty');
end

function testSinglePredictor(testCase)
% Only 1 numeric predictor -> Should return empty VIF
tbl = table(randn(20,1), randn(20,1), 'VariableNames', {'y', 'x1'});
tbl.Name = repmat({'A'; 'B'}, 10, 1);
frml = 'y ~ x1 + (1|Name)';

[~, ~, lmeInfo] = lme_analyse(tbl, frml, 'verbose', false);

verifyTrue(testCase, isfield(lmeInfo, 'vif'));
verifyTrue(testCase, isempty(lmeInfo.vif), 'VIF should be empty for single predictor');
end
