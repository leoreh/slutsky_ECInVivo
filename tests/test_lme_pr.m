%% Test lme_pr
% Verification script for lme_pr function

clear; close all; clc;
addpath(genpath('d:\Code\slutsky_ECInVivo'));

% 1. create synthetic data
rng(42);
n = 100;

% Make x1 Z-scored (simulate transform)
rawX1 = randn(n, 1) * 10 + 50; % Mean 50, Std 10
mu = mean(rawX1);
sigma = std(rawX1);
x1 = (rawX1 - mu) / sigma;

x2 = randn(n, 1);

% Fix grp creation: Ensure it is 100x1
grp = [repmat({'Control'}, n/2, 1); repmat({'Treated'}, n/2, 1)];
grp = grp(randperm(n));

subj = repmat((1:10)', 10, 1);
subj = subj(1:n);

% True model: y = 2*x1 + 3*x2 + (grp=='Treated')*1 + noise
noise = randn(n, 1) * 0.5;
y = 2*x1 + 3*x2 + strcmp(grp, 'Treated') + noise;

tbl = table(y, x1, x2, grp, subj);

% Create Fake transParams
transParams = struct();
transParams.varsTrans.x1.flgZ = true;
transParams.varsTrans.x1.flgNorm = false;
transParams.varsTrans.x1.stats.Mean = mu;
transParams.varsTrans.x1.stats.SD = sigma;
transParams.varsTrans.x1.logBase = [];
transParams.varsTrans.x1.offset = 0;
transParams.varsGrp = {};

% 2. Fit Full Model
% Interaction: grp * x1
frml = 'y ~ x1 * grp + x2 + (1|subj)';
try
    mdl = fitlme(tbl, frml);
    fprintf('Model fitted successfully.\n');
catch me
    warning('Failed to fit model: %s', me.message);
    return;
end

%% Test 1: 1 Variable (Continuous) without Back-Transform
% Plot should show Z-scored x1 (around 0)
fprintf('\n--- Test 1: lme_pr(mdl, {''x1''}) [Z-Scored] ---\n');
try
    [hFig1, tblRes1] = lme_pr(mdl, 'x1', 'verbose', true);
    title(gca, 'Test 1: Z-Scored x1');
    fprintf('Test 1 passed.\n');
    close(hFig1);
catch me
    fprintf('Test 1 FAILED: %s\n', me.message);
    disp(getReport(me));
end

%% Test 2: 1 Variable (Continuous) WITH Back-Transform
% Plot should show Original x1 (around 50) and straight line (since raw relationship is linear)
fprintf('\n--- Test 2: lme_pr(mdl, {''x1''}) [Back-Transformed] ---\n');
try
    [hFig2, tblRes2] = lme_pr(mdl, 'x1', 'transParams', transParams, 'verbose', true);
    title(gca, 'Test 2: Original x1 (Mean ~50)');
    fprintf('Test 2 passed.\n');

    % Verify back-transformation
    % tblRes2 is returned from lme_pr.
    % Wait! lme_pr (new version) does NOT return the back-transformed table in output argument
    % unless we change the output arg.
    % The function returns [hFig, tblRes].
    % In the new code: tblPlot = tbl_trans(...) but tblRes is NOT updated except for fitting loop.
    % Actually, let's check lme_pr output.
    % "tblFit = tblRes; ... tblPlot = tbl_trans(...)".
    % The function outputs `tblRes` (the one passed in + Resid).
    % It does NOT output `tblPlot`.
    % So tblRes2 should still be Z-scored!

    if abs(mean(tblRes2.x1)) < 0.5
        fprintf('   -> Output Table is Transformed (correct logic, Mean=%.1f)\n', mean(tblRes2.x1));
    else
        fprintf('   -> Output Table is Back-Transformed (unexpected? Mean=%.1f)\n', mean(tblRes2.x1));
    end

    close(hFig2);

catch me
    fprintf('Test 2 FAILED: %s\n', me.message);
    disp(getReport(me));
end

fprintf('\nVerification Complete.\n');
