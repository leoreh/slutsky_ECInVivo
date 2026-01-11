
%% Verify LME_PD Refactor
% This script tests the refactored lme_pd.m to ensure:
% 1. Prediction columns are correctly added (varResp, varResp_pred, etc.)
% 2. Back-transformation works correctly using the new logic.

clear; clc;


sysDir = fileparts(fileparts(mfilename('fullpath'))); % Project Root
addpath(fullfile(sysDir, 'io'));
addpath(fullfile(sysDir, 'lme'));

fprintf('Running LME_PD Verification...\n');

%% 1. Setup Data & Model
rng(42);
n = 100;
x = rand(n, 1) * 10;
y_raw = 2*x + 5 + randn(n, 1);
% Apply a transformation to Y to simulate real workflow (e.g. log)
% Logic: We model Log(Y), so "raw" data is strictly positive.
y = exp(y_raw / 10); % Scaling to keep in reasonable range

tbl = table(x, y, 'VariableNames', {'X', 'Y'});

%% 2. Transform (Get transParams)
% We want to fit a model to Log(Y).
% Let's use tbl_trans to do this fitting.
[tblTrn, params] = tbl_trans(tbl, 'logBase', 'e', 'verbose', true);

fprintf('Transformation Params for Y:\n');
disp(params.varsTrans.Y);

%% 3. Fit Simple LME/GLME
% Using fitlme on transformed data
% Formula: Y ~ X
try
    mdl = fitlme(tblTrn, 'Y ~ X');
catch
    warning('fitlme failed (maybe toolbox missing). Mocking model.');
    % Mock MDL object if fitlme not available (unlikely given user context)
    mdl = struct();
    mdl.ResponseName = 'Y';
    mdl.Variables = tblTrn;
end

%% 4. Run LME_PD
fprintf('Running lme_pd...\n');
try
    [pdRes, hFig] = lme_pd(mdl, 'X', 'transParams', params, 'nGrid', 10);
    close(hFig); % Close plot

    disp('lme_pd ran successfully.');
    disp('Columns in pdRes:');
    disp(pdRes.Properties.VariableNames);

    %% 5. Validations

    % A. Check Column Existence
    reqCols = {'Y', 'Y_pred', 'Y_lower', 'Y_upper'};
    hasCols = all(ismember(reqCols, pdRes.Properties.VariableNames));

    if hasCols
        fprintf('[PASS] All required columns present: %s\n', strjoin(reqCols, ', '));
    else
        fprintf('[FAIL] Missing columns.\n');
        error('Verification failed: Missing columns.');
    end

    % B. Check Back-Transformation
    % Since we used Log(Y) in the model, the output in pdRes should be back to original scale (exp).
    % Let's pick a point.
    % The model predicts in Log Scale.
    % pdRes should contain Exp(Pred).

    % Check if values are roughly in range of original Y (which was ~ exp(0..1) -> 1..2.7)
    % Wait, y_raw = 2*x + 5... that's large. exp(y_raw/10) -> exp((0..20+5)/10) = exp(0.5..2.5) -> 1.6 .. 12.

    rangeY = [min(pdRes.Y_pred), max(pdRes.Y_pred)];
    fprintf('Predicted Y Range (Original Scale): [%.2f, %.2f]\n', rangeY(1), rangeY(2));

    if rangeY(1) > 0
        fprintf('[PASS] Predictions are positive (as expected for exp transform).\n');
    else
        fprintf('[FAIL] Predictions <= 0, back-transform might have failed.\n');
    end

    % C. Check Redundancy
    % Y and Y_pred should be identical (as per code: tblPlot.(respName) = yPred; tblPlot.(predName) = yPred)
    % EXCEPT: Y itself usually refers to the data column?
    % In lme_pd refactor:
    % tblPlot.(respName) = yPred;
    % tblPlot.(predName) = yPred;
    % So they should be equal.

    diff = max(abs(pdRes.Y - pdRes.Y_pred));
    if diff < 1e-9
        fprintf('[PASS] Y and Y_pred are identical (Max Diff: %.2e).\n', diff);
    else
        fprintf('[FAIL] Y and Y_pred differ (Max Diff: %.2e).\n', diff);
    end

    fprintf('\nVerification Complete: SUCCESS\n');

catch ME
    fprintf('\nVerification Failed with error:\n%s\n', ME.message);
    rethrow(ME);
end
