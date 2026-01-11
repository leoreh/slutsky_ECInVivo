% test_tbl_transform.m
% Comprehensive test suite for io/tbl_transform.m
% Verifies: Fit Mode, Apply Mode, Inverse Mode, Logic Flows, and Edge Cases.

root = 'd:\Code\slutsky_ECInVivo';
addpath(fullfile(root, 'io'));

fprintf('===================================================\n');
fprintf('TEST: tbl_transform.m (Comprehensive)\n');
fprintf('===================================================\n');

nFail = 0;

%% 1. Scenario: Basic Log + Z-Score (Fit -> Apply -> Inverse)
fprintf('\n--- Scenario 1: Basic Log + Z-Score (Fit -> Apply -> Inverse) ---\n');
try
    % Generate Data (Lognormal)
    rng(1);
    dataTrn = exp(randn(100, 1) * 0.5 + 2); % Mean ~ exp(2.125)
    dataTst = exp(randn(50, 1) * 0.5 + 2);

    T_Trn = table(dataTrn, 'VariableNames', {'X'});
    T_Tst = table(dataTst, 'VariableNames', {'X'});

    % A. FIT
    fprintf('Step A: Fitting Training Data...\n');
    % Force Log by setting skewThr to -Inf (or 0)
    [T_Trn_Out, params] = tbl_transform(T_Trn, 'logBase', 'e', 'flgZ', true, 'skewThr', 0, 'verbose', false);

    % Check Z-scoring on Logged Data
    logData = log(T_Trn.X);
    mu = mean(logData);
    sigma = std(logData);

    expectedZ = (logData - mu) / sigma;
    diff = max(abs(T_Trn_Out.X - expectedZ));

    if diff < 1e-10
        fprintf('[PASS] Fit: Output matches manual Z-score of Log data.\n');
    else
        fprintf('[FAIL] Fit: Output mismatch (Diff = %.4e)\n', diff); nFail=nFail+1;
    end

    % B. APPLY
    fprintf('Step B: Applying to Test Data...\n');
    [T_Tst_Out, ~] = tbl_transform(T_Tst, 'template', params);

    % Should use TRAINING mu/sigma
    logDataTst = log(T_Tst.X);
    expectedZTst = (logDataTst - mu) / sigma; % Uses TRN mu/sigma

    diffTst = max(abs(T_Tst_Out.X - expectedZTst));
    if diffTst < 1e-10
        fprintf('[PASS] Apply: Test data transformed using Training stats.\n');
    else
        fprintf('[FAIL] Apply: Output mismatch (Diff = %.4e)\n', diffTst); nFail=nFail+1;
    end

    % C. INVERSE (Test Data)
    fprintf('Step C: Inverse Transform (Recovery)...\n');
    T_Tst_Rec = tbl_transform(T_Tst_Out, 'template', params, 'flgInv', true);

    recDiff = max(abs(T_Tst.X - T_Tst_Rec.X));
    if recDiff < 1e-10
        fprintf('[PASS] Inverse: Recovered original data (Diff = %.4e).\n', recDiff);
    else
        fprintf('[FAIL] Inverse: Recovery failed (Diff = %.4e).\n', recDiff); nFail=nFail+1;
    end

catch ME
    fprintf('[ERROR] Scenario 1 Crashed: %s\n', ME.message);
    disp(ME.stack(1));
    nFail=nFail+1;
end


%% 2. Scenario: Logit Transform (Probabilities)
fprintf('\n--- Scenario 2: Logit Transform ---\n');
try
    % Data: Probabilities [0, 1]
    % Include edge cases close to 0/1 to test clamping
    X = [0.01; 0.5; 0.99; 0.000001; 0.999999];
    T = table(X);

    [T_Out, params] = tbl_transform(T, 'logBase', 'logit', 'flgZ', false);

    % Check logic
    % logit(p) = log(p/(1-p))
    expected = log(X ./ (1 - X));
    diff = max(abs(T_Out.X - expected));

    if diff < 1e-10
        fprintf('[PASS] Logit: Transformation correct.\n');
    else
        fprintf('[FAIL] Logit: Mismatch (Diff = %.4e)\n', diff); nFail=nFail+1;
    end

    % Inverse
    T_Rec = tbl_transform(T_Out, 'template', params, 'flgInv', true);
    diffInv = max(abs(T.X - T_Rec.X));

    if diffInv < 1e-10
        fprintf('[PASS] Logit Inverse: Correct recovery.\n', diffInv);
    else
        fprintf('[FAIL] Logit Inverse: Recovery failed (Diff = %.4e)\n', diffInv); nFail=nFail+1;
    end

catch ME
    fprintf('[ERROR] Scenario 2 Crashed: %s\n', ME.message); nFail=nFail+1;
end

%% 3. Scenario: Offset Handling (Log(0))
fprintf('\n--- Scenario 3: Offset Handling (Log(0)) ---\n');
try
    X = [0; 2; 10];
    T = table(X);

    % Auto-detect skew? No, force log.
    % Should detect 0 and add offset.
    [T_Out, params] = tbl_transform(T, 'logBase', 'e', 'flgZ', false, 'skewThr', 0);

    c = params.varList.X.offset;
    if c > 0
        fprintf('[PASS] Offset: Detected 0, added offset %.4f.\n', c);
    else
        fprintf('[FAIL] Offset: Did NOT add offset to 0-containing data.\n'); nFail=nFail+1;
    end

    % Check Values: log(X + c)
    expected = log(X + c);
    diff = max(abs(T_Out.X - expected));
    if diff < 1e-10
        fprintf('[PASS] Offset+Log: Values match.\n');
    else
        fprintf('[FAIL] Offset+Log: Mismatch.\n'); nFail=nFail+1;
    end

    % Inverse
    T_Rec = tbl_transform(T_Out, 'template', params, 'flgInv', true);
    diffInv = max(abs(T.X - T_Rec.X));
    if diffInv < 1e-10
        fprintf('[PASS] Offset Inverse: Correct recovery.\n');
    else
        fprintf('[FAIL] Offset Inverse: Recovery failed (Diff = %.4e).\n', diffInv); nFail=nFail+1;
    end

catch ME
    fprintf('[ERROR] Scenario 3 Crashed: %s\n', ME.message); nFail=nFail+1;
end

%% 4. Scenario: Grouped Normalization
fprintf('\n--- Scenario 4: Grouped Normalization ---\n');
try
    Group = {'A'; 'A'; 'A'; 'B'; 'B'; 'B'};
    Val   = [10; 10; 10; 20; 20; 20];
    T = table(Group, Val);
    T.Group = categorical(T.Group);

    % 1. Error Check: varNorm included in varsGrp
    try
        tbl_transform(T, 'varsGrp', 'Group', 'varNorm', 'Group');
        fprintf('[FAIL] Did not catch invalid varsGrp/varNorm overlap.\n'); nFail=nFail+1;
    catch
        fprintf('[PASS] Correctly errored on overlapping varsGrp/varNorm (Conflict).\n');
    end

    % 2. Correct Usage: Global Norm relative to Group A
    [T_Norm2, pNorm2] = tbl_transform(T, 'varsGrp', [], 'varNorm', 'Group', 'flgZ', false);

    if all(T_Norm2.Val(1:3) == 100) && all(T_Norm2.Val(4:6) == 200)
        fprintf('[PASS] Normalization (Global) works as expected.\n');
    else
        fprintf('[FAIL] Normalization failed. A=%.1f (Exp 100), B=%.1f (Exp 200).\n', ...
            mean(T_Norm2.Val(1:3)), mean(T_Norm2.Val(4:6)));
        nFail=nFail+1;
    end

    % INVERSE Norm
    T_Rec2 = tbl_transform(T_Norm2, 'template', pNorm2, 'flgInv', true);
    if max(abs(T_Rec2.Val - T.Val)) < 1e-10
        fprintf('[PASS] Normalization Inverse works.\n');
    else
        fprintf('[FAIL] Norm Inverse failed.\n'); nFail=nFail+1;
    end

catch ME
    fprintf('[ERROR] Scenario 4 Crashed: %s\n', ME.message); nFail=nFail+1;
end

%% 5. Scenario: Norm + Z-Score (Conflict Resolution)
fprintf('\n--- Scenario 5: Norm + Z-Score (Conflict Resolution) ---\n');
try
    % Data: [10, 20]. Mean 15, SD 7.07.
    % Requested: Norm + Z.
    % Expectation: Norm disabled. Z applied to raw data.

    Val = [10; 20];
    T = table(Val);

    [T_Out, params] = tbl_transform(T, 'varNorm', 'Dummy', 'flgZ', true, 'verbose', false);

    % Check Norm is disabled
    if isfield(params.varList.Val, 'flgNorm') && ~params.varList.Val.flgNorm
        fprintf('[PASS] Conflict: Norm was correctly disabled when Z is enabled.\n');
    else
        fprintf('[FAIL] Conflict: Norm was NOT disabled.\n'); nFail=nFail+1;
    end

    % Check Z values
    mu = 15; sigma = std([10; 20]); % 7.0711 creates z=[-0.7071, 0.7071]
    expected = ([10;20] - mu) / sigma;

    diff = max(abs(T_Out.Val - expected));
    if diff < 1e-10
        fprintf('[PASS] Z-Score applied correctly on raw data.\n');
    else
        fprintf('[FAIL] Z-Score mismatch.\n'); nFail=nFail+1;
    end

catch ME
    fprintf('[ERROR] Scenario 5 Crashed: %s\n', ME.message); nFail=nFail+1;
end

%% Final Report
fprintf('\n==================\n');
if nFail == 0
    fprintf('ALL TESTS PASSED.\n');
else
    fprintf('%d TESTS FAILED.\n', nFail);
end
fprintf('==================\n');

