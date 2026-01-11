function test_tbl_trans_filtering()

% Add IO path
rootDir = fileparts(fileparts(mfilename('fullpath'))); % Go up one level from tests/
addpath(fullfile(rootDir, 'io'));

fprintf('Running varsInc/varsExc Filtering Tests (tbl_trans)...\n');

% Setup Data
T = table();
T.A = rand(10,1);
T.B = rand(10,1) * 10; % Bigger scale
T.C = rand(10,1);

% 1. Fit Mode (Test varsInc/Exc in Fit)
fprintf('Test 1: Fit Mode Filtering... ');
[T_Fit, P_Fit] = tbl_trans(T, 'flgZ', true, 'varsInc', {'A', 'B'});

assert(isfield(P_Fit.varsTrans, 'A'), 'A should be transformed');
assert(isfield(P_Fit.varsTrans, 'B'), 'B should be transformed');
assert(~isfield(P_Fit.varsTrans, 'C'), 'C should NOT be transformed');

% Check Z-scoring happened
assert(abs(mean(T_Fit.A)) < 1e-10, 'A should be centered');
assert(abs(mean(T_Fit.B)) < 1e-10, 'B should be centered');
assert(abs(mean(T_Fit.C) - mean(T.C)) < 1e-10, 'C should be untouched');

fprintf('PASS\n');

% 2. Inverse Mode with varsInc
% We want to reverse only A, leaving B transformed.
fprintf('Test 2: Inverse Mode with varsInc... ');

[T_Inv, ~] = tbl_trans(T_Fit, 'template', P_Fit, 'flgInv', true, 'varsInc', {'A'});

% A should be back to original
assert(max(abs(T_Inv.A - T.A)) < 1e-10, 'A should be recovered');

% B should still be transformed (Z-scored)
assert(max(abs(T_Inv.B - T_Fit.B)) < 1e-10, 'B should remain transformed');

fprintf('PASS\n');

% 3. Inverse Mode with varsExc
% Exclude A, so only B should be reversed.
fprintf('Test 3: Inverse Mode with varsExc... ');

[T_Inv2, ~] = tbl_trans(T_Fit, 'template', P_Fit, 'flgInv', true, 'varsExc', {'A'});

% A should match T_Fit (transformed)
assert(max(abs(T_Inv2.A - T_Fit.A)) < 1e-10, 'A should remain transformed');

% B should match T (recovered)
assert(max(abs(T_Inv2.B - T.B)) < 1e-10, 'B should be recovered');

fprintf('PASS\n');

% 4. Apply Mode with varsInc
% Apply transform only to A
fprintf('Test 4: Apply Mode with varsInc... ');

T_New = T;
% Make B different to ensure we notice if it gets transformed
T_New.B = T.B + 100;

[T_App, ~] = tbl_trans(T_New, 'template', P_Fit, 'varsInc', {'A'});

% A should be transformed (using stats from T)
muA = P_Fit.varsTrans.A.stats.Mean(1);
sigA = P_Fit.varsTrans.A.stats.SD(1);
expectedA = (T.A - muA) / sigA;
assert(max(abs(T_App.A - expectedA)) < 1e-10, 'A should be transformed');

% B should be UNTOUCHED
assert(max(abs(T_App.B - T_New.B)) < 1e-10, 'B should be untouched');

fprintf('PASS\n');

fprintf('ALL FILTERING TESTS PASSED.\n');
end
