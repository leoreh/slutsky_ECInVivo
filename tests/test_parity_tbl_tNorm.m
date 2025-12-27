
% Test Parity between tbl_tNorm and tbl_tNorm_original

% Modify tbl_tNorm_original function name to run it
% (Assuming we just copy the file, the function name inside is still tbl_tNorm)
% We need to dynamically call it or rename it.
% Easier: Rename the file content too.

try
    content = fileread('d:\Code\slutsky_ECInVivo\io\tbl_tNorm_original.m');
    content = strrep(content, 'function tblOut = tbl_tNorm(tbl, varargin)', 'function tblOut = tbl_tNorm_original(tbl, varargin)');
    fid = fopen('d:\Code\slutsky_ECInVivo\io\tbl_tNorm_original.m', 'w');
    fwrite(fid, content);
    fclose(fid);
catch ME
    disp('Error renaming function in original file copy:');
    disp(ME.message);
end

% Create Test Data
nRows = 20;
nCols = 100;
tbl = table();
tbl.Group = string(repmat({'A'; 'B'; 'C'; 'D'}, nRows/4, 1));
tbl.Mouse = string(repmat({'M1'; 'M2'}, nRows/2, 1));
tbl.LFP = randn(nRows, nCols);
tbl.FR = rand(nRows, nCols) * 10;
tbl.Scalar = randn(nRows, 1); % Should be ignored or handled?

disp('--- Running Parity Tests ---');

% Test 1: Simple percentage (default)
disp('Test 1: Simple percentage (default)...');
T1 = tbl_tNorm_original(tbl, 'varsInc', {'LFP', 'FR'});
T2 = tbl_tNorm(tbl, 'varsInc', {'LFP', 'FR'});
assert(isequaln(T1.LFP, T2.LFP), 'Test 1 LFP Mismatch');
assert(isequaln(T1.FR, T2.FR), 'Test 1 FR Mismatch');

% Test 2: Group Z-Score with Window
disp('Test 2: Group Z-Score (Group, Mouse) with Window [10 50]...');
T1 = tbl_tNorm_original(tbl, 'varsInc', {'LFP'}, 'varsGrp', {'Group', 'Mouse'}, 'Method', 'zscore', 'winNorm', [10 50]);
T2 = tbl_tNorm(tbl, 'varsInc', {'LFP'}, 'varsGrp', {'Group', 'Mouse'}, 'Method', 'zscore', 'winNorm', [10 50]);
assert(isequaln(T1.LFP, T2.LFP), 'Test 2 LFP Mismatch');

% Test 3: Range Method
disp('Test 3: Range Method...');
T1 = tbl_tNorm_original(tbl, 'varsInc', {'FR'}, 'varsGrp', {'Group'}, 'Method', 'range');
T2 = tbl_tNorm(tbl, 'varsInc', {'FR'}, 'varsGrp', {'Group'}, 'Method', 'range');
assert(isequaln(T1.FR, T2.FR), 'Test 3 FR Mismatch');

% Test 4: Center Method
disp('Test 4: Center Method...');
T1 = tbl_tNorm_original(tbl, 'varsInc', {'FR'}, 'Method', 'center');
T2 = tbl_tNorm(tbl, 'varsInc', {'FR'}, 'Method', 'center');
assert(isequaln(T1.FR, T2.FR), 'Test 4 FR Mismatch');

disp('ALL PARITY TESTS PASSED!');
