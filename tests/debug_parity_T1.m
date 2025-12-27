
% Verify parity logic for Test 1
nRows = 20;
nCols = 100;
tbl = table();
tbl.Group = string(repmat({'A'; 'B'; 'C'; 'D'}, nRows/4, 1));
tbl.Mouse = string(repmat({'M1'; 'M2'}, nRows/2, 1));
tbl.LFP = randn(nRows, nCols);

% Original
content = fileread('d:\Code\slutsky_ECInVivo\io\tbl_tNorm_original.m');
content = strrep(content, 'function tblOut = tbl_tNorm(tbl, varargin)', 'function tblOut = tbl_tNorm_original(tbl, varargin)');
fid = fopen('d:\Code\slutsky_ECInVivo\io\tbl_tNorm_original.m', 'w');
fwrite(fid, content);
fclose(fid);

T1 = tbl_tNorm_original(tbl, 'varsInc', {'LFP'});
T2 = tbl_tNorm(tbl, 'varsInc', {'LFP'});

% Compare Means
diffLFP = T1.LFP - T2.LFP;
maxDiff = max(abs(diffLFP(:)));
disp(['Max Diff LFP: ' num2str(maxDiff)]);

if maxDiff > 1e-10
    disp('Significant difference found.');
    % Inspect first row
    row1_orig = T1.LFP(1, 1:5);
    row1_new = T2.LFP(1, 1:5);
    disp('Original Row 1 (first 5):');
    disp(row1_orig);
    disp('New Row 1 (first 5):');
    disp(row1_new);
else
    disp('Difference is negligible.');
end
