function test_lme_ablation_metrics()
% TEST_LME_ABLATION_METRICS Verifies correct calculation of Pooled RMSE and R2

rng(42); % Reproducibility

% 1. Create Synthetic Data (Regression)
n = 100;
x = randn(n, 1);
group = repmat({'A'; 'B'}, n/2, 1); % Ensure column cell array
name = string((1:n)'); % Ensure column string array
y = 2*x + 5 + 0.5*randn(n, 1); % Linear relation

tbl = table(x, group, name, y);
frml = 'y ~ x + (1|group)';

% 2. Run lme_ablation
% Use small folds for speed
fprintf('Running lme_ablation...\n');
res = lme_ablation(tbl, frml, 'nReps', 2, 'flgPlot', false);

% 3. Verify Pooled Metrics
% Re-calculate manually from the stored sums in res
sumSSE = nansum(res.sse, 1);
sumSST = nansum(res.sst, 1);
totalN = sum(res.foldSz, 'omitnan');

calcRMSE = sqrt(sumSSE ./ totalN);
calcR2   = 1 - (sumSSE ./ sumSST);

% Check RMSE
assert(all(abs(res.pRMSE - calcRMSE) < 1e-10), 'Pooled RMSE mismatch');
fprintf('Pooled RMSE check passed.\n');

% Check R2
assert(all(abs(res.pR2 - calcR2) < 1e-10), 'Pooled R2 mismatch');
fprintf('Pooled R2 check passed.\n');

% 4. Verify Importance Logic
% Prim: (Reduced - Full) / Full * 100 (Percentage Increase)
impPrim = (res.pRMSE(2:end) - res.pRMSE(1)) ./ res.pRMSE(1) * 100;
assert(all(abs(res.dRMSE - impPrim) < 1e-10), 'Importance Primary mismatch');
fprintf('Importance Primary (RMSE) check passed.\n');

% Sec: Full - Reduced (Difference)
impSec = res.pR2(1) - res.pR2(2:end);
assert(all(abs(res.dR2 - impSec) < 1e-10), 'Importance Secondary mismatch');
fprintf('Importance Secondary (R2) check passed.\n');

fprintf('\nALL TESTS PASSED.\n');

end
