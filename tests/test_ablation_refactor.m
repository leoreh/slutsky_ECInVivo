% Test LME Ablation Refactor

disp(['Current Dir: ' pwd]);
root = 'd:\Code\slutsky_ECInVivo';
addpath(fullfile(root, 'lme'));
addpath(fullfile(root, 'io'));
% addpath('tests');

% 1. Create Dummy Data
n = 40;
tbl = table();
tbl.Name = repmat({'M1';'M2';'M3';'M4'}, n/4, 1);
tbl.Group = repmat({'Ctrl';'Ctrl';'KO';'KO'}, n/4, 1);
tbl.Fr = exp(randn(n, 1)+2); % Log-normal-ish
tbl.Burst = rand(n, 1);
tbl.uRcv = exp(randn(n, 1)); % Response (Continuous for Log-Normal test)

% Introduces correlation for ablation
tbl.uRcv(tbl.Fr > 10) = 1;

frml = 'uRcv ~ Fr + Burst + (1|Name)';

% 2. Run Ablation (Normal / Binomial)
fprintf('Running Ablation Check (Log-Normal)...\n');
try
    % Use Log-Normal to test Response Transformation logic
    res = lme_ablation(tbl, frml, 'dist', 'Log-Normal', 'nReps', 2, 'flgPlot', false);
    fprintf('SUCCESS: Ablation finished.\n');
    disp(res.vars);
    disp('Primary Perf (RMSE):');
    disp(res.perfPrim);
    disp(res.impPrim);
catch ME
    fprintf('FAILURE: %s\n', ME.message);
    disp(ME.stack(1));
end
