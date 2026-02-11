
%% test_lme_mediation_units.m
% Validates unit consistency in lme_mediation.m
% specifically regarding distM='log-normal' and predictor transformations.

clear; clc;

% 1. Simulate Data
% Structure: X -> log(M) -> Y
rng(42);
N = 100;
Name = repmat({'A';'B';'C';'D';'E'}, N/5, 1);
X = randn(N, 1); % Independent Var

% M is Log-Normal
% log(M) ~ 0.5*X + noise
logM = 0.5 * X + 0.2 * randn(N, 1);
M = exp(logM); 

% Y is Normal
% Y ~ 2.0*log(M) + X + noise
% Direct effect of X (c') = 1.0
% Indirect effect path b = 2.0
% Total indirect = a*b = 0.5 * 2.0 = 1.0
Y = 1.0 * X + 2.0 * logM + 0.2 * randn(N, 1);

tbl = table(Name, X, M, Y);

% 2. Run lme_mediation with distM='log-normal'
% Path A: M ~ X (should use log-normal -> log(M))
% Path B: Y ~ X + M (M is predictor)
fprintf('\n--- Running lme_mediation ---\n');
res = lme_mediation(tbl, 'Y ~ X + (1|Name)', 'X', 'M', 'distM', 'log-normal', 'verbose', false);

betaA_pkg = res.paths.Estimate(1); % Path A
betaB_pkg = res.paths.Estimate(2); % Path B (from table, adjusted)
ind_pkg  = res.paths.Estimate(5); % Indirect

fprintf('Package Results:\n');
fprintf('Beta A (X->M): %.4f\n', betaA_pkg);
fprintf('Beta B (M->Y): %.4f\n', betaB_pkg);
fprintf('Indirect:      %.4f\n', ind_pkg);


% 3. Manual Check (The "Truth" assuming log-log consistency)
% Path A: log(M) ~ X
tbl.logM = log(tbl.M);
mdlA_man = fitlme(tbl, 'logM ~ X + (1|Name)');
betaA_man = mdlA_man.Coefficients.Estimate(strcmp(mdlA_man.Coefficients.Name, 'X'));

% Path B: Y ~ X + log(M)
% Note: We use log(M) as predictor because that's what Path A predicts
mdlB_man = fitlme(tbl, 'Y ~ X + logM + (1|Name)');
betaB_man = mdlB_man.Coefficients.Estimate(strcmp(mdlB_man.Coefficients.Name, 'logM'));

indirect_man = betaA_man * betaB_man;

fprintf('\nManual "Consistent" Results:\n');
fprintf('Beta A (X->logM): %.4f\n', betaA_man);
fprintf('Beta B (logM->Y): %.4f\n', betaB_man);
fprintf('Indirect:         %.4f\n', indirect_man);


% 4. Analyze Discrepancy
fprintf('\n--- Analysis ---\n');
fprintf('Diff A: %.4f\n', betaA_pkg - betaA_man);
fprintf('Diff B: %.4f\n', betaB_pkg - betaB_man);
