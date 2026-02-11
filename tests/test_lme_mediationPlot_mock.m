function test_lme_mediationPlot_mock()
% TEST_LME_MEDIATIONPLOT_MOCK Tests lme_mediationPlot with synthetic data.

%% Test 1: Numeric X (Continuous Treatment)
fprintf('Test 1: Numeric X...\n');

N = 100;
X = randn(N, 1);
M = 0.5 * X + randn(N, 1) * 0.5;
Y = 0.3 * X + 0.6 * M + randn(N, 1) * 0.5;

% Mock res structure
res = struct();
RowNames = {'Path A (X->M)'; 'Path B (M->Y|X)'; 'Path C (Total X->Y)'; 'Path C'' (Direct X->Y|M)'; 'Mediation (Sobel)'};
Variable = {'X'; 'M'; 'X'; 'X'; 'Indirect'};
Estimate = [0.5; 0.6; 0.6; 0.3; 0.3];
SE = [0.1; 0.1; 0.1; 0.1; 0.1];
PValue = [0.001; 0.001; 0.001; 0.05; 0.01];

res.paths = table(Variable, Estimate, SE, PValue, 'RowNames', RowNames);

res.plot = struct();
res.plot.X = X;
res.plot.M = M;
res.plot.Y = Y;
res.plot.Y_part_M = Y; % Dummy
res.plot.Y_part_X = Y; % Dummy

try
    lme_mediationPlot(res);
    fprintf('  PASS: Numeric X plotted successfully.\n');
    close(gcf);
catch ME
    fprintf('  FAIL: %s\n', ME.message);
end

%% Test 2: Categorical X (Group Treatment)
fprintf('\nTest 2: Categorical X...\n');

X_cat = repmat({'Ctrl', 'Treat'}, N/2, 1);
X_cat = X_cat(randperm(N));
% Convert to categorical often used in tables
X_cat = categorical(X_cat);

res.plot.X = X_cat;

try
    lme_mediationPlot(res);
    fprintf('  PASS: Categorical X plotted successfully.\n');
    close(gcf);
catch ME
    fprintf('  FAIL: %s\n', ME.message);
    disp(ME.stack(1));
end

end
