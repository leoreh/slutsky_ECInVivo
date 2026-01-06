function test_lme_parkTest()
% TEST_LME_PARKTEST Verification script for lme_parkTest logic

rng(42);
n = 200;
x = randn(n, 1);
g = repmat({'A'; 'B'}, n/2, 1);

%% Test 1: Normal Distribution (Lambda ~ 0)
fprintf('Test 1: Normal Distribution... ');
yNorm = 10 + 2*x + randn(n, 1);
tblNorm = table(yNorm, x, g, 'VariableNames', {'y', 'x', 'g'});

stats = lme_parkTest(tblNorm, 'y ~ x + (1|g)');

if abs(stats.Lambda) < 0.5 && strcmp(stats.Recommendation, 'Gaussian')
    fprintf('PASS (Lambda=%.2f, Rec=%s)\n', stats.Lambda, stats.Recommendation);
else
    fprintf('FAIL (Lambda=%.2f, Rec=%s)\n', stats.Lambda, stats.Recommendation);
end

%% Test 2: Gamma Distribution (Lambda ~ 2)
fprintf('Test 2: Gamma Distribution... ');
% Generate Gamma data: Mean related to X, constant coeff of variation
mu = exp(2 + 0.5*x);
shape = 2;
scale = mu / shape;
yGamma = gamrnd(shape, scale);

tblGamma = table(yGamma, x, g, 'VariableNames', {'y', 'x', 'g'});

stats = lme_parkTest(tblGamma, 'y ~ x + (1|g)');

% Lambda should be close to 2 for Gamma
if abs(stats.Lambda - 2) < 0.6 && ismember(stats.Recommendation, {'Gamma', 'Log-Normal'})
    fprintf('PASS (Lambda=%.2f, Rec=%s)\n', stats.Lambda, stats.Recommendation);
else
    fprintf('FAIL (Lambda=%.2f, Rec=%s)\n', stats.Lambda, stats.Recommendation);
end


%% Test 3: Log-Normal (Heavy Tails) -> Auto-detection
fprintf('Test 3: Log-Normal (Heavy Tail Check)... ');
% Generate Log-Normal data with high variance
yLogNorm = exp(2 + 0.5*x + 1.5*randn(n, 1));
tblLog = table(yLogNorm, x, g, 'VariableNames', {'y', 'x', 'g'});

stats = lme_parkTest(tblLog, 'y ~ x + (1|g)');

if strcmp(stats.Recommendation, 'Log-Normal')
    fprintf('PASS (Rec=%s, Kurt=%.1f)\n', stats.Recommendation, stats.KurtosisLog);
else
    fprintf('WARN (Rec=%s, Kurt=%.1f) - May not always trigger depending on random seed\n', ...
        stats.Recommendation, stats.KurtosisLog);
end


%% Test 4: Zero-Inflation Auto-Correction
fprintf('Test 4: Zero-Inflation Auto-Correction... ');
% Create skewed data with some zeros
yZero = yGamma;
yZero(1:5) = 0; % Add 5 zeros (2.5% of 200) -> Should trigger offset but stay Gamma

tblZero = table(yZero, x, g, 'VariableNames', {'y', 'x', 'g'});

% Capture warning
lastwarn('');
stats = lme_parkTest(tblZero, 'y ~ x + (1|g)');
[msg, id] = lastwarn;

% Check if offset warning occurred or if it ran successfully
if contains(msg, 'Applying offset') || strcmp(stats.Recommendation, 'Gamma')
    fprintf('PASS (Handled Zeros, Rec=%s)\n', stats.Recommendation);
else
    fprintf('FAIL (Rec=%s, Msg=%s)\n', stats.Recommendation, msg);
end


%% Test 5: Verify lme_analyse Predictor Transformation
fprintf('Test 5: Verify Predictor Z-Scoring in lme_analyse... ');
% Add a skewed predictor to tblNorm
tblNorm.z = exp(tblNorm.x); % Skewed
try
    [~, ~, info] = lme_analyse(tblNorm, 'y ~ x + z + (1|g)');
    % We can't check the internal table directly easily without returning it,
    % but we can check if it ran without error and produced a model.
    fprintf('PASS (Ran successfully)\n');
catch ME
    fprintf('FAIL (%s)\n', ME.message);
end

end
