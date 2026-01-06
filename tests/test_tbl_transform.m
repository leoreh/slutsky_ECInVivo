
function test_tbl_transform()

% Create sample data
rng(42);
n = 100;
Group = [repmat({'A'}, n/2, 1); repmat({'B'}, n/2, 1)];
Day = [repmat({'D1'}, n/4, 1); repmat({'D2'}, n/4, 1); repmat({'D1'}, n/4, 1); repmat({'D2'}, n/4, 1)];
ValNormal = randn(n, 1) + 10; % Normal distribution
ValSkewed = exp(randn(n, 1)); % Log-normal (skewed)
ValLogit = rand(n, 1); % Proportion [0, 1]

% Add some NaNs
ValNormal(1) = NaN;
ValSkewed(5) = NaN;
ValLogit(10) = NaN;

tbl = table(Group, Day, ValNormal, ValSkewed, ValLogit);
tbl.Group = categorical(tbl.Group);
tbl.Day = categorical(tbl.Day);

fprintf('Running tests for tbl_transform...\n');

%% Test 1: Normalization Logic
% Normalization should make the reference group mean = 100
fprintf('\nTest 1: Normalization Logic... ');
try
    % Old: 'flgNorm', true, 'varNorm', 'Group'
    % New: 'varNorm', 'Group' (implies normalization)
    tblNorm = tbl_transform(tbl, 'varNorm', 'Group', 'varsGrp', {'Day'}, 'flgZ', false, 'logBase', []);

    % Check if Group A (Ref) has mean 100 within Day D1
    idxRef = tblNorm.Group == 'A' & tblNorm.Day == 'D1';
    meanRef = mean(tblNorm.ValNormal(idxRef), 'omitnan');

    if abs(meanRef - 100) < 1e-4
        fprintf('PASS\n');
    else
        fprintf('FAIL (Mean is %.4f, expected 100)\n', meanRef);
    end
catch ME
    fprintf('ERROR: %s\n', ME.message);
end

%% Test 2: Log Transformation (Base 10)
fprintf('Test 2: Log10 Transformation... ');

try
    % Force highly skewed data
    tbl.ValSkewed = 10.^randn(n,1);

    % Old: 'flgLog', true
    % New: 'logBase', 10
    tblLog = tbl_transform(tbl, 'logBase', 10, 'flgZ', false, 'skewThr', 0.5, 'varsInc', {'ValSkewed'});

    % Check if log was applied (values should be small)
    if max(tblLog.ValSkewed) < 15
        fprintf('PASS (Log likely applied)\n');
    else
        fprintf('FAIL (Max value %.4f too high, log not applied)\n', max(tblLog.ValSkewed));
    end
catch ME
    fprintf('ERROR: %s\n', ME.message);
end

%% Test 3: Log Transformation (Base e)
fprintf('Test 3: Log_e Transformation... ');
try
    tbl.ValSkewed = exp(randn(n, 1));
    tblLogE = tbl_transform(tbl, 'logBase', 'e', 'flgZ', false, 'skewThr', -Inf, 'varsInc', {'ValSkewed'});

    % If x = e^y, then ln(x) = y. Mean of y should be close to 0, std close to 1
    if abs(mean(tblLogE.ValSkewed, 'omitnan')) < 0.5
        fprintf('PASS\n');
    else
        fprintf('FAIL (Mean=%.4f) \n', mean(tblLogE.ValSkewed, 'omitnan'));
    end
catch ME
    fprintf('ERROR: %s\n', ME.message);
end

%% Test 4: Logit Transformation
fprintf('Test 4: Logit Transformation... ');
try
    % Old: 'flgLogit', true
    % New: 'logBase', 'logit'
    tblLogit = tbl_transform(tbl, 'logBase', 'logit', 'flgZ', false, 'varsInc', {'ValLogit'});

    % Check range (logits map [0,1] to [-inf, inf])
    % Just check if values changed and are roughly symmetric
    if min(tblLogit.ValLogit) < 0 && max(tblLogit.ValLogit) > 0
        fprintf('PASS\n');
    else
        fprintf('FAIL (Range=[%.2f, %.2f]) \n', min(tblLogit.ValLogit), max(tblLogit.ValLogit));
    end
catch ME
    fprintf('ERROR: %s\n', ME.message);
end

%% Test 5: Z-scoring
fprintf('Test 5: Z-scoring... ');
tblZ = tbl_transform(tbl, 'flgZ', true, 'logBase', [], 'varNorm', '');
meanZ = mean(tblZ.ValNormal, 'omitnan');
stdZ = std(tblZ.ValNormal, 'omitnan');

if abs(meanZ) < 1e-4 && abs(stdZ - 1) < 1e-4
    fprintf('PASS\n');
else
    fprintf('FAIL (Mean=%.4f, Std=%.4f)\n', meanZ, stdZ);
end

%% Test 6: Z-scoring per Group
fprintf('Test 6: Z-scoring per Group... ');
tblZG = tbl_transform(tbl, 'flgZ', true, 'varsGrp', {'Group'}, 'logBase', []);

idxA = tblZG.Group == 'A';
meanA = mean(tblZG.ValNormal(idxA), 'omitnan');
stdA = std(tblZG.ValNormal(idxA), 'omitnan');

if abs(meanA) < 1e-4 && abs(stdA - 1) < 1e-4
    fprintf('PASS\n');
else
    fprintf('FAIL (Group A Mean=%.4f, Std=%.4f)\n', meanA, stdA);
end

%% Test 7: NaN Handling
fprintf('Test 7: NaN Handling... ');
if any(isnan(tblZ.ValNormal))
    fprintf('PASS (NaNs preserved)\n');
else
    fprintf('FAIL (NaNs lost or filled)\n');
end

end
