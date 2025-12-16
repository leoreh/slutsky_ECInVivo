
function test_tbl_transform()

% Create sample data
rng(42);
n = 100;
Group = [repmat({'A'}, n/2, 1); repmat({'B'}, n/2, 1)];
Day = [repmat({'D1'}, n/4, 1); repmat({'D2'}, n/4, 1); repmat({'D1'}, n/4, 1); repmat({'D2'}, n/4, 1)];
ValNormal = randn(n, 1) + 10; % Normal distribution
ValSkewed = exp(randn(n, 1)); % Log-normal (skewed)

% Add some NaNs
ValNormal(1) = NaN;
ValSkewed(5) = NaN;

tbl = table(Group, Day, ValNormal, ValSkewed);
tbl.Group = categorical(tbl.Group);
tbl.Day = categorical(tbl.Day);

fprintf('Running tests for tbl_transform...\n');

%% Test 1: Normalization Logic
% Normalization should make the reference group mean = 100
fprintf('\nTest 1: Normalization Logic... ');
try
    tblNorm = tbl_transform(tbl, 'flgNorm', true, 'varNorm', 'Group', 'varsGrp', {'Day'}, 'flgZ', false, 'flgLog', false);

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

%% Test 2: Global vs Local Log Transformation
% Check if skewness check is global. Create data where one group is skewed
% but pooled is not (or vice versa), or simply ensure consistency across groups.
fprintf('Test 2: Log Transformation... ');

% Force highly skewed data
tbl.ValSkewed = 10.^randn(n,1);

tblLog = tbl_transform(tbl, 'flgLog', true, 'flgZ', false, 'flgNorm', false, 'skewThr', 0.5);

% Check if log was applied (values should be small)
if max(tblLog.ValSkewed) < 10
    fprintf('PASS (Log likely applied)\n');
else
    fprintf('FAIL (Max value %.4f too high, log not applied)\n', max(tblLog.ValSkewed));
end

%% Test 3: Z-scoring
fprintf('Test 3: Z-scoring... ');
tblZ = tbl_transform(tbl, 'flgZ', true, 'flgLog', false, 'flgNorm', false);
meanZ = mean(tblZ.ValNormal, 'omitnan');
stdZ = std(tblZ.ValNormal, 'omitnan');

if abs(meanZ) < 1e-4 && abs(stdZ - 1) < 1e-4
    fprintf('PASS\n');
else
    fprintf('FAIL (Mean=%.4f, Std=%.4f)\n', meanZ, stdZ);
end

%% Test 4: Z-scoring with Grouping
fprintf('Test 4: Z-scoring per Group... ');
tblZG = tbl_transform(tbl, 'flgZ', true, 'varsGrp', {'Group'}, 'flgLog', false);

idxA = tblZG.Group == 'A';
meanA = mean(tblZG.ValNormal(idxA), 'omitnan');
stdA = std(tblZG.ValNormal(idxA), 'omitnan');

if abs(meanA) < 1e-4 && abs(stdA - 1) < 1e-4
    fprintf('PASS\n');
else
    fprintf('FAIL (Group A Mean=%.4f, Std=%.4f)\n', meanA, stdA);
end

%% Test 5: NaN Handling
fprintf('Test 5: NaN Handling... ');
if any(isnan(tblZ.ValNormal))
    fprintf('PASS (NaNs preserved)\n');
else
    fprintf('FAIL (NaNs lost or filled)\n');
end

end
