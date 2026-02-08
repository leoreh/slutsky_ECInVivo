%% Test match_qntl Aggregation Logic (Updated)

% Create synthetic data
rng(42);
n = 1000;
fr = lognrnd(1, 1, n, 1); % Log-normal distribution
fr(1:10) = 0; % Add some zeros to test robustness

tbl = table();
tbl.Group = repmat({'TestGroup'}, n, 1);
tbl.Name = repmat({'TestAnimal'}, n, 1);
tbl.Day = [repmat({'BSL'}, n/2, 1); repmat({'BAC3'}, n/2, 1)];
tbl.fr = fr;
tbl.frBspk = fr; % Just reuse for testing
tbl.frSspk = fr;
tbl.pBspk = rand(n, 1);
tbl.Group = categorical(tbl.Group);
tbl.Name = categorical(tbl.Name);
tbl.Day = categorical(tbl.Day);

nBins = 10;

% Test Mean (Default Sort: 'fr')
fprintf('Testing MEAN...\n');
tblMean = match_qntl(tbl, nBins, 'flgPool', true, 'var', 'fr', 'avgType', 'mean');
disp(head(tblMean));

% Test Median (Default Sort: 'fr')
fprintf('Testing MEDIAN...\n');
tblMedian = match_qntl(tbl, nBins, 'flgPool', true, 'var', 'fr', 'avgType', 'median');
disp(head(tblMedian));

% Test Geometric Mean
fprintf('Testing GEOMEAN...\n');
tblGeomean = match_qntl(tbl, nBins, 'flgPool', true, 'var', 'fr', 'avgType', 'geomean');
disp(head(tblGeomean));

% Verify Differences
diffv = tblMean.fr(1) - tblMedian.fr(1);
fprintf('Difference between Mean and Median (Bin 1): %f\n', diffv);

if diffv > 0
    fprintf('SUCCESS: Mean > Median for log-normal-like data.\n');
else
    fprintf('WARNING: Mean <= Median. Check data distribution.\n');
end

% Verify Auto-Processing of Variables
if ismember('frBspk', tblMean.Properties.VariableNames)
    fprintf('SUCCESS: frBspk was automatically processed.\n');
else
    fprintf('ERROR: frBspk was NOT processed.\n');
end
