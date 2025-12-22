% Test brst_isiValley

%% Generate Synthetic Data
% 10 units
nUnits = 10;
spktimes = cell(1, nUnits);

for i = 1:nUnits
    % Burst mode: High freq (ISI ~ 0.005s)
    % Tonic mode: Low freq (ISI ~ 0.1s)

    nSpikes = 2000;
    isis = [];

    % Mix 30% burst ISIs and 70% tonic ISIs
    nBurst = round(0.3 * nSpikes);
    nTonic = nSpikes - nBurst;

    % Log-normal distributions
    burst_isis = 10.^(randn(nBurst, 1)*0.2 + log10(0.005));
    tonic_isis = 10.^(randn(nTonic, 1)*0.3 + log10(0.1));

    spktimes{i} = cumsum([burst_isis; tonic_isis]);
end

%% Run Detection
fprintf('Running brst_isiValley...\n');
[isiValley, stats] = brst_isiValley(spktimes, 'flgPlot', false);

fprintf('Detected Valley: %.4f s (%.1f ms)\n', isiValley, isiValley*1000);

%% Validate
expected_valley_range = [0.01, 0.05];
if isiValley > expected_valley_range(1) && isiValley < expected_valley_range(2)
    fprintf('SUCCESS: Valley is within expected range.\n');
else
    fprintf('FAILURE: Valley %.4f is outside expected range [%.4f, %.4f]\n', ...
        isiValley, expected_valley_range(1), expected_valley_range(2));
end
