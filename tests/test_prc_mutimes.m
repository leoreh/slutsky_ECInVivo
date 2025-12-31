function test_prc_mutimes()
% TEST_PRC_MUTIMES Verifies that mutimes are correctly integrated into prc_calc

% 1. Setup synthetic data
% -----------------------
duration = 10; % seconds
nUnits = 5;
nMutimes = 3; % 3 tetrodes of MUA

% Create some random spike train for units
spktimes = cell(nUnits, 1);
all_unit_spikes = [];
for i = 1:nUnits
    spktimes{i} = sort(rand(200, 1) * duration);
    all_unit_spikes = [all_unit_spikes; spktimes{i}];
end

% Create mutimes:
% - Contains copies of unit spikes (simulating unsorted spikes being present)
% - Contains unique "background" spikes
mutimes = cell(nMutimes, 1);
extra_spikes_count = 0;

for i = 1:nMutimes
    % Add some spikes from random units (to be filtered out)
    unit_subset = spktimes{randi(nUnits)};

    % Add some unique background noise
    noise_spikes = sort(rand(50, 1) * duration);

    % Combine
    mutimes{i} = sort([unit_subset; noise_spikes]);

    extra_spikes_count = extra_spikes_count + length(noise_spikes);
end

fprintf('Total unit spikes: %d\n', length(all_unit_spikes));
fprintf('Total MUA unique background spikes: %d\n', extra_spikes_count);

% 2. Run WITHOUT mutimes
% ----------------------
fprintf('\nRunning baseline PRC (no mutimes)...\n');
prc_base = prc_calc(spktimes, 'winLim', [0 duration], 'nShuffles', 10, 'flgSave', false);

% 3. Run WITH mutimes
% -------------------
fprintf('\nRunning enhanced PRC (with mutimes)...\n');
prc_mu = prc_calc(spktimes, 'mutimes', mutimes, 'winLim', [0 duration], 'nShuffles', 10, 'flgSave', false);

% 4. Verify Results
% -----------------

% Check 1: Info field
if isfield(prc_mu.info, 'nMutimes')
    fprintf('[PASS] nMutimes field present in info\n');
else
    fprintf('[WARN] nMutimes field MISSING in info\n');
end

% Check 2: Population Rate Impact
% We can't directly check the internal popRate, but we can verify that the
% underlying raster used for calculation was likely denser if the results change.
% However, prc0 might go up or down depending on the correlations.
% A better check would be if we could return the popRate or check raster size,
% but for now let's just ensure it runs and produces different results.

diff_prc = max(abs(prc_base.prc0 - prc_mu.prc0));
if diff_prc > 1e-6
    fprintf('[PASS] Results differ significantly (max diff: %.4f), implying mutimes were used.\n', diff_prc);
else
    fprintf('[WARN] Results are identical! Mutimes might have been ignored or filtered completely.\n');
end

% Check 3: nUnits should be same (reference units shouldn't increase)
if prc_base.info.nUnits == prc_mu.info.nUnits
    fprintf('[PASS] nUnits is consistent (%d)\n', prc_base.info.nUnits);
else
    fprintf('[FAIL] nUnits changed! Base: %d, Mu: %d\n', prc_base.info.nUnits, prc_mu.info.nUnits);
end

fprintf('\nDone.\n');
end
