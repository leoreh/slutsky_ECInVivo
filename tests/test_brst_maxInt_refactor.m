function test_brst_maxInt_refactor()
% Test script for brst_maxInt refactor

% 1. Generate Synthetic Data
dur = 10; % seconds
fs = 20000;

% Unit 1: Sparse random spikes (Poisson-ish)
% Rate = 2 Hz
u1 = sort(rand(1, 20) * dur);

% Unit 2: Bursts
% Bursts at 1s, 3s, 5s, 7s, 9s.
% Burst: 5 spikes, ISI = 0.005s (200 Hz)
burst_starts = 1:2:9;
u2 = [];
for t = burst_starts
    b_spikes = t + (0:4)*0.005;
    u2 = [u2, b_spikes];
end
% Add some noise spikes
noise = sort(rand(1, 10) * dur);
u2 = sort([u2, noise]);

spktimes = {u1, u2};

% 2. Run Function
fprintf('Running brst_maxInt with plotting...\n');
try
    brst = brst_maxInt(spktimes, 'flgPlot', true, 'flgForce', true);

    % 3. Verify
    detected_u2 = brst.detect(2);
    expected_u2 = length(burst_starts);

    fprintf('Unit 2 Detected Bursts: %d (Expected: %d)\n', detected_u2, expected_u2);

    if detected_u2 == expected_u2
        fprintf('SUCCESS: Burst count matches.\n');
    else
        fprintf('FAILURE: Burst count mismatch.\n');
    end

    % Check Unit 1 (should be 0 or very few)
    fprintf('Unit 1 Detected Bursts: %d\n', brst.detect(1));

    fprintf('Please check the generated figure for red burst overlays.\n');

catch ME
    fprintf('ERROR: %s\n', ME.message);
    fprintf('Stack trace:\n');
    for k = 1:length(ME.stack)
        fprintf('  %s: line %d\n', ME.stack(k).name, ME.stack(k).line);
    end
end

end
