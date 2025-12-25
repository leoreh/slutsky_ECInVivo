%% Test spk2cCa
% Verification script for the spk2cCa function.

% 1. Setup Synthetic Data
dt = 0.001; % 1 ms
T  = 60;    % 60 seconds (1 minute) to match default binSize
t  = 0:dt:T;

% Unit 1: Sparse firing (2 Hz)
lambda1 = 2;
spk1 = [];
currT = 0;
while currT < T
    isi = -log(rand) / lambda1;
    currT = currT + isi;
    if currT < T
        spk1(end+1) = currT;
    end
end

% Unit 2: Bursty (Train of 10 spikes at 100 Hz every second)
spk2 = [];
for i = 1:60 % 60 bursts (1 per sec)
    startT = i - 0.5;
    for j = 0:9
        spk2(end+1) = startT + j * 0.01; % 10ms ISI = 100 Hz
    end
end

spktimes = {spk1, spk2};

% 2. Run spk2cCa
fprintf('Running spk2cCa (Binned)...\n');
% Disable plotting for automated logic check
mCa = spk2cCa(spktimes, 'dt', dt, 'binSize', 60, 'flgPlot', false);

% 3. Verify Dimensions
% Should be 1 bin for 60s since binSize=60? Or 2 if edges push over?
% n2chunks logic: if T=60, chunksize=60 -> [1 60]. 1 chunk.
fprintf('Number of bins: %d\n', length(mCa.time));
assert(size(mCa.cyto, 1) == 2, 'Cyto matrix rows mismatch');

% 4. Verify Logic
fprintf('Total Accumulation (Psi) in first bin:\n');
% Check if we have data in first bin
if ~isempty(mCa.total)
    fprintf('Unit 1 (Sparse): %.4f\n', mCa.total(1, 1));
    fprintf('Unit 2 (Bursty): %.4f\n', mCa.total(2, 1));

    % Expect Unit 2 to be much higher than Unit 1
    if mCa.total(2, 1) > mCa.total(1, 1) * 5
        fprintf('PASS: Bursty unit accumulation is significantly higher.\n');
    else
        error('FAIL: Bursty unit accumulation is not significantly higher.');
    end
end

fprintf('Test spk2cCa COMPLETED SUCCESSFULLY.\n');
