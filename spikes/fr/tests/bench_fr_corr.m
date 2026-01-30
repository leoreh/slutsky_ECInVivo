
%% Benchmark fr_corr

% Config
nUnits = 50;
nTime = 9000; % 15 mins at 0.1s bin
nShuffles = 500;

% Generate Random Spike Matrix (0 or 1, mostly 0)
spkMat = rand(nUnits, nTime) > 0.95;

fprintf('Running Benchmark: Units=%d, Time=%d, Shuffles=%d\n', nUnits, nTime, nShuffles);

% Run
tic;
cc = fr_corr(spkMat, 'nShuffles', nShuffles);
elapsed = toc;

fprintf('Elapsed Time: %.2f seconds\n', elapsed);
fprintf('Time per shuffle: %.4f seconds\n', elapsed / nShuffles);
