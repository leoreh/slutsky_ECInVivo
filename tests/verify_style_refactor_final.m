% Verify refactored functions (Final Check)

% 1. Test n2chunks style and logic
n = 130;
chunksize = 60;
c = n2chunks('n', n, 'chunksize', chunksize, 'lastChunk', 'exclude');
assert(size(c, 1) == 2, 'n2chunks exclude failed size');
assert(c(2, 2) == 120, 'n2chunks exclude failed end');

% 2. Test times2rate
spktimes = { [10; 70; 125] };
[r, edges, cents] = times2rate(spktimes, 'binSize', 60);

% Should have 2 bins: 0-60, 60-120.
% Rate in bin 1 (0-60): 1 spike (at 10). Rate = 1/60 = 0.0167
% Rate in bin 2 (60-120): 1 spike (at 70). Rate = 1/60.
assert(size(r, 2) == 2, 'times2rate failed bin count');
assert(edges(end, 2) == 120, 'times2rate failed last edge');
assert(abs(r(1,1) - 1/60) < 1e-6, 'Rate calculation failed');

disp('Style Refactor Final Verification Passed');
