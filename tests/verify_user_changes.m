
% Verify user simplified functions

% 1. Test n2chunks
n = 130;
chunksize = 60;
% 'exclude' -> [1 60; 61 120]
c = n2chunks('n', n, 'chunksize', chunksize, 'lastChunk', 'exclude');
assert(size(c, 1) == 2, 'n2chunks exclude failed size');
assert(c(2, 2) == 120, 'n2chunks exclude failed end');

% 'extend' -> [1 60; 61 130]
c_ext = n2chunks('n', n, 'chunksize', chunksize, 'lastChunk', 'extend');
assert(size(c_ext, 1) == 2, 'n2chunks extend failed size');
assert(c_ext(2, 2) == 130, 'n2chunks extend failed end');

% 2. Test times2rate
spktimes = { [10; 70; 125] };
% binSize 60. Max time 125.
% Exclude -> bins 0-60, 60-120. (120-125 dropped)
[r, edges, cents] = times2rate(spktimes, 'binSize', 60);

assert(size(r, 2) == 2, 'times2rate failed bin count');
assert(edges(end, 2) == 120, 'times2rate failed last edge');

disp('User Simplifications Verified.');
