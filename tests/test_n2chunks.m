
% Test n2chunks new parameter

n = 130;
chunksize = 60;

% 1. keep (default)
c_keep = n2chunks('n', n, 'chunksize', chunksize, 'lastChunk', 'keep');
% Expect: [1 60; 61 120; 121 130]
assert(size(c_keep, 1) == 3);
assert(c_keep(3, 2) == 130);

% 2. exclude
c_excl = n2chunks('n', n, 'chunksize', chunksize, 'lastChunk', 'exclude');
% Expect: [1 60; 61 120]
assert(size(c_excl, 1) == 2);
assert(c_excl(2, 2) == 120);

% 3. extend
c_extend = n2chunks('n', n, 'chunksize', chunksize, 'lastChunk', 'extend');
% Expect: [1 60; 61 130]
assert(size(c_extend, 1) == 2);
assert(c_extend(2, 2) == 130);

% 4. Case where n divides perfectly
n2 = 120;
c_perf = n2chunks('n', n2, 'chunksize', chunksize);
% Expect 2 chunks
assert(size(c_perf, 1) == 2);
assert(c_perf(end, 2) == 120);

c_perf_ex = n2chunks('n', n2, 'chunksize', chunksize, 'lastChunk', 'exclude');
assert(size(c_perf_ex, 1) == 2);

c_perf_ext = n2chunks('n', n2, 'chunksize', chunksize, 'lastChunk', 'extend');
assert(size(c_perf_ext, 1) == 2);

% 5. Case where n < chunksize
n3 = 40;
c_short_keep = n2chunks('n', n3, 'chunksize', 60, 'lastChunk', 'keep');
assert(size(c_short_keep, 1) == 1);
assert(c_short_keep(1, 2) == 40);

c_short_excl = n2chunks('n', n3, 'chunksize', 60, 'lastChunk', 'exclude');
assert(isempty(c_short_excl));

disp('ALL TESTS PASSED for n2chunks');
