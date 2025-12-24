
% Test script for mea_tAlign.m

clear; clc;

% --- 1. Synthesize Data ---
% Create a mock 'v' struct array with 2 files.
% File 1: idxPert=10, nBins=20. (9 pre, 10 post)
% File 2: idxPert=5,  nBins=20. (4 pre, 15 post)

% Intersection (flgEdge=true):
% minPre = min(9, 4) = 4
% minPost = min(10, 15) = 10
% Result: 4 Pre, 10 Post. Total 15 bins.

% Union (flgEdge=false):
% maxPre = max(9, 4) = 9
% maxPost = max(10, 15) = 15
% Result: 9 Pre, 15 Post. Total 25 bins.

binSize = 1;
nUnits = 2;
nBins = 20;

% File 1
idxPert1 = 10;
v(1).fr.fr = ones(nUnits, nBins); % Use ones/twos to verify data not shuffled
v(1).fr.info.idxPert = idxPert1;
v(1).fr.info.binSize = binSize;
v(1).dyn.rate = ones(nUnits, nBins);

% File 2
idxPert2 = 5;
v(2).fr.fr = ones(nUnits, nBins) * 2;
v(2).fr.info.idxPert = idxPert2;
v(2).fr.info.binSize = binSize;
v(2).dyn.rate = ones(nUnits, nBins) * 2;

% Define varMap
varMap.fr = 'fr.fr';
varMap.rate = 'dyn.rate';
pertMap = 'fr.info.idxPert';

% --- Test 1: flgEdge = true (Default) -> INTERSECTION ---
fprintf('\n--- TEST 1: flgEdge = true (INTERSECTION) ---\n');
try
    [vOut1, tGlobal1] = mea_tAlign(v, varMap, pertMap, true);

    % Expected: minPre=4, minPost=10. Total 15 bins.
    expectedBins = 15;

    fprintf('Expected Bins: %d\n', expectedBins);
    fprintf('Actual Bins:   %d\n', length(tGlobal1));

    if length(tGlobal1) == expectedBins
        fprintf('SUCCESS: Correct intersection length.\n');
    else
        fprintf('FAILURE: Incorrect length. Got %d\n', length(tGlobal1));
    end

    % Verify no NaNs (since it's intersection of valid data regions, assuming input is full)
    if ~any(isnan(vOut1(1).fr.fr), 'all') && ~any(isnan(vOut1(2).fr.fr), 'all')
        fprintf('SUCCESS: No NaNs in output.\n');
    else
        fprintf('FAILURE: NaNs detected in intersection mode.\n');
    end

catch ME
    fprintf('FAILURE: Error running mea_tAlign: %s\n', ME.message);
end


% --- Test 2: flgEdge = false -> UNION ---
fprintf('\n--- TEST 2: flgEdge = false (UNION) ---\n');
try
    [vOut2, tGlobal2] = mea_tAlign(v, varMap, pertMap, false);

    % Expected: maxPre=9, maxPost=15. Total 25 bins.
    expectedBins = 25;

    fprintf('Expected Bins: %d\n', expectedBins);
    fprintf('Actual Bins:   %d\n', length(tGlobal2));

    if length(tGlobal2) == expectedBins
        fprintf('SUCCESS: Correct union length.\n');
    else
        fprintf('FAILURE: Incorrect length. Got %d\n', length(tGlobal2));
    end

    % Verify NaNs present (padding)
    % File 1 (short post) should have NaNs at end.
    % File 2 (short pre) should have NaNs at start.
    if any(isnan(vOut2(1).fr.fr(:, end-2:end)), 'all')
        fprintf('SUCCESS: File 1 padded at end.\n');
    else
        fprintf('FAILURE: File 1 NOT padded at end.\n');
    end

    if any(isnan(vOut2(2).fr.fr(:, 1:2)), 'all')
        fprintf('SUCCESS: File 2 padded at start.\n');
    else
        fprintf('FAILURE: File 2 NOT padded at start.\n');
    end

catch ME
    fprintf('FAILURE: Error running mea_tAlign: %s\n', ME.message);
end
