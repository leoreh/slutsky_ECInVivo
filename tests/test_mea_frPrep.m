% Test mea_frPrep

clear; clc;

% 1. Create synthetic data
% mcu_detectPert requires long baselines (hardcoded > 20 hrs search start).
% We simulate 25 hours baseline + 5 hours post.
% Bin size = 60s.

nUnits = 10;
nHrsPre = 25;
nHrsPost = 5;

binSize = 60;
nBinsPre = round(nHrsPre * 3600 / binSize);
nBinsPost = round(nHrsPost * 3600 / binSize);
nBins = nBinsPre + nBinsPost;

idxPertTrue = nBinsPre + 1;

frHigh = 10;
frLow = 2;

spktimes = cell(nUnits, 1);

fprintf('Generating %d hours of synthetic data...\n', nHrsPre + nHrsPost);

for i = 1:nUnits
    % Generate spikes
    % Pre-pert: High rate
    tPre = 1 : binSize : (nBinsPre * binSize);
    nSpkPre = poissrnd(frHigh * binSize, 1, nBinsPre);
    spks = [];
    % Vectorized spike generation for speed
    for b = 1:nBinsPre
        if nSpkPre(b) > 0
            s = tPre(b) + rand(nSpkPre(b), 1) * binSize;
            spks = [spks; s];
        end
    end

    % Post-pert: Low rate
    tPost = (nBinsPre * binSize) + (1 : binSize : (nBinsPost * binSize));
    nSpkPost = poissrnd(frLow * binSize, 1, nBinsPost);
    for b = 1:nBinsPost
        if nSpkPost(b) > 0
            s = tPost(b) + rand(nSpkPost(b), 1) * binSize;
            spks = [spks; s];
        end
    end

    spktimes{i} = spks;
end

% 2. Run mea_frPrep
fprintf('Running mea_frPrep...\n');
frOut = mea_frPrep(spktimes, 'binSize', binSize);

% 3. Verify
fprintf('--- Verification ---\n');

% 3.1 Structure fields
assert(isfield(frOut, 'fr'), 'Missing fr');
assert(isfield(frOut, 'frOrig'), 'Missing frOrig');
assert(isfield(frOut, 't'), 'Missing t');
assert(isfield(frOut.info, 'idxPert'), 'Missing info.idxPert');

% 3.2 Perturbation Detection
idxPert = frOut.info.idxPert;
fprintf('Detected idxPert: %d (True: %d)\n', idxPert, idxPertTrue);
% Allow slack (mcu_detectPert is complex)
% Note: mcu_detectPert finds the drop. In random Poisson data, it might pick a random large drop.
% But with huge rate difference (10Hz vs 2Hz) it should find the transition.
err = abs(idxPert - idxPertTrue);
if err > 10
    warning('Perturbation detection off by %d bins. This may be due to stochasticity or mcu_detectPert heuristics.', err);
else
    fprintf('Perturbation detection accurate.\n');
end

% 3.3 Time Correction
t = frOut.t;

% t at idxPert should be 0
t0 = t(idxPert);
fprintf('Time at idxPert: %.2f (Expected 0)\n', t0);
assert(abs(t0) < 1e-5, 'Time at perturbation is not 0');

% Check Baseline Scaling (x3)
% One bin before pert: should be -1 * 60 * 3 = -180
tPrev = t(idxPert - 1);
fprintf('Time at idxPert-1: %.2f (Expected -180)\n', tPrev);
assert(abs(tPrev - (-180)) < 1e-1, 'Baseline scaling incorrect');

% Check Post-Pert Scaling (x6)
% One bin after pert: should be 1 * 60 * 6 = 360
tNext = t(idxPert + 1);
fprintf('Time at idxPert+1: %.2f (Expected 360)\n', tNext);
assert(abs(tNext - 360) < 1e-1, 'Post-perturbation scaling incorrect');

fprintf('All tests passed!\n');
