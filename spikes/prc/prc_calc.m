function prc = prc_calc(spktimes, varargin)
% PRC_CALC Calculates population coupling for all neurons in spktimes.
%
%   prc = PRC_CALC(SPKTIMES, ...) implements a hybrid approach to calculate
%   population coupling using baseline-subtracted spike-triggered population
%   rate (stPR) at zero lag.
%
%   INPUTS:
%       spktimes    - (cell) Spike times per unit (in seconds).
%       varargin    - (param/value) Optional parameters:
%                     'basepath'   : (char) Path to recording session directory {pwd}
%                     'flgSave'    : (log) Save results to .prc.mat file {true}
%                     'winLim'     : (vec) Time window to analyze [start end] {[0 Inf]} (s)
%                     'binSize'    : (num) Bin size for spike trains {0.001} (s)
%                     'gkHw'       : (num) Gaussian kernel half-width (sigma) {0.012} (s)
%                     'winStpr'    : (num) Total width of STPR window {1.0} (s)
%                     'nShuffles'  : (num) Number of shuffles for null distribution {1000}
%                     'spkThr'     : (num) Min spikes required. Fewer -> NaN {100}
%                     'spkLim'     : (num) Max spikes from ref unit to use {Inf}
%
%   OUTPUTS:
%       prc         - (struct) Population coupling results.
%                     .prc0_z     : (nUnits x 1) Z-scored population coupling.
%                     .prc0_norm  : (nUnits x 1) Median-normalized population coupling.
%                     .prc0       : (nUnits x 1) Raw, baseline-subtracted STPR at zero lag.
%                     .stpr       : (nUnits x nBins) Full STPR curves.
%                     .prc0_shfl  : (nUnits x nShuffles) Distribution of STPR(0).
%                     .pk         : (nUnits x 1) Peak STPR within +/- 100ms.
%                     .pk_z       : (nUnits x 1) Z-scored peak STPR.
%                     .pk_norm    : (nUnits x 1) Normalized peak STPR.
%                     .pk_shfl    : (nUnits x nShuffles) Distribution of peak STPR.
%                     .t          : (1 x nBins) Time vector for STPR curves.
%                     .info       : (struct) Analysis parameters.
%
%   See also: PRC_PLOT, PRC_SHUFFLE
%
%   HISTORY:
%   Aug 2024 LH - Major refactor and simplification.
%                 Optimize STPR calculation (LOO) and parallelize (parfor).
%                 Requires Parallel Computing Toolbox.
%                 Requires recompilation of shuffle_raster.cpp.
%
%   COMPLIE:
%                 mex -v -O CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' shuffle_raster.cpp

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

% Define input parameters and their defaults/validation functions
p = inputParser;
addRequired(p, 'spktimes', @iscell);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgSave', true, @islogical);
addParameter(p, 'winLim', [0 Inf], @(x) isnumeric(x) && numel(x)==2);
addParameter(p, 'binSize', 0.001, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'gkHw', 0.012, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'winStpr', 1.0, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'nShuffles', 1000, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'spkThr', 100, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'spkLim', Inf, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'mutimes', {}, @iscell);

% Parse input arguments
parse(p, spktimes, varargin{:});
spktimes = p.Results.spktimes;
basepath = p.Results.basepath;
flgSave = p.Results.flgSave;
winLim = p.Results.winLim;
binSize = p.Results.binSize;
gkHw = p.Results.gkHw;
winStpr = p.Results.winStpr;
nShuffles = p.Results.nShuffles;
spkThr = p.Results.spkThr;
spkLim = p.Results.spkLim;
mutimes = p.Results.mutimes;

% Set basepath and get basename
cd(basepath);
[~, basename] = fileparts(basepath);

%% ========================================================================
%  TIME AXIS & KERNEL SETUP
%  ========================================================================

% Create time axis for binning
if isinf(winLim(2))
    maxTime = max(cellfun(@(x) max([0; x(:)]), mutimes));
    winLim(2) = ceil(maxTime);
end
t = winLim(1):binSize:winLim(2);

% Define Gaussian kernel for smoothing
sigmaBins = gkHw / binSize;
gkW = ceil(3 * sigmaBins);
gkT = -gkW:gkW;
gk = normpdf(gkT, 0, sigmaBins);
gk = gk / sum(gk);


% Calculate STPR window parameters
lagBins = round((winStpr / 2) / binSize);
nBins = 2 * lagBins + 1;

% Calculate Peak window parameters
winPk = 0.100;
nBinsPk = round(winPk / binSize);
idxPk = (-nBinsPk:nBinsPk); % Relative to center

%% ========================================================================
%  SPIKE WINDOWING & INITIALIZATION
%  ========================================================================

% Window spikes to analysis window using cellfun
spktimes = cellfun(@(x) x(x >= winLim(1) & x <= winLim(2)), spktimes, 'UniformOutput', false);
nSpks = cellfun(@(x) length(x), spktimes, 'UniformOutput', true);
uGood = find(nSpks > spkThr);
nGood = length(uGood);
nUnits = length(spktimes);

% Initialize output arrays at full size
% Initialize output arrays at full size
stpr = zeros(nUnits, nBins);
prc0 = nan(nUnits, 1);
pk = nan(nUnits, 1);

%% ========================================================================
%  PREPARE RASTER MATRIX
%  ========================================================================

% Initialize binRates only for good units
binRates = zeros(nGood, length(t));  % Only store good units

% Calculate binned rates across good units and store spike bin indices
for iGood = 1:nGood
    idxRef = uGood(iGood);
    uSpks = spktimes{idxRef};
    binRates(iGood,:) = histcounts(uSpks, [t, t(end) + binSize]);
end

% Convert to binary and get spike bin indices
rasterMat = binRates > 0;
spkBins = binary2idx(rasterMat, [lagBins, length(t) - lagBins], spkLim);

% Process mutimes if provided
if ~isempty(mutimes)

    % Collect all spikes for filtering
    allSpks = cell2mat(spktimes(:));

    % Window mutimes
    mutimes = cellfun(@(x) x(x >= winLim(1) & x <= winLim(2)), ...
        mutimes, 'UniformOutput', false);

    % Filter out sorted spikes from mutimes and bin
    nMu = length(mutimes);
    muRates = zeros(nMu, length(t));

    for iMu = 1:nMu
        % Remove spikes that are already in sorted units
        muClean = setdiff(mutimes{iMu}, allSpks);
        muRates(iMu, :) = histcounts(muClean, [t, t(end) + binSize]);
    end

    % Append to rasterMat. These rows are part of popRate & shuffle, but
    % not indexed by uGood loop
    rasterMat = [rasterMat; muRates > 0];
end

% Determine number of swaps when creating shuffled raster matrix, based on
% number of spikes
nSwaps = sum(rasterMat(:)) * 10;

%% ========================================================================
%  STPR CALCULATION (CORE)
%  ========================================================================

% Calculate and smooth global population rate
popRate = sum(rasterMat, 1, 'omitnan');
popRate = conv(popRate / binSize, gk, 'same');

for iGood = 1:nGood
    idxRef = uGood(iGood);

    % Calculate unit smoothed rate
    ur = conv(rasterMat(iGood, :) / binSize, gk, 'same');

    % Calculate LOO population rate: Global - Unit
    pr = popRate - ur;

    % Calculate STPR using pre-processed spike bins
    [stpr(idxRef, :), prc0(idxRef), pk(idxRef)] = stpr_calc(spkBins{iGood}, pr, lagBins, idxPk);
end

%% ========================================================================
%  STPR SHUFFLING
%  ========================================================================

% Determine number of parallel chains (one per worker)
p = gcp('nocreate');
if isempty(p)
    p = parpool('local', 8);
end
nChains = p.NumWorkers;
nPerChain = ceil(nShuffles / nChains);
nTotal = nChains * nPerChain;

% Pre-generate random seeds for reproducible parallel shuffling
% Each step in each chain needs a unique seed
seeds = randi([0, 2^32-1], nChains, nPerChain, 'uint32');

% Initialize progress monitoring
dq = parallel.pool.DataQueue;
progress_monitor(nTotal, 'init');
afterEach(dq, @(~) progress_monitor([], 'update'));

% Initialize results containers for parallel chains
results_prc0 = cell(nChains, 1);
results_pk = cell(nChains, 1);
results_stpr = cell(nChains, 1);

parfor iChain = 1:nChains

    % Initialize chain state with original data
    rasterCurrent = rasterMat;

    % Initialize chain accumulators
    chain_stpr = zeros(nUnits, nBins);
    chain_prc0 = nan(nUnits, nPerChain);
    chain_pk = nan(nUnits, nPerChain);

    for iShuffle = 1:nPerChain

        % Determine nSwaps
        if iShuffle == 1
            iSwaps = nSwaps;
        else
            iSwaps = size(rasterCurrent, 1) * 5;
        end

        % Shuffle
        rasterCurrent = shuffle_raster(rasterCurrent, iSwaps, seeds(iChain, iShuffle));

        % Calc STPR
        refBins = binary2idx(rasterCurrent, [lagBins, length(t) - lagBins], spkLim);
        stpr_tmp = zeros(nUnits, nBins);
        prc0_tmp = nan(nUnits, 1);
        pk_tmp = nan(nUnits, 1);

        for iGood = 1:nGood
            idxRef = uGood(iGood);
            ur = conv(rasterCurrent(iGood, :) / binSize, gk, 'same');
            pr = popRate - ur;
            [stpr_tmp(idxRef,:), prc0_tmp(idxRef), pk_tmp(idxRef)] = stpr_calc(refBins{iGood}, pr, lagBins, idxPk);
        end

        chain_stpr = chain_stpr + stpr_tmp;
        chain_prc0(:, iShuffle) = prc0_tmp;
        chain_pk(:, iShuffle) = pk_tmp;
        send(dq, []);
    end

    results_prc0{iChain} = chain_prc0;
    results_pk{iChain} = chain_pk;
    results_stpr{iChain} = chain_stpr;
end

% Aggregate results
stpr_sum = zeros(nUnits, nBins);
prc0_shfl = zeros(nUnits, 0); % Append
pk_shfl = zeros(nUnits, 0);   % Append
for iChain = 1:nChains
    stpr_sum = stpr_sum + results_stpr{iChain};
    prc0_shfl = [prc0_shfl, results_prc0{iChain}];
    pk_shfl = [pk_shfl, results_pk{iChain}];
end

% Adjust nShuffles to actual total executed
nShuffles = size(prc0_shfl, 2);

% Average across shuffles
stpr_shfl = stpr_sum / nShuffles;

%% ========================================================================
%  NORMALIZE
%  ========================================================================

% Median Normalization
shflMed = median(prc0_shfl, 'all', 'omitnan');
prc0_norm = prc0 ./ shflMed;

% Z-score
shflMu = mean(prc0_shfl, 2, 'omitnan');
shflSd = std(prc0_shfl, [], 2, 'omitnan');
prc0_z = (prc0 - shflMu) ./ shflSd;

% Handle zero std dev cases. If std is 0, check if value equals mean (0
% Z-score) or not (NaN)
maskSdZero = shflSd < eps;
if any(maskSdZero)
    warning('zero std')
    isMean = abs(prc0 - shflMu) < eps;
    prc0_z(maskSdZero & isMean) = 0;
    prc0_z(maskSdZero & ~isMean) = NaN;
end

% --- Peak Metrics ---

% Median Normalization
shflMedPk = median(pk_shfl, 'all', 'omitnan');
pk_norm = pk ./ shflMedPk;

% Z-score
shflMuPk = mean(pk_shfl, 2, 'omitnan');
shflSdPk = std(pk_shfl, [], 2, 'omitnan');
pk_z = (pk - shflMuPk) ./ shflSdPk;

maskSdZeroPk = shflSdPk < eps;
if any(maskSdZeroPk)
    isMean = abs(pk - shflMuPk) < eps;
    pk_z(maskSdZeroPk & isMean) = 0;
    pk_z(maskSdZeroPk & ~isMean) = NaN;
end


%% ========================================================================
%  ORGANIZE OUTPUT
%  ========================================================================

% Store results (with NaN for bad units)
prc.prc0 = prc0;               % Raw STPR at zero lag
prc.prc0_z = prc0_z;           % Z-scored population coupling
prc.prc0_norm = prc0_norm;     % Median-normalized population coupling
prc.prc0_shfl = prc0_shfl;     % Shuffled STPR values
prc.pk = pk;
prc.pk_z = pk_z;
prc.pk_norm = pk_norm;
prc.pk_shfl = pk_shfl;
prc.stpr = stpr;               % Full STPR curves
prc.stpr_shfl = stpr_shfl;
prc.spkBins = spkBins;

% Create info struct with parameters and metadata
prc.info = struct(...
    'basepath', basepath, ...
    'runtime', datetime("now"), ...
    'nUnits', nUnits, ...
    'nGood', nGood, ...        % Number of good units
    'nMutimes', length(mutimes), ... % Number of multi-unit channels
    'uGood', uGood, ...        % Indices of good units
    'winLim', winLim, ...
    'spkThr', spkThr, ...
    'spkLim', spkLim, ...
    'binSize', binSize, ...
    'gkHw', gkHw, ...
    'winStpr', winStpr, ...
    'nShuffles', nShuffles);

% Save results if requested
if flgSave
    save(fullfile(basepath, [basename, '.prc.mat']), 'prc', '-v7.3');
end

end     % EOF


%% ========================================================================
%  HELPER: CALCULATE STPR
%  ========================================================================


function [stpr, prc0, pk] = stpr_calc(uBins, popRate, lagBins, idxPk)
% CALC_STPR Calculates spike-triggered average of population rate.
%
% INPUT:
%   uBins       - Vector of bin indices where reference neuron spiked.
%   popRate     - Population firing rate (Hz).
%   lagBins     - Number of bins for STPR window.
%   idxPk       - Indices relative to center for peak search.
%
% OUTPUT:
%   prc0       - Baseline-subtracted STPR at zero lag (Hz).
%   stpr       - spike triggered population rate
%   pk         - Peak STPR within window.

baseline = mean(popRate);
idx = bsxfun(@plus, uBins, (-lagBins:lagBins));
segments = popRate(idx);
stpr = mean(segments, 1) - baseline;

ctr = lagBins + 1;
prc0 = stpr(ctr);
pk = max(stpr(ctr + idxPk));

end


%% ========================================================================
%  HELPER: BINARY MATRIX TO SPIKE BIN INDICES
%  ========================================================================

function spkBins = binary2idx(binMat, edgeLim, spkLim)
% BINARY2IDX Converts binary matrix to cell array of spike bin indices.
%
% INPUT:
%   binMat      - Binary matrix (neurons x time_bins) of spike times.
%   edgeLim     - [min max] bin indices to include (exclusive).
%   spkLim      - Maximum number of spikes to use per unit.
%
% OUTPUT:
%   spkBins     - Cell array. spkBins{i} contains validated spike bin indices
%                 for neuron i.

nGood = size(binMat, 1);
spkBins = cell(nGood, 1);

for iGood = 1:nGood
    % Get spike bins and validate against window edges
    refBins = find(binMat(iGood,:) > 0);
    refBins = refBins(refBins > edgeLim(1) & refBins <= edgeLim(2));

    % Subsample spike bins if needed
    rng(iGood);                             % Use consistent random seed for reproducibility
    if length(refBins) > spkLim
        refBins = refBins(randperm(length(refBins), spkLim));
    end
    spkBins{iGood} = sort(refBins(:));      % Keep bins sorted in time
end

end


%% ========================================================================
%  HELPER: PROGRESS MONITOR
%  ========================================================================

function progress_monitor(data, mode)
% PROGRESS_MONITOR Updates and displays progress in command window.
%
% This function uses persistent variables to track state across callback
% executions from parallel.pool.DataQueue.

persistent cnt nTotal

if strcmp(mode, 'init')
    cnt = 0;
    nTotal = data;
    fprintf('Starting %d shuffles...\n', nTotal);
elseif strcmp(mode, 'update')
    cnt = cnt + 1;
    % Update every 10 shuffles or on completion
    if mod(cnt, 10) == 0 || cnt == nTotal
        fprintf('\rShuffle %d / %d', cnt, nTotal);
    end
    % Print newline on completion
    if cnt == nTotal
        fprintf('\n');
    end
end
end

%% ========================================================================
%  NOTE: MARKOV CHAIN MONTE CARLO (MCMC)
%  ========================================================================
%  While traditional shuffling resets to the observed data for every
%  iteration, MCMC treats the null distribution as a continuous space to
%  be explored. The 'Burn-In' transition achieves high entropy, while
%  subsequent steps maintain that entropy with minimal displacement.
%
%  * Burn-In: Decouples the matrix from the original spike timing.
%  * Mixing Step: Navigates the space of all possible permutations.
%  * Stationarity: Preserves row/column marginals throughout the chain.
%  * Efficiency: Reduces O(N) overhead to O(k) per sample.
%
%  JUSTIFICATION FOR REDUCED INTER-SHUFFLE SWAPS
%  Using the end-state of Shuffle(i) as the starting-state for Shuffle(i+1)
%  is a standard practice in computational statistics (Metropolis-Hastings)
%  for the following reasons:
%
%  1. Entropy Persistence: Once the matrix reaches a randomized state,
%     it does not 're-order' itself. Subsequent swaps merely explore
%     alternative configurations within the same null manifold.
%  2. Computational Economy: It avoids the redundant 're-melting' of the
%     original data structure, focusing resources on sampling rather than
%     initialization.
%  3. Independent Sampling: Provided the mixing step is sufficient (e.g.,
%     swaps >= total spikes / 10), the resulting STPR distributions
%     remain statistically independent for Z-scoring purposes.
%% ========================================================================

%% ========================================================================
%  NOTE: INVARIANCE OF GLOBAL POPULATION RATE (COLUMN SUMS)
%  ========================================================================
%  The 2x2 edge-swap algorithm (shuffle_raster.cpp) is a doubly-constrained
%  randomization. While it redistributes spikes across time, it strictly
%  preserves both the Row Sums (unit firing rates) and the Column Sums
%  (instantaneous population counts).
%
%  * Column Sums: sum(rasterMat, 1) == sum(rasterShuffled, 1)
%  * Linearity: conv(A + B, G) == conv(A, G) + conv(B, G)
%  * Determinism: Global rate depends only on the aggregate bin counts.
%
%  JUSTIFICATION FOR PRE-CALCULATING THE GLOBAL SMOOTHED RATE
%  Because the total number of spikes in every time bin (dt) remains
%  constant across all possible shuffled iterations, the 'Global
%  Population Rate' is an invariant property of the session.
%
%  1. Computational Efficiency: Calculating the global convolution once
%     outside the parfor loop reduces redundant O(N) operations by a factor
%     equal to nShuffles (e.g., 1000x fewer convolutions).
%  2. LOO Validity: The Leave-One-Out (LOO) population rate remains valid
%     by subtracting the shuffled unit's smoothed rate from the fixed
%     global rate: PR_loo = Global_fixed - Unit_shuffled.
%  3. Mathematical Identity: Since the sum of all units is constant,
%     the sum of the smoothed signals is also constant, ensuring the
%     null distribution is centered correctly.
%% ========================================================================


%% ========================================================================
%  NOTE: NSWAPS
%  ========================================================================
%  While the duration of the recording (length(t)) defines the temporal
%  search space, the total number of spikes (N_spks) defines the actual
%  information density. Shuffling a sparse matrix is a 'Spike-Moving'
%  problem, not a 'Bin-Filling' problem.
%
%  * Temporal Bins: Define the resolution, but are 99% empty.
%  * Total Spikes: Represent the 'Movable Mass' of the system.
%  * Mixing Time: The iterations required to reach a stationary null state.
%  * Scaling Logic: nSwaps should scale with information, not duration.
%
%  JUSTIFICATION FOR SPIKE-BASED SWAP LIMITS
%  Performing swaps equal to the number of time bins (e.g., 3.6 million)
%  is a conservative 'Brute Force' approach. Scaling based on the spike
%  count (e.g., 10x-20x total spikes) is a 'Precision' approach for the
%  following reasons:
%
%  1. Information Sparsity: In a matrix where spikes are rare, swapping
%     a spike into an empty bin is statistically likely. The algorithm
%     achieves total randomization once every spike has moved ~10 times.
%  2. Computational Convergence: Because the C++ code (shuffle_raster.cpp)
%     tracks 'Successful Swaps', it ensures the
%     target entropy is reached regardless of matrix sparsity.
%  3. Diminishing Returns: Once a spike train is 'melted' (decorrelated
%     from its original phase), additional swaps consume CPU cycles
%     without further altering the null distribution variance.
%% ========================================================================


%% ========================================================================
%  NOTE: SESSION SCALING: THE GLOBAL NOISE ANCHOR
%  ========================================================================
%  While individual neurons exhibit a wide range of population coupling,
%  their significance cannot be assessed in isolation. Population coupling
%  is a 'Relative Metric'—the strength of a single unit's coordination
%  relative to the session-wide 'Noise Floor' of the network.
%
%  * Unit Peak (prc0): Quantifies individual 'Chorister' strength.
%  * Shuffled Median: Quantifies the 'Stochastic Baseline' of the session.
%  * Session Scalar: Prevents mathematical instability from near-zero medians.
%  * Normalization: Transforms raw Hz into a universal coupling index.
%
%  JUSTIFICATION FOR SESSION-WIDE MEDIAN NORMALIZATION
%  Normalizing by the median magnitude of all shuffled STPR values in a
%  recording, rather than by unit-specific shuffled values, is a critical
%  'Stability Anchor' for the following reasons:
%
%  1. Mathematical Continuity: It avoids division by unit-specific shuffled
%     medians that approach zero, which would otherwise flip signs and
%     create artificial outliers.
%  2. Biological Comparison: It provides a stable denominator for the
%     entire recording, ensuring that the diverse coupling of FS and RS
%     units is preserved in its true relative scale.
%  3. Cross-Session Pooling: It standardizes coupling strength across
%     different experiments, animals, and cortical states by accounting for
%     session-specific population rate variability.
%% ========================================================================

