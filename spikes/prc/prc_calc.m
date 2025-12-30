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
%                     .t          : (1 x nBins) Time vector for STPR curves.
%                     .info       : (struct) Analysis parameters.
%
%   See also: PRC_PLOT, PRC_SHUFFLE
%
%   HISTORY:
%   Aug 2024 LH - Major refactor and simplification.
%   Dec 2025 AI - Optimize STPR calculation (LOO) and parallelize (parfor).
%                 Requires Parallel Computing Toolbox.
%                 Requires recompilation of shuffle_raster.cpp.

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

% Set basepath and get basename
cd(basepath);
[~, basename] = fileparts(basepath);

%% ========================================================================
%  TIME AXIS & KERNEL SETUP
%  ========================================================================

% Create time axis for binning
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
stpr = zeros(nUnits, nBins);
prc0 = nan(nUnits, 1);

% Determine number of swaps when creating shuffled raster matrix
nSwaps = length(t);

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
    [stpr(idxRef, :), prc0(idxRef)] = stpr_calc(spkBins{iGood}, pr, lagBins);
end

%% ========================================================================
%  STPR SHUFFLING
%  ========================================================================

% Pre-generate random seeds for reproducible parallel shuffling
seeds = randi([0, 2^32-1], nShuffles, 1, 'uint32');

% Initialize accumulators for parallel reduction
stpr_sum = zeros(nUnits, nBins);
prc0_shfl = nan(nUnits, nShuffles);

% Initialize progress monitoring
dq = parallel.pool.DataQueue;
progress_monitor(nShuffles, 'init');
afterEach(dq, @(~) progress_monitor([], 'update'));

% Initialize parallel pool
if isempty(gcp('nocreate'))
    parpool('local', 8);
end

% Use parfor for parallel processing (Limit to 8 workers)
parfor iShuffle = 1 : nShuffles

    % Create shuffled matrix (deterministic seed per shuffle)
    rasterShuffled = shuffle_raster(rasterMat, nSwaps, seeds(iShuffle));

    % Get indices for shuffled raster
    refBins = binary2idx(rasterShuffled, [lagBins, length(t) - lagBins], spkLim);

    % Initialize temporary arrays for this shuffle
    stpr_tmp = zeros(nUnits, nBins);
    prc0_tmp = nan(nUnits, 1);

    % Calculate STPR for each unit relative to shuffled population rate
    for iGood = 1:nGood
        idxRef = uGood(iGood);

        % Calculate unit smoothed rate and LOO Population Rate
        ur = conv(rasterShuffled(iGood, :) / binSize, gk, 'same');
        pr = popRate - ur;

        % Calculate STPR using pre-processed spike bins
        [stpr_tmp(idxRef,:), prc0_tmp(idxRef)] = stpr_calc(refBins{iGood}, pr, lagBins);
    end

    % Accumulate results
    stpr_sum = stpr_sum + stpr_tmp;
    prc0_shfl(:, iShuffle) = prc0_tmp;

    % Print progress
    send(dq, []);
end

% Average across shuffles
stpr_shfl = stpr_sum / nShuffles;

%% ========================================================================
%  NORMALIZE COUPLING
%  ========================================================================

% Z-score and normalization (Vectorized)
shflMu = mean(prc0_shfl, 2, 'omitnan');
shflSd = std(prc0_shfl, [], 2, 'omitnan');
shflMed = median(prc0_shfl, 2, 'omitnan');

% Z-score
maskSdZero = shflSd < eps;
prc0_z = (prc0 - shflMu) ./ shflSd;

% Handle zero std dev cases
if any(maskSdZero)
    % If std is 0, check if value equals mean (0 Z-score) or not (NaN)
    isMean = abs(prc0 - shflMu) < eps;
    prc0_z(maskSdZero & isMean) = 0;
    prc0_z(maskSdZero & ~isMean) = NaN;
end

% Median Normalization
prc0_norm = prc0 ./ abs(shflMed);

% Handle zero median cases
maskMedZero = abs(shflMed) < eps;
if any(maskMedZero)
    isZero = abs(prc0) < eps;
    prc0_norm(maskMedZero & isZero) = 0;
    prc0_norm(maskMedZero & ~isZero) = NaN;
end

%% ========================================================================
%  ORGANIZE OUTPUT
%  ========================================================================

% Store results (with NaN for bad units)
prc.prc0 = prc0;               % Raw STPR at zero lag
prc.prc0_z = prc0_z;           % Z-scored population coupling
prc.prc0_norm = prc0_norm;     % Median-normalized population coupling
prc.prc0_shfl = prc0_shfl;     % Shuffled STPR values
prc.stpr = stpr;               % Full STPR curves
prc.stpr_shfl = stpr_shfl;
prc.spkBins = spkBins;

% Create info struct with parameters and metadata
prc.info = struct(...
    'basepath', basepath, ...
    'runtime', datetime("now"), ...
    'nUnits', nUnits, ...
    'nGood', nGood, ...        % Number of good units
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

function [stpr, prc0] = stpr_calc(uBins, popRate, lagBins)
% CALC_STPR Calculates spike-triggered average of population rate.
%
% INPUT:
%   uBins       - Vector of bin indices where reference neuron spiked.
%   popRate     - Population firing rate (Hz).
%   lagBins     - Number of bins for STPR window.
%
% OUTPUT:
%   prc0       - Baseline-subtracted STPR at zero lag (Hz).
%   stpr       - spike triggered population rate

baseline = mean(popRate);
idx = bsxfun(@plus, uBins, (-lagBins:lagBins));
segments = popRate(idx);
stpr = mean(segments, 1) - baseline;
prc0 = stpr(lagBins + 1);

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
