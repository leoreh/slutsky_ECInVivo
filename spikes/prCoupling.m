function prc = prCoupling(spktimes, varargin)
% PRCOUPLING Calculates population coupling for all neurons in spktimes.
%
% SUMMARY:
% This function implements a hybrid approach to calculate population coupling:
%   1. Core metric: baseline-subtracted spike-triggered population rate (stPR)
%      at zero lag, following Okun et al. (2015). 
%   2. Shuffling: Two methods available:
%      a) Circular temporal shifts of individual population neuron spike trains
%         (Dorman et al., 2023)
%      b) Raster marginals model preserving row and column sums (Okun et al., 2015)
%   3. Normalization: Two methods available:
%      a) Z-score against shuffled distribution (Dorman et al., 2023)
%      b) Median normalization (Okun et al., 2015)
%
% INPUT (Required):
%   spktimes      - Cell array. spktimes{i} contains spike times (s) for neuron i.
%
% INPUT (Optional Key-Value Pairs):
%   basepath      - Path to recording session directory {pwd}.
%   flgSave       - Logical flag to save results to .prc.mat file {true}.
%   winLim        - [start end] time window to analyze [s] {[0 Inf]}.
%   binSize       - Bin size for spike trains [s] {0.001}.
%   gkHw          - Half-width (sigma) of Gaussian kernel for rate smoothing [s] {0.012}.
%   winStpr       - Total width of STPR window [s] {1.0}.
%   nShuffles     - Number of shuffles for null distribution {1000}.
%   spkThr        - Minimum number of spikes required per unit {300}. Units with
%                   fewer spikes will have NaN values in all output fields.
%   flg_par       - Logical flag to use parfor for shuffling {false}.
%   spkLim        - Maximum number of spikes to use from reference unit {Inf}.
%                   If a unit has more spikes, a random subset is used to
%                   reduce computation time while maintaining statistical power.
%   shuffleMet    - Method for shuffling spike trains {'raster'}. Options:
%                   'circshift': circular temporal shifts (Dorman et al., 2023)
%                   'raster': raster marginals model (Okun et al., 2015)
%
% OUTPUT:
%   prc           - Structure containing population coupling results:
%     .prc0_z     - Z-scored population coupling for each neuron. NaN for units
%                   with fewer than spkThr spikes.
%     .prc0_norm  - Median-normalized population coupling (Okun et al. 2015).
%                   NaN for units with fewer than spkThr spikes.
%     .prc0       - Raw, baseline-subtracted STPR at zero lag. NaN for bad units.
%     .stpr       - Full STPR curves for each neuron. NaN for bad units.
%     .sptr0shfl  - Distribution of STPR(0) from shuffles. NaN for bad units.
%     .t          - Time vector for STPR curves [s].
%     .info       - Analysis parameters and metadata.
% 
% NOTES:
%   - Units with fewer than spkThr spikes are excluded from analysis and
%     have NaN values in all output fields.
%   - For each reference unit, we calculate the population rate by summing
%     all other units' spike trains, then convolve with a Gaussian kernel.
%     This ensures the reference unit doesn't contribute to its own STPR
%     while maintaining correct population coupling calculation.
%
% DEPENDENCIES:
%   None
%
% HISTORY:
%   Aug 2024 LH - Major refactor and simplification.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING & INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

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
addParameter(p, 'flg_par', false, @islogical);
addParameter(p, 'spkLim', Inf, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'shuffleMet', 'raster', @(x) ismember(x, {'raster', 'circshift'}));

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
flg_par = p.Results.flg_par;
spkLim = p.Results.spkLim;
shuffleMet = p.Results.shuffleMet;

% Set basepath and get basename
cd(basepath);
[~, basename] = fileparts(basepath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME AXIS & KERNEL SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPIKE WINDOWING & INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Window spikes to analysis window using cellfun
spktimes = cellfun(@(x) x(x >= winLim(1) & x <= winLim(2)), spktimes, 'UniformOutput', false);
nSpks = cellfun(@(x) length(x), spktimes, 'UniformOutput', true);
uGood = find(nSpks > spkThr);
nGood = length(uGood);
nUnits = length(spktimes);

% Initialize output arrays at full size
stpr = zeros(nUnits, nBins);
prc0_z = nan(nUnits, 1);      % Renamed from z0
prc0 = nan(nUnits, 1);
prc0_shfl = nan(nUnits, nShuffles);
stpr_shfl = zeros(nUnits, nBins);
prc0_norm = nan(nUnits, 1);   % Median-normalized coupling

% Generate all random shifts once. popShifts(i,j) is the shift for the i-th
% good unit in the j-th shuffl
popShifts = randi(length(t), nGood, nShuffles);

% Determine number of swaps when creating shuffled raster matrix
nSwaps = length(t) * nUnits;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARE RASTER MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize binRates only for good units
binRates = zeros(nGood, length(t));  % Only store good units
spkBins = cell(nGood, 1);              

% Calculate binned rates across good units and store spike bin indices
for iGood = 1:nGood
    idxRef = uGood(iGood);
    uSpks = spktimes{idxRef};
    binRates(iGood,:) = histcounts(uSpks, [t, t(end) + binSize]);
    spkBins{iGood} = find(binRates(iGood,:) > 0);
    
    % Subsample spike bins if needed
    rng(idxRef);                                    % Use consistent random seed for reproducibility
    if length(spkBins{iGood}) > spkLim        
        spkBins{iGood} = spkBins{iGood}(randperm(length(spkBins{iGood}), spkLim));
        spkBins{iGood} = sort(spkBins{iGood});      % Keep bins sorted in time
    end
end

% Convert to binary
rasterMat = double(binRates > 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALC STPR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iGood = 1:nGood
    idxRef = uGood(iGood);  
    unitMask = true(nGood, 1);
    unitMask(iGood) = false;

    % Calculate population rate with leave-one-out and convolve
    pr = sum(rasterMat(unitMask), 1, 'omitnan') / binSize;
    pr = conv(pr, gk, 'same');

    % Calculate STPR using pre-processed spike bins
    [stpr(idxRef, :), prc0(idxRef)] = stpr_calc(spkBins{iGood}, pr, lagBins);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALC STPR TO SHUFFLED POPULATION RATES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pre-allocate cell arrays for parallel results
stpr_cell = cell(nShuffles, 1);
prc0_cell = cell(nShuffles, 1);

% Initialize progress indicator
fprintf('Shuffling progress: 0/%d', nShuffles);

% Parallel processing
for iShuffle = 1 : nShuffles
    % Create shuffled matrix using specified method
    if strcmp(shuffleMet, 'raster')
        rasterShuffled = shuffle_raster(rasterMat, nSwaps);
    
    elseif strcmp(shuffleMet, 'circshift')
        rasterShuffled = zeros(size(rasterMat));
        for iGood = 1:nGood
            shift = popShifts(iGood, iShuffle);
            rasterShuffled(iGood, :) = circshift(binRates(iGood,:), [0, shift]);
        end
    end
    
    % Initialize temporary arrays for this shuffle
    stpr_tmp = zeros(nUnits, nBins);
    prc0_tmp = zeros(nUnits, 1);
    
    % Calculate STPR for each unit relative to shuffled population rate
    for iGood = 1:nGood
        idxRef = uGood(iGood);
        unitMask = true(nGood, 1);
        unitMask(iGood) = false;

        % Calculate shuffled population rate with leave-one-out
        prShuffeld = sum(rasterShuffled(unitMask, :), 1, 'omitnan') / binSize;
        prShuffeld = conv(prShuffeld, gk, 'same');

        % Calculate STPR using pre-processed spike bins
        [stpr_tmp(idxRef,:), prc0_tmp(idxRef)] = stpr_calc(spkBins{iGood}, prShuffeld, lagBins);
    end
    
    % Store results in cell arrays
    stpr_cell{iShuffle} = stpr_tmp;
    prc0_cell{iShuffle} = prc0_tmp;
    
    % Print progress every 10 shuffles
    if mod(iShuffle, 1) == 0
        fprintf('\rShuffling progress: %d/%d', iShuffle, nShuffles);
    end
end
fprintf('\n'); % New line after completion

% Aggregate results
stpr_shfl = cell2padmat(stpr_cell, 3, 0);      % Concatenate along 3rd dim
stpr_shfl = mean(stpr_shfl, 3);                % Average across shuffles
[prc0_shfl, ~] = cell2padmat(prc0_cell, 2, 0); % Concatenate along 2nd dim

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NORMALIZE COUPLING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iGood = 1:nGood
    refIdx = uGood(iGood);

    % Z-score normalization (Dorman et al., 2023)
    shflMu = mean(prc0_shfl(refIdx,:), 2, 'omitnan');
    shflSd = std(prc0_shfl(refIdx,:), [], 2, 'omitnan');
    shflMed = median(prc0_shfl(refIdx,:), 2, 'omitnan');

    % Calculate Z-score (Dorman et al., 2023)
    if shflSd < eps
        if abs(prc0(refIdx) - shflMu) < eps
            prc0_z(refIdx) = 0;
        else
            prc0_z(refIdx) = NaN;
        end
    else
        prc0_z(refIdx) = (prc0(refIdx) - shflMu) / shflSd;
    end

    % Calculate median normalization (Okun et al., 2015)
    if abs(shflMed) < eps
        if abs(prc0(refIdx)) < eps
            prc0_norm(refIdx) = 0;
        else
            prc0_norm(refIdx) = NaN;
        end
    else
        prc0_norm(refIdx) = prc0(refIdx) / abs(shflMed);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORGANIZE OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    'nShuffles', nShuffles, ...
    'shuffleMet', shuffleMet); % Shuffling method used

% Save results if requested
if flgSave
    save(fullfile(basepath, [basename, '.prc.mat']), 'prc', '-v7.3');
end

end     % EOF


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION: CALCULATE STPR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

nBins = 2 * lagBins + 1;
segments = zeros(length(uBins), nBins);
validSpks = 0;

for iSpk = 1:length(uBins)
    spikeBin = uBins(iSpk);
    startBin = spikeBin - lagBins;
    endBin = spikeBin + lagBins;
    
    if startBin >= 1 && endBin <= length(popRate)
        segments(iSpk,:) = popRate(startBin:endBin);
        validSpks = validSpks + 1;
    end
end

% Calculate baseline and remove from stpr
baseline = mean(popRate);

stpr = mean(segments(1:validSpks,:), 1) - baseline;
prc0 = stpr(lagBins + 1);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION: SHUFFLE RASTER MATRIX (replaced with mex c++ version)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function shuffledMat = shuffle_raster(rasterMat, nSwaps)
% % SHUFFLE_RASTER Shuffles binary spike matrix while preserving marginal sums.
% %
% % INPUT:
% %   rasterMat     - Binary matrix (neurons x time_bins) of spike times.
% %   nSwaps        - Number of swap attempts. Defaults to nBins if empty.
% %
% % OUTPUT:
% %   shuffledMat   - Shuffled spike matrix with preserved marginal sums.
% %
% % NOTES:
% %   Implements the "raster marginals model" from Okun et al. (2015) by swapping
% %   2x2 submatrices that preserve both row (neuron rates) and column (population
% %   activity) sums. Specifically swaps patterns [[1,0],[0,1]] ↔ [[0,1],[1,0]].
% 
% % Get matrix dimensions and set default swaps
% [nUnits, nBins] = size(rasterMat);
% if nargin < 2 || isempty(nSwaps)
%     nSwaps = nBins * nUnits;
% end
% 
% % Initialize output
% shuffledMat = rasterMat;
% 
% % Pre-generate all random indices for rows and columns
% r1_all = randi(nUnits, nSwaps, 1);
% r2_all = randi(nUnits-1, nSwaps, 1);
% r2_all = r2_all + (r2_all >= r1_all);  % Ensure f2 != f1
% c1_all = randi(nBins-1, nSwaps, 1);
% c2_all = randi(nBins-1, nSwaps, 1);
% c2_all = c2_all + (c2_all >= c1_all);  % Ensure t2 != t1
% 
% % Perform swaps
% for iSwap = 1:nSwaps
%     % Get current indices
%     r1 = r1_all(iSwap);
%     r2 = r2_all(iSwap);
%     c1 = c1_all(iSwap);
%     c2 = c2_all(iSwap);
% 
%     % Get values directly
%     v11 = shuffledMat(r1, c1);
%     v12 = shuffledMat(r1, c2);
%     v21 = shuffledMat(r2, c1);
%     v22 = shuffledMat(r2, c2);
% 
%     % Swap if pattern matches [1,0;0,1] or [0,1;1,0]
%     if (v11 == 1 && v12 == 0 && v21 == 0 && v22 == 1)
%         shuffledMat(r1, c1) = 0;
%         shuffledMat(r1, c2) = 1;
%         shuffledMat(r2, c1) = 1;
%         shuffledMat(r2, c2) = 0;
%     elseif (v11 == 0 && v12 == 1 && v21 == 1 && v22 == 0)
%         shuffledMat(r1, c1) = 1;
%         shuffledMat(r1, c2) = 0;
%         shuffledMat(r2, c1) = 0;
%         shuffledMat(r2, c2) = 1;
%     end
% end
% 
% end     % EOF



