function prc = prCoupling(spktimes, varargin)
% PRCOUPLING Calculates population coupling for all neurons in spktimes.
%
% SUMMARY:
% This function implements a hybrid approach to calculate population coupling:
%   1. Core metric: baseline-subtracted spike-triggered population rate (stPR)
%      at zero lag, following Okun et al. (2015).
%   2. Shuffling: circular temporal shifts of individual population neuron
%      spike trains (Dorman et al., 2023). Reference neuron is not shuffled.
%   3. Normalization: Z-score against shuffled distribution mean and std.
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
%
% OUTPUT:
%   prc           - Structure containing population coupling results:
%     .z0         - Z-scored population coupling for each neuron. NaN for units
%                   with fewer than spkThr spikes.
%     .stpr0      - Raw, baseline-subtracted STPR at zero lag. NaN for bad units.
%     .stpr       - Full STPR curves for each neuron. NaN for bad units.
%     .sptr0shfl  - Distribution of STPR(0) from shuffles. NaN for bad units.
%     .t          - Time vector for STPR curves [s].
%     .info       - Analysis parameters and metadata.
% 
% NOTES:
%   - Units with fewer than spkThr spikes are excluded from analysis and
%     have NaN values in all output fields.
%   - The population rate is computed as the sum of all good units' rates,
%     excluding the reference unit for each calculation.
%   - Shuffling is performed by circularly shifting each unit's spike train
%     independently, preserving the reference unit's spikes.
%   - The Gaussian kernel is used to smooth both the population rate and
%     individual unit rates before computing STPR.
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
stpr = nan(nUnits, nBins);
z0 = nan(nUnits, 1);
stpr0 = nan(nUnits, 1);
sptr0shfl = nan(nUnits, nShuffles);

% Initialize binRates only for good units
binRates = zeros(nGood, length(t));  % Only store good units
spkBins = cell(nGood, 1);              

% Calculate binned rates across good units and store spike bin indices
for iGood = 1:nGood
    uIdx = uGood(iGood);
    uSpks = spktimes{uIdx};
    binRates(iGood,:) = histcounts(uSpks, [t, t(end) + binSize]);
    spkBins{iGood} = find(binRates(iGood,:) > 0);
end

% Pre-compute total population rate (all good units)
popRateAll = sum(binRates, 1, 'omitnan');  % Sum across all good units
popRateAll = conv(popRateAll, gk, 'same');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARE SHUFFLED RATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pre-compute shuffled population rates for all units
popShifts = randi(length(t), nGood, nShuffles);
popRateShfl = zeros(nShuffles, length(t));

% Pre-compute all shuffled population rates
for iShfl = 1:nShuffles
    if mod(iShfl, 50) == 0
        fprintf('Pre-computing shuffled PR %d/%d\n', iShfl, nShuffles);
    end
    
    % Shift rates and sum
    tmpRate = zeros(1, length(t));
    for iUnitPop = 1:nGood
        shift = popShifts(iUnitPop, iShfl);
        % Get shifted indices using modulo arithmetic
        shiftedInds = mod((1:length(t)) + shift - 1, length(t)) + 1;
        tmpRate = tmpRate + binRates(iUnitPop, shiftedInds);
    end
    
    % Smooth shuffled population rate
    popRateShfl(iShfl,:) = conv(tmpRate, gk, 'same');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALC STPR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iGood = 1:nGood
    uIdx = uGood(iGood);  % Get the actual unit index
    fprintf('Analyzing unit %d/%d\n', uIdx, nUnits);
    
    % Get reference unit's spike bins and subsample if needed
    uBins = spkBins{iGood};
    if length(uBins) > spkLim
        % Use consistent random seed for reproducibility. Different for
        % each unit but consistent across runs
        rng(uIdx); 
        uBins = uBins(randperm(length(uBins), spkLim));
    end
    
    % Calculate population rate by subtracting reference unit's rate
    refRate = conv(binRates(iGood,:), gk, 'same');
    popRate = popRateAll - refRate;
    
    % --- Calculate actual STPR at zero lag ---
    [stpr(uIdx, :), stpr0(uIdx)] = stpr_calc(uBins, popRate, lagBins);
    
    % --- Shuffling procedure ---
    % Initialize array for shuffled STPR values
    unitShfl = nan(1, nShuffles);

    % Calculate STPR for each pre-computed shuffled population rate
    for iShfl = 1:nShuffles
        % Get pre-computed shuffled population rate and subtract reference unit
        popRate = popRateShfl(iShfl,:) - refRate;
        
        % Calculate STPR for this shuffle
        [~, unitShfl(iShfl)] = stpr_calc(uBins, popRate, lagBins);
    end
    
    sptr0shfl(uIdx,:) = unitShfl;
    
    % --- Z-score normalization ---
    shflMu = mean(unitShfl, 'omitnan');
    shflSd = std(unitShfl, 0, 'omitnan');
    
    if shflSd < eps
        if abs(stpr0(uIdx) - shflMu) < eps
            z0(uIdx) = 0;
        else
            z0(uIdx) = NaN;
        end
    else
        z0(uIdx) = (stpr0(uIdx) - shflMu) / shflSd;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORGANIZE OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Store results (already in full size with NaN for bad units)
prc.stpr = stpr;               % Full STPR curves
prc.stpr0 = stpr0;             % Raw STPR at zero lag
prc.sptr0shfl = sptr0shfl;     % Shuffled STPR values
prc.z0 = z0;                   % Z-scored population coupling
prc.t = t;              

% Create info struct with parameters and metadata
prc.info = struct(...
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

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION: CALCULATE STPR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [stpr, stpr0] = stpr_calc(uBins, popRate, lagBins)
% CALC_STPR Calculates spike-triggered average of population rate.
%
% INPUT:
%   uBins  - Vector of bin indices where reference neuron spiked.
%   popRate     - Population firing rate (Hz).
%   lagBins     - Number of bins for STPR window.
%   nBins       - Total number of bins in STPR window.
%
% OUTPUT:
%   stpr0       - Baseline-subtracted STPR at zero lag (Hz).
%   uBins  - Indices of reference spikes in binned time.


% Calculate baseline
baseline = mean(popRate);

% Calculate STPR at zero lag
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

stpr = mean(segments(1:validSpks,:), 1);
stpr0 = stpr(lagBins + 1) - baseline;


end
