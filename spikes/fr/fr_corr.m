function cc = fr_corr(spkmat, varargin)
% FR_CORR Computes pairwise correlations of single-unit FR trajectories.
%
%   cc = FR_CORR(spkmat, ...) calculates functional connectivity metrics
%   between all pairs of single-unit FR trajectories. It computes TWO
%   complementary metrics simultaneously:
%
%   1. SHUFFLE Z-SCORE (Reliability):
%       - Z-scores the raw correlation using a per-pair noise distribution
%         derived from circular shuffling.
%       - Metric: Standard Deviations above chance (Z).
%       - Good for: Filtering out burst-driven artifacts.
%
%   2. FISHER Z (Magnitude):
%       - Uses a Global Noise Threshold (Mean + K*Std).
%       - Masks correlations below threshold.
%       - Metric: Pearson Correlation (r).
%       - Good for: Quantifying coupling strength of established links.
%
%   INPUTS:
%       spkmat      - (matrix) Activity Matrix (Neurons x Time).
%       varargin    - (param/value) Optional parameters:
%                     'nShuffles'     : (int) Number of shuffles for noise est {10}
%                     'thrSig'        : (num) SD multiplier for global threshold {2}
%                     'flgPlot'       : (log) Plot correlation matrix {false}
%
%   OUTPUTS:
%       cc          - (struct) Result structure with fields:
%                     .raw      : Raw Pearson correlations (diag removed).
%                     .noise    : Noise statistics (mean, std, limit).
%                     .shuffle  : Shuffle Z-score metrics (zcc, mcc, funcon).
%                     .fisher   : Fisher Masked metrics (r, mcc, funcon, mask).
%                     .params   : Input parameters.
%
%   See also DIM_CALC, FR_NETWORK

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'spkmat', @(x) isnumeric(x) || islogical(x));
addParameter(p, 'nShuffles', 10, @isnumeric);
addParameter(p, 'thrSig', 2, @isnumeric);
addParameter(p, 'flgPlot', false, @islogical);

parse(p, spkmat, varargin{:});
spkmat       = double(p.Results.spkmat);
nShuffles    = p.Results.nShuffles;
thrSig       = p.Results.thrSig;
flgPlot      = p.Results.flgPlot;


%% ========================================================================
%  PRE-PROCESSING
%  ========================================================================

% Initialize Output
cc.raw      = struct('mcc', NaN, 'funcon', [], 'ccmat', []);
cc.noise    = struct('mean', NaN, 'std', NaN, 'limit', NaN);
cc.shuffle  = struct('mcc', NaN, 'funcon', [], 'ccmat', []);
cc.fisher   = struct('mcc', NaN, 'funcon', [], 'ccmat', [], 'mask', []);
cc.params   = p.Results;

% Remove silent neurons
nUnitsRaw = size(spkmat, 1);
mY = mean(spkmat, 2);
validIdx = mY > 0;
nUnits = sum(validIdx);

% Z-score valid units (allows using cov() for Pearson r)
spkz = zscore(spkmat(validIdx, :), 0, 2);

%% ========================================================================
%  RAW CORRELATIONS
%  ========================================================================

% cov expects observations in rows, variables in cols.
% SPKZ is Neurons x Time. We want Neuron-Neuron corr.
r = cov(spkz');

% Remove diagonal (auto-correlations)
r = r - diag(diag(r));

% Raw Statistics
rOffDiag = r(eye(nUnits) == 0);

cc.raw.mcc = mean(rOffDiag);
funcon = sum(abs(r), 2) ./ (nUnits - 1);

% Expand Matrix Maps
cc.raw.ccmat = nan(nUnitsRaw);
cc.raw.ccmat(validIdx, validIdx) = r;
cc.raw.funcon = nan(nUnitsRaw, 1);
cc.raw.funcon(validIdx) = funcon;


%% ========================================================================
%  NOISE ESTIMATION
%  ========================================================================

sumR = zeros(nUnits);
sumR2 = zeros(nUnits);
allNoiseVals = [];

for iShuff = 1:nShuffles
    Z_shuff = zeros(size(spkz));
    nTime = size(spkz, 2);

    % Circular shift each neuron independently
    for iUnit = 1:nUnits
        shift = randi(nTime);
        Z_shuff(iUnit, :) = circshift(spkz(iUnit, :), shift, 2);
    end

    % Compute shuffle correlations
    rShuff = cov(Z_shuff');
    rShuff = rShuff - diag(diag(rShuff));

    % Accumulate
    sumR = sumR + rShuff;
    sumR2 = sumR2 + rShuff.^2;

    % Store off-diagonal for global noise distribution
    allNoiseVals = [allNoiseVals; rShuff(eye(nUnits) == 0)];
end

% Per-Pair Noise Stats
noiseAvg = sumR / nShuffles;
noiseVar = (sumR2 - (sumR.^2)/nShuffles) / (nShuffles - 1);
noiseVar = max(noiseVar, 0);
noiseStd = sqrt(noiseVar);

% Regularization
% Clamp standard deviations to the "physical floor" of noise variance in
% this network. Prevents Z-score explosion due to under-sampling
% noise. Should not occur with high nShuffles (> 500).
noiseFloor = prctile(noiseStd(eye(nUnits)==0), 5);
noiseStd(noiseStd < noiseFloor) = noiseFloor;

cc.noise.mean  = noiseAvg;
cc.noise.std   = noiseStd;
cc.noise.limit = mean(allNoiseVals) + std(allNoiseVals) * thrSig;   % Global noise

%% ========================================================================
%  METRIC 1: SHUFFLE Z-SCORE
%  ========================================================================
%  Z = (r - mu_pair) / sigma_pair

z = (r - noiseAvg) ./ noiseStd;
z = z - diag(diag(z));

% Save Results
maskOffDiag = true(nUnits) & ~eye(nUnits);
cc.shuffle.mcc = mean(z(maskOffDiag));
cc.shuffle.funcon = sum(abs(z), 2) ./ (nUnits - 1);

% Expand
cc.shuffle.ccmat = nan(nUnitsRaw);
cc.shuffle.ccmat(validIdx, validIdx) = z;


%% ========================================================================
%  METRIC 2: FISHER MASKED
%  ========================================================================
% KEEP MAGNITUDE (r), zero out insignificant

% maskSig = abs(r) > cc.noise.limit;
maskSig = abs(z) > thrSig;
rSig = r .* maskSig; 

% Compute Mean
valsSig = r(maskSig);
if isempty(valsSig)
    mccFisher = 0;
else
    % Fisher Transform -> Mean -> Inverse
    valsSig = max(min(valsSig, 0.999), -0.999);
    zVals = atanh(valsSig);
    mccFisher = tanh(mean(zVals));
end

% Save Results
cc.fisher.mcc = mccFisher;
cc.fisher.funcon = sum(abs(rSig), 2) ./ (nUnits - 1);

% Expand
cc.fisher.ccmat = nan(nUnitsRaw);
cc.fisher.ccmat(validIdx, validIdx) = rSig;
cc.fisher.mask = false(nUnitsRaw);
cc.fisher.mask(validIdx, validIdx) = maskSig;


%% ========================================================================
%  PLOT
%  ========================================================================

if flgPlot
    figure('Color', 'w', 'Position', [100 100 1400 500]);

    % Raw + Noise Dist
    subplot(1, 4, 1);
    hold on;
    histogram(allNoiseVals, 50, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'DisplayName', 'Noise', 'LineWidth', 1.5);
    histogram(rOffDiag, 50, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'DisplayName', 'Data', 'LineWidth', 1.5);
    xline(cc.noise.limit, 'r--', 'Thr'); xline(-cc.noise.limit, 'r--', 'Thr');
    legend('Location', 'best');
    axis square;
    title(['Distributions (Thr: ' num2str(thrSig) '\sigma)']);
    xlabel('CV (Pearson r)');

    % Raw Matrix
    subplot(1, 4, 2);
    imagesc(cc.raw.ccmat, [-0.5 0.5]);
    title(['Raw (mcc: ' num2str(cc.raw.mcc, '%.2f') ')']);
    axis square; colorbar;

    % Shuffle Z Matrix
    subplot(1, 4, 3);
    imagesc(cc.shuffle.ccmat, [-5 10]); % Z-score range
    title(['Shuffle Z (mcc: ' num2str(cc.shuffle.mcc, '%.2f') ')']);
    axis square; colorbar;

    % Fisher Masked Matrix
    subplot(1, 4, 4);
    imagesc(cc.fisher.ccmat, [-0.5 0.5]);
    title(['Fisher Z (mcc: ' num2str(cc.fisher.mcc, '%.2f') ')']);
    axis square; colorbar;

    sgtitle('Functional Connectivity Analysis');
end

end     % EOF
