function cc = fr_corr(Y, varargin)
% FR_CORR Computes pairwise correlations of single-unit FR trajectories.
%
%   cc = FR_CORR(Y) calculates the Pearson correlation coefficients between
%   all pairs of single-unit FR trajectories. It employs a shuffling
%   procedure to determine a significance threshold (noise floor) and
%   computes metrics based on both valid (significant) and raw correlations.
%
%   METHODOLOGY:
%   1.  Omit silent neurons.
%   2.  Z-score each unit (temporally).
%   3.  Compute raw covariance matrix (Pearson correlation).
%   4.  Estimate noise floor via circular shuffling (shift-predictor):
%       - Each unit's spike train is circularly shifted by a random amount.
%       - Correlations are re-computed on shifted data.
%       - This preserves auto-correlation but destroys cross-correlation.
%   5.  Define threshold = mean(Noise) + K * std(Noise).
%   6.  Identify "significant" correlations as those where |raw| > threshold.
%   7.  Compute average correlation from these significant pairs.
%
%   INPUTS:
%       Y           - (matrix) Activity Matrix (Neurons x Time).
%       varargin    - (param/value) Optional parameters:
%                     'nShuffles'     : (int) Number of shuffles for noise est {10}
%                     'thresholdVal'  : (num) SD multiplier for threshold {2}
%                     'flgPlot'       : (log) Plot correlation matrix {false}
%
%   OUTPUTS:
%       cc          - (struct) Result structure:
%                     .mcc      : Mean Signed correlation of significant pairs.
%                     .mccAbs   : Mean Absolute correlation of significant pairs.
%                     .mccRaw   : Mean correlation of all pairs (diag removed).
%                     .cc       : Raw correlation matrix (diag removed).
%                     .ccSig    : significant correlation matrix (others 0).
%                     .mask     : Logical mask of significant pairs.
%                     .noise    : Noise statistics (mean, std, limit).
%                     .params   : Input parameters.
%
%   See also DIM_CALC, FR_NETWORK

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'Y', @isnumeric);
addParameter(p, 'nShuffles', 10, @isnumeric);
addParameter(p, 'thresholdVal', 2, @isnumeric);
addParameter(p, 'flgPlot', false, @islogical);
addParameter(p, 'flgLog', true, @islogical);

parse(p, Y, varargin{:});
nShuffles    = p.Results.nShuffles;
thresholdVal = p.Results.thresholdVal;
flgPlot      = p.Results.flgPlot;
flgLog       = p.Results.flgLog;


%% ========================================================================
%  PRE-PROCESSING
%  ========================================================================

% Initialize Output
cc.mcc    = NaN;
cc.mccAbs = NaN;
cc.mccRaw = NaN;
cc.cc     = [];
cc.ccSig  = [];
cc.ccNorm = []; % Degree-Normalized Matrix
cc.mask   = [];
cc.noise  = struct('mean', NaN, 'std', NaN, 'limit', NaN);
cc.params = p.Results;

% Remove silent or constant neurons (prevent NaNs in zscore)
nUnitsRaw = size(Y, 1);
mY = mean(Y, 2);
sY = std(Y, 0, 2);
validIdx = mY > 0 & sY > eps; % Must be active AND vary
nUnits = sum(validIdx);

if nUnits < 2
    cc.cc = nan(nUnitsRaw);
    return;
end

% Z-score only valid units
% If flgLog is true, apply log transformation first to handle lognormal rates
Y_valid = Y(validIdx, :);
if flgLog
    Y_valid = log(Y_valid + eps);
end

% Normalizing allows using cov() to get Pearson correlation directly
Z = zscore(Y_valid, 0, 2);


%% ========================================================================
%  RAW CORRELATIONS
%  ========================================================================

% cov expects observations in rows, variables in cols.
% Z is Neurons x Time. We want Neuron-Neuron corr.
% So we treat Time as observations (rows) and Neurons as variables (cols).
ccValid = cov(Z');

% Remove diagonal (auto-correlations)
ccValid = ccValid - diag(diag(ccValid));

% Raw Mean Correlation (average of all off-diagonal elements)
% Since matrix is symmetric and 0 on diag, sum(all)/(N^2-N) is mean.
valOffDiag = ccValid(eye(nUnits) == 0);
cc.mccRaw = mean(valOffDiag);


%% ========================================================================
%  NOISE ESTIMATION
%  ========================================================================

valsShuffle = [];

for iShuff = 1:nShuffles
    Z_shuff = zeros(size(Z));
    nTime = size(Z, 2);

    % Circular shift each neuron independently
    for iU = 1:nUnits
        shift = randi(nTime);
        Z_shuff(iU, :) = circshift(Z(iU, :), shift, 2);
    end

    % Compute shuffle correlations
    ccShuff = cov(Z_shuff');
    ccShuff = ccShuff - diag(diag(ccShuff));

    % Collect off-diagonal values
    valsShuffle = [valsShuffle; ccShuff(eye(nUnits) == 0)]; %#ok<AGROW>
end

% Noise Statistics
noiseMean = mean(valsShuffle);
noiseStd  = std(valsShuffle);
limit = noiseStd * thresholdVal;

cc.noise.mean  = noiseMean;
cc.noise.std   = noiseStd;
cc.noise.limit = limit;

%% ========================================================================
%  CALC: SIGNIFICANT CORRELATIONS
%  ========================================================================

% Identify significant pairs
% Using absolute value for thresholding (strong positive OR negative corr)
maskValid = abs(ccValid) > limit;

% Apply mask (keep sign of original correlation)
ccSigValid = ccValid;
ccSigValid(~maskValid) = 0;

% Compute Mean of significant correlations
valsSig = ccValid(maskValid);

if isempty(valsSig)
    cc.mcc    = 0;
    cc.mccAbs = 0;
else
    % Fisher Z-Transformation for correct averaging
    % Clamp values to prevent Inf at +/- 1
    valsClamped = max(min(valsSig, 1-eps), -1+eps);
    zVals = atanh(valsClamped);
    cc.mcc = tanh(mean(zVals));

    % For Absolute mean, we take mean of abs values directly or handle Z of abs?
    % Standard practice for "Magnitude of Synchrony" often just averages Abs(r).
    % However, to be rigorous, one could average Z(|r|).
    % Let's stick to simple mean(abs) for mccAbs as it's a descriptive stat,
    % but the mcc (signed) is the critical one for bias.
    cc.mccAbs = mean(abs(valsSig));
end

% Degree Normalization (Geometric Mean Normalization)
% Helps remove "Rich Club" artifacts driven by high firing rates.
% Formula: F_ij = W_ij / sqrt(d_i * d_j), where d is weighted degree (strength).

W_abs = abs(ccSigValid);
d = sum(W_abs, 2); % Weighted degree (Strength)
d(d==0) = 1;       % Avoid division by zero for isolated nodes

[Di, Dj] = meshgrid(d, d);
ccNorm = ccSigValid ./ sqrt(Di .* Dj);

% Mean Functional Connectivity per Unit (from Normalized Matrix)
% Using the normalized matrix for "funcon" to allow comparison between
% high/low rate neurons without rate bias.
funcon = sum(abs(ccNorm), 2);

%% ========================================================================
%  EXPAND
%  ========================================================================

% Reconstruct full matrices (keeping NaNs for silent units)
cc.cc = nan(nUnitsRaw);
cc.cc(validIdx, validIdx) = ccValid;

cc.ccSig = nan(nUnitsRaw);
cc.ccSig(validIdx, validIdx) = ccSigValid;

cc.ccNorm = nan(nUnitsRaw);
cc.ccNorm(validIdx, validIdx) = ccNorm;

cc.mask = false(nUnitsRaw);
cc.mask(validIdx, validIdx) = maskValid;

cc.funcon = nan(nUnitsRaw, 1);
cc.funcon(validIdx) = funcon;

%% ========================================================================
%  PLOT
%  ========================================================================

if flgPlot
    figure('Color', 'w', 'Position', [100 100 1200 400]);

    % 1. Raw Matrix
    subplot(1, 3, 1);
    imagesc(cc.cc);
    colorbar; axis square;
    title(['Raw (mcc: ' num2str(cc.mccRaw, '%.3f') ')']);
    xlabel('Neuron ID'); ylabel('Neuron ID');

    % 2. Histogram of Values vs Noise
    subplot(1, 3, 2);
    hold on;
    histogram(valsShuffle, 50, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'DisplayName', 'Noise');
    histogram(valOffDiag, 50, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'DisplayName', 'Data');
    xline(limit, 'r--', 'Limit');
    xline(-limit, 'r--', 'Limit');
    legend('Location', 'best');
    title(['Threshold: ' num2str(thresholdVal) ' \sigma']);
    axis square;

    % 3. Significant Matrix
    subplot(1, 3, 3);
    imagesc(cc.ccSig);
    colorbar; axis square;
    title(['Significant (mcc: ' num2str(cc.mcc, '%.3f') ')']);
    xlabel('Neuron ID'); ylabel('Neuron ID');

    sgtitle('Pairwise Correlation Analysis');
end

end     % EOF



