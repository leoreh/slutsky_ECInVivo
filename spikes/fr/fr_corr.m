function cc = fr_corr(spkmat, varargin)
% FR_CORR Computes pairwise correlations of single-unit FR trajectories.
%
%   cc = FR_CORR(spkmat, ...) calculates functional connectivity metrics
%   between all pairs of single-unit FR trajectories. It supports two distinct
%   normalization/thresholding methodologies: 'shuffle' (Z-score) and 'fisher'
%   (Significance Masking).
%
%   METHODOLOGY:
%   1.  Omit silent neurons.
%   2.  Z-score each unit (temporally) to allow easy Pearson calc.
%   3.  Compute raw covariance matrix (Pearson correlation).
%   4.  Estimate noise via circular shuffling (shift-predictor):
%       - Each unit's spike train is circularly shifted by a random amount.
%       - Correlations are re-computed on shifted data.
%   5.  Calculate Metrics based on 'zMet':
%       - 'shuffle' (Default):
%           * Computes mean and std of noise for EACH PAIR.
%           * Z-scores the raw correlation using its specific noise stats.
%           * Metric is the Z-score itself (Standardized Connectivity).
%       - 'fisher':
%           * Computes GLOBAL noise mean and std.
%           * Threshold = Mean + K * Std.
%           * Masks correlations where |r| < Threshold.
%           * Fisher Z-transforms significant r values, averages them, then
%             inverse transforms back to r-space.
%
%   INPUTS:
%       spkmat      - (matrix) Activity Matrix (Neurons x Time).
%       varargin    - (param/value) Optional parameters:
%                     'zMet'          : (char) 'shuffle' {default} or 'fisher'.
%                     'nShuffles'     : (int) Number of shuffles {10}.
%                     'thresholdVal'  : (num) SD multiplier for threshold {2}.
%                     'flgPlot'       : (log) Plot correlation matrix {false}.
%
%   OUTPUTS:
%       cc          - (struct) Result structure:
%                     .mcc      : Mean connectivity (Z-score or r, depends on zMet).
%                     .mccRaw   : Mean raw correlation (diag removed).
%                     .cc       : Raw correlation matrix (Pearson r).
%                     .zcc    : Processed matrix (Z-score or Masked r).
%                     .mask     : Logical mask of significant pairs (All true for shuffled).
%                     .funcon   : Node degree/strength metric.
%                     .noise    : Noise statistics.
%                     .params   : Input parameters.
%
%   See also DIM_CALC, FR_NETWORK

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'spkmat', @(x) isnumeric(x) || islogical(x));
addParameter(p, 'zMet', 'shuffle', @(x) any(validatestring(x, {'shuffle', 'fisher'})));
addParameter(p, 'nShuffles', 10, @isnumeric);
addParameter(p, 'thresholdVal', 2, @isnumeric);
addParameter(p, 'flgPlot', false, @islogical);

parse(p, spkmat, varargin{:});
spkmat       = double(p.Results.spkmat); % Ensure double for calculations
zMet         = validatestring(p.Results.zMet, {'shuffle', 'fisher'});
nShuffles    = p.Results.nShuffles;
thresholdVal = p.Results.thresholdVal;
flgPlot      = p.Results.flgPlot;


%% ========================================================================
%  PRE-PROCESSING
%  ========================================================================

% Initialize Output
cc.mcc    = NaN;
cc.mccRaw = NaN;
cc.cc     = [];
cc.zcc  = [];
cc.mask   = [];
cc.noise  = struct('mean', NaN, 'std', NaN, 'limit', NaN);
cc.params = p.Results;

% Remove silent neurons
nUnitsRaw = size(spkmat, 1);
mY = mean(spkmat, 2);
validIdx = mY > 0;
nUnits = sum(validIdx);

if nUnits < 2
    cc.cc = nan(nUnitsRaw);
    return;
end

% Z-score valid units
% Normalizing allows using cov() to get Pearson correlation directly
spkz = zscore(spkmat(validIdx, :), 0, 2);


%% ========================================================================
%  RAW CORRELATIONS
%  ========================================================================

% cov expects observations in rows, variables in cols.
% SPKZ is Neurons x Time. We want Neuron-Neuron corr.
% So we treat Time as observations (rows) and Neurons as variables (cols).
r = cov(spkz');

% Remove diagonal (auto-correlations)
r = r - diag(diag(r));

% Raw mean correlation (average off-diagonal elements)
rOffDiag = r(eye(nUnits) == 0);
cc.mccRaw = mean(rOffDiag);

% Raw functinoal connectivity
cc.funconRaw = sum(abs(r), 2) ./ (nUnits - 1);

%% ========================================================================
%  NOISE ESTIMATION
%  ========================================================================

sumR = zeros(nUnits);
sumR2 = zeros(nUnits);
allNoiseVals = []; % For global stats / plotting

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

    % Store off-diagonal for noise distribution (for plotting)
    allNoiseVals = [allNoiseVals; rShuff(eye(nUnits) == 0)]; 
end

% Calculate Per-Pair Noise Stats
noiseAvg = sumR / nShuffles;
noiseVar = (sumR2 - (sumR.^2)/nShuffles) / (nShuffles - 1);
noiseVar = max(noiseVar, 0); % SAFEGUARD: Prevent negative variance due to precision error
noiseStd = sqrt(noiseVar);

% Regularization
% Clamp standard deviations to the "physical floor" of noise variance in
% this network. This prevents Z-score explosion due to under-sampling
% noise. This should not occur with high nShuffles (> 500).
noiseFloor = prctile(noiseStd(eye(nUnits)==0), 5);            
noiseStd(noiseStd < noiseFloor) = noiseFloor;

% Store noise metadata 
cc.noise.mean = noiseAvg;
cc.noise.std  = noiseStd;


%% ========================================================================
%  METRIC CALCULATION
%  ========================================================================

cc.cc = nan(nUnitsRaw);
cc.cc(validIdx, validIdx) = r;

% Initialize Internal Vars
z = zeros(nUnits);
maskValid = false(nUnits);

% Calculate Z-Scores (Standardized Connectivity)
z = (r - noiseAvg) ./ noiseStd;
z = z - diag(diag(z)); 

switch zMet
    case 'shuffle'
        % -----------------------------------------------------------------
        % SHUFFLE Z-SCORE (Standardized Connectivity)
        % -----------------------------------------------------------------
        % Z = (r - mu_null) / sigma_null

        maskValid = true(nUnits) & ~eye(nUnits); % All off-diagonal are valid

        % Mean Connectivity: Average Z-score
        cc.mcc = mean(z(maskValid));

        % Funcon: Mean Absolute Z-score per node
        % normalized by N-1
        funcon = sum(abs(z), 2) ./ (nUnits - 1);

    case 'fisher'
        % -----------------------------------------------------------------
        % FISHER Z + SIGNIFICANCE MASKING
        % -----------------------------------------------------------------
        % Threshold determines significance (e.g., p < 0.05 -> Z > 1.96)
        % Significant r values are Fisher Transformed -> Averaged -> Inverse.

        % Identify significant pairs
        maskValid = abs(z) > thresholdVal;

        % Apply Mask to RAW correlations
        % We keep the MAGNITUDE (r), but only for valid pairs.
        valsSig = r(maskValid);
        z = r .* maskValid; % Zero out non-significant

        % Compute Mean
        if isempty(valsSig)
            cc.mcc = 0;
        else
            % Fisher Transform: atanh(r)
            % Clip to Avoid Inf at r=1/-1 (though unlikely with noise)
            valsSig = max(min(valsSig, 0.999), -0.999);
            zVals = atanh(valsSig);

            % Average Z
            meanZ = mean(zVals);

            % Inverse Fisher: tanh(meanZ)
            cc.mcc = tanh(meanZ);
        end

        % Funcon: Node Degree logic (sum of significant weights / N-1)
        funcon = sum(abs(z), 2) ./ (nUnits - 1);

end


%% ========================================================================
%  EXPAND OUTPUT
%  ========================================================================

cc.zcc = nan(nUnitsRaw);
cc.zcc(validIdx, validIdx) = z;

cc.mask = false(nUnitsRaw);
cc.mask(validIdx, validIdx) = maskValid;

cc.funcon = nan(nUnitsRaw, 1);
cc.funcon(validIdx) = funcon;


%% ========================================================================
%  PLOT
%  ========================================================================

if flgPlot
    figure('Color', 'w', 'Position', [100 100 1200 400]);

    % 1. Raw Matrix (Pearson r)
    subplot(1, 3, 1);
    imagesc(cc.cc, [-0.5 0.5]); % Fixed range for r usually good
    colorbar; axis square;
    title(['Raw Pearson (mccRaw: ' num2str(cc.mccRaw, '%.3f') ')']);
    xlabel('Neuron ID'); ylabel('Neuron ID');

    % 2. Histogram / Metric Dist
    subplot(1, 3, 2);
    hold on;

    if strcmp(zMet, 'fisher')
        % Histogram of Raw Correlations vs Noise
        histogram(allNoiseVals, 50, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', ...
            'DisplayName', 'Noise (Shuff)');
        histogram(rOffDiag, 50, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', ...
            'DisplayName', 'Data (Raw)');
        xline(cc.noise.limit, 'r--', 'Limit');
        xline(-cc.noise.limit, 'r--', 'Limit');
        xlabel('Pearson r');
        title(['Fisher Method (Thr: ' num2str(thresholdVal) '\sigma)']);
    else
        % Histogram of Z-Scores
        zValsFlat = z(eye(nUnits)==0);
        histogram(zValsFlat, 50, 'Normalization', 'pdf', 'DisplayName', 'Data (Z)');
        xline(0, 'k--');
        xlabel('Z-Score');
        title('Shuffle Z-Score Method');
    end
    legend('Location', 'best');
    axis square;

    % 3. Sig/Processed Matrix
    subplot(1, 3, 3);
    imagesc(cc.zcc);
    colorbar; axis square;

    if strcmp(zMet, 'fisher')
        title(['Significant r (mccFisher: ' num2str(cc.mcc, '%.3f') ')']);
    else
        title(['Z-Scores (mccZ: ' num2str(cc.mcc, '%.3f') ')']);
    end
    xlabel('Neuron ID'); ylabel('Neuron ID');

    sgtitle(['Functional Connectivity: ' zMet]);
end

end     % EOF



%% ========================================================================
%  NOTE: METRIC INTERPRETATION
%  ========================================================================
%  This function supports two distinct frameworks for defining "connectivity."
%  The resulting metrics (mcc, funcon) have fundamentally different meanings
%  and units depending on the selected 'zMet'.
%
%  OPTION A: zMet = 'fisher' (MAGNITUDE)
%  -------------------------------------
%  * Unit: Pearson Correlation Coefficient (r), range [-1, 1].
%  * Meaning: Describes the STRENGTH of the coupling.
%  * Interpretation: A value of 0.8 implies that the firing rate of Unit A 
%    linearly predicts Unit B with high accuracy. 
%  * Bias Risk: Highly sensitive to firing rate artifacts. Two bursty units 
%    will have high 'fisher' connectivity simply due to chance overlap of 
%    bursts, even if they are not functionally coupled.
%
%  OPTION B: zMet = 'shuffle' (RELIABILITY)
%  ----------------------------------------
%  * Unit: Z-Score (Standard Deviations), range ~[-5, 20+].
%  * Meaning: Describes the PROBABILITY of the coupling.
%  * Interpretation: A value of 5.0 implies that the observed synchrony is 
%    5 standard deviations above the noise floor. It answers: "How unlikely 
%    is this synchronization to have occurred by chance?"
%  * Advantage: Inherently normalizes for burstiness. Bursty units have 
%    high noise variance, which penalizes their score. This separates 
%    "true reliability" from "rate-driven coincidence."
%  ========================================================================
