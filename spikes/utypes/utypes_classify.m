function uTbl = utypes_classify(varargin)

% UTYPES_CLASSIFY classifies units into Regular Spiking (RS) and Fast
% Spiking (FS) types using a Gaussian Mixture Model (GMM). It returns a
% table containing the classification results and the features used for
% classification.
%
% INPUT (Optional Key-Value Pairs):
%   basepaths    (cell array) Full paths to recording folders. If empty,
%                uses current directory.
%   fetSelect    (numeric or cell) Features to use for GMM classification.
%                If numeric:
%                   1 - {'asym', 'hpk', 'tp', 'lidor', 'mfr'}
%                   2 - {'asym', 'hpk', 'tp'}
%                If cell array: List of feature names to use.
%                If empty: No GMM refinement is performed (returns initial).
%                {1}
%   flgSave      (logical) Flag to save classification results
%                {false}
%   flgPlot      (logical) Flag to plot classification results
%                {false}
%   fetTbl       (table) Table containing features.
%   rsPrior      (numeric) Prior probability of RS unit. {0.9}
%   regVal       (numeric) Regularization value for GMM covariance. {0.01}
%
% OUTPUT:
%   uTbl         (table) Table containing the features used for classification
%                and the final 'UnitType' column.
%
% HISTORY:
%   LH - Aug 2024
%   Based in part on Oghazian et al., Front. Biomedical Tech., 2021

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addOptional(p, 'fetTbl', table(), @istable);
addOptional(p, 'fetSelect', 1, @(x) iscell(x) || isnumeric(x));
addOptional(p, 'flgSave', false, @islogical);
addOptional(p, 'flgPlot', false, @islogical);
addOptional(p, 'basepaths', {}, @(x) iscell(x));
addOptional(p, 'rsPrior', 0.9, @isnumeric);
addOptional(p, 'regVal', 0.01, @isnumeric);

parse(p, varargin{:});
fetTbl = p.Results.fetTbl;
fetSelectInput = p.Results.fetSelect;
flgSave = p.Results.flgSave;
flgPlot = p.Results.flgPlot;
basepaths = p.Results.basepaths;
rsPrior = p.Results.rsPrior;
regVal = p.Results.regVal;

%% ========================================================================
%  PREPARATIONS
%  ========================================================================

% get number of units
initType = fetTbl.initType;

% Identify good units (those that were classified)
idxGood = ~fetTbl.bad;

% Convert categorical initType to numeric codes 1=RS, 2=FS
UnitType = zeros(size(initType));
UnitType(initType == 'RS') = 1;
UnitType(initType == 'FS') = 2;

% Handle fetSelect presets
if isnumeric(fetSelectInput)
    switch fetSelectInput
        case 1
            fetSelect = {'asym', 'hpk', 'tp', 'lidor', 'mfr'};
        case 2
            fetSelect = {'asym', 'hpk', 'tp'};
        otherwise
            fetSelect = {};
    end
elseif iscell(fetSelectInput)
    fetSelect = fetSelectInput;
else
    fetSelect = {};
end

%% ========================================================================
%  REFINE CLASSIFICATION
%  ========================================================================

% Refine classification using GMM with threshold results as starting point
if ~isempty(fetSelect)

    % Prepare feature matrix for GMM
    fet = fetTbl{idxGood, fetSelect};

    % Get refined classification using GMM (on good units only)
    UnitType(idxGood) = gm2units(fet, UnitType(idxGood), rsPrior, regVal);
end

% Assign to table
fetTbl.UnitType = categorical(UnitType, [0, 1, 2], {'Other', 'RS', 'FS'});

%% ========================================================================
%  PLOT & SAVE
%  ========================================================================

if flgPlot
    mcuCfg = mcu_cfg;
    cfg.xVar = 'tp';
    cfg.yVar = 'lidor';
    cfg.szVar = 'mfr';
    cfg.grpVar = 'UnitType';
    cfg.clr = mcuCfg.clr.unit([3 : -1 : 1], :);
    cfg.alpha = 0.4;

    plot_utypes('uTbl', fetTbl, 'cfg', cfg, 'flgSave', flgSave, 'basepaths', basepaths);
end

% Create output table with only used features and UnitType
varsToKeep = unique([fetSelect(:)', {'UnitType'}], 'stable');
uTbl = fetTbl(:, varsToKeep);


end     % EOF

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function typeOut = gm2units(fet, typeIn, rsPrior, regVal)
% GM2UNITS Refines unit classification using Gaussian Mixture Model
%
% This helper function performs GMM clustering to refine the initial unit
% classification. It uses the initial classification as a starting point
% for the GMM to improve convergence and interpretability.
%
% INPUT:
%   fet      (matrix) Feature matrix for GMM [nUnits x nFeatures]
%   typeIn   (vector) Initial classification [nUnits x 1]
%                     1 = Regular Spiking (RS)
%                     2 = Fast Spiking (FS)
%
% OUTPUT:
%   typeOut  (vector) Refined classification [nUnits x 1]

% Learn the cluster shapes (Means and Covariances) from the data, using
% 'typeIn' as starting point
gm = fitgmdist(fet, 2, 'Start', typeIn, 'RegularizationValue', regVal,...
    'CovarianceType', 'diagonal');

% Assign clusters to RS and FS based on the initial classification
clusterLabels = cluster(gm, fet);
rsCounts = zeros(1, 2);
for iUnit = 1:2
    rsCounts(iUnit) = sum(clusterLabels == iUnit & typeIn == 1);
end
[~, rsIdx] = max(rsCounts);
fsIdx = 3 - rsIdx;

% Apply Bayesian Priors and Reconstruct the GMM. This shifts the decision
% boundary without distorting the cluster shapes.
priors = zeros(1, 2);
priors(rsIdx) = rsPrior;
priors(fsIdx) = 1 - rsPrior;
gm = gmdistribution(gm.mu, gm.Sigma, priors);

% Final Classification.
% Cluster assignment is now: Maximize (Likelihood * FixedPrior)
clusterLabels = cluster(gm, fet);

% Map internal GMM components back to output format (1=RS, 2=FS)
typeOut = ones(size(typeIn));          % Initialize all as RS
typeOut(clusterLabels ~= rsIdx) = 2;   % Set non-RS cluster to FS

end



%% ========================================================================
%  NOTE: BAYESIAN PRIORS & POPULATION RATIO ENFORCEMENT
%  ========================================================================
%  The classification utilizes a Gaussian Mixture Model (GMM) with
%  informed Bayesian Priors to address potential sampling biases in
%  neuronal recordings.
%
%  1. The Bayesian Framework:
%     Standard GMM classification assigns units to clusters by maximizing
%     the Posterior Probability: P(Type|Waveform). According to Bayes'
%     Theorem, this is proportional to:
%         P(Waveform|Type)  * P(Type)
%         [Likelihood]         [Prior]
%
%  2. Likelihood (The Shape):
%     P(Waveform|Type) represents how well a unit's metrics fit the
%     geometric cluster (mean and variance) of that cell type. This constitutes
%     the "objective" fit of the data to the cluster distribution.
%
%  3. The Prior (The Population Count):
%     P(Type) represents the expected abundance of that cell type in the
%     population. In standard unsupervised clustering, this is estimated
%     directly from the recording data. However, electrophysiological
%     recordings often suffer from sampling bias (e.g., under-sampling
%     small interneurons), leading to incorrect data-driven priors.
%
%  4. Informed Priors vs. Forcing:
%     By fixing the Prior (e.g., setting FS = 0.2), we do not "force"
%     dissimilar units into the FS cluster. We merely adjust the
%     decision boundary. A unit with a distinct FS waveform will
%     still be classified as such because its Likelihood term is high.
%     However, for "ambiguous" units lying between clusters, the fixed Prior
%     shifts the classification probability toward the more common cell type
%     (RS). This reduces False Positives for the rarer class and
%     aligns the results with established histological ratios (e.g., ~80/20
%     RS/FS split).
%
%  5. In the mouse CA1 hippocampus, the actual biological ratio is
%     approximately 85% to 90% Pyramidal. Still, a prior of 80% RS units
%     is conservative, given electrodes often have a sampling bias toward
%     Interneurons because they have higher firing rates (making them
%     easier to detect/sort) compared to silent Pyramidal cells.
%  ========================================================================

%% ========================================================================
%  NOTE: FEATURE TRANSFORMATION & POOLING STRATEGY
%  ========================================================================
%  1. Firing rates and burst indices typically follow
%     Log-Normal distributions. We apply a Log10 transform to normalize
%     these inputs, satisfying the Gaussian assumption required by the GMM.
%
%  2. Units are pooled across all animals prior to
%     clustering. This ensures sufficient sample size for stable covariance
%     estimation and relies on the fact that spike duration is a conserved
%     biophysical constant across subjects.
%
%  3. We explicitly avoid Z-scoring within animals.
%     In datasets with sampling bias (e.g., an animal with only RS
%     cells), local Z-scoring incorrectly shifts the data mean to 0,
%     causing RS cells to be misclassified as FS. Pooling
%     preserves the absolute decision boundaries.
%  ========================================================================

%% ========================================================================
%  NOTE: REGULARIZATION (NUMERICAL STABILITY)
%  ========================================================================
%  An 'RegularizationValue' can be added to the GMM to ensure
%  mathematical stability during the Expectation-Maximization (EM) algorithm.
%
%  1. The Problem (Singularity):
%     GMMs require inverting the covariance matrix to calculate the
%     Mahalanobis distance. If features are highly correlated (redundant)
%     or if a cluster contains very few points (collapse), the covariance
%     matrix becomes "singular" (non-invertible), causing the code to crash.
%
%  2. The Solution:
%     Regularization adds a small constant (epsilon) to the diagonal of the
%     covariance matrix: Sigma_new = Sigma + epsilon * IdentityMatrix.
%     This artificially inflates the variance slightly, ensuring the matrix
%     is always Positive Definite and invertible.
%
%  3. The Value:
%     This specific value is likely empirical. Standard values typically
%     range from 0.0001 to 0.1 depending on the scale of the data.
%     If the value is too high, it "blurs" the clusters (artificially
%     widening them). If too low, the code may error on outliers.