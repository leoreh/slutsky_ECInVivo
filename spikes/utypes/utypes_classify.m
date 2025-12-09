function unitType = utypes_classify(varargin)

% UTYPES_CLASSIFY classifies neuronal units into Regular Spiking (RS) and
% Fast Spiking (FS) types using a two-step approach:
% 1. Initial classification using trough-to-peak (tp) threshold.
% 2. Optional refinement using GMM clustering starting from the threshold-based
%    classification.
%
% METHODOLOGY:
% Method 1 (altClassify = 1): Uses simple threshold-based classification on tp.
% Method 2 (altClassify = 2): Uses GMM clustering initialized with threshold-based
%                            classification results.
%
% INPUT (Optional Key-Value Pairs):
%   basepaths    (cell array) Full paths to recording folders. If empty,
%                uses current directory.
%   altClassify  (numeric) Classification method to use:
%                1 - Threshold-based classification using tp
%                2 - GMM-based classification using multiple metrics, initialized
%                    with threshold-based results
%                {1}
%   flgSave      (logical) Flag to save classification results
%                {false}
%   v            (struct) Data structure containing metrics. If empty,
%                data will be loaded from basepaths. {[]}
%   fetTbl       (table) Table containing features, if already loaded in memory
%
% OUTPUT:
%   unitType     (vector) Vector of unit classifications where:
%                1 = Regular Spiking (RS) - Typically Pyramidal
%                2 = Fast Spiking (FS) - Typically Interneuron
%
% DEPENDENCIES:
%   basepaths2vars, catfields
%
% HISTORY:
%   LH - Aug 2024
%   Based in part on Oghazian et al., Front. Biomedical Tech., 2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments using inputParser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments using inputParser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addOptional(p, 'fetTbl', table(), @istable);
addOptional(p, 'altClassify', 1, @(x) ismember(x, [1, 2, 3]));
addOptional(p, 'flgSave', false, @islogical);
addOptional(p, 'basepaths', {}, @(x) iscell(x));
addOptional(p, 'fetSelect', {}, @iscell);

parse(p, varargin{:});
fetTbl = p.Results.fetTbl;
altClassify = p.Results.altClassify;
flgSave = p.Results.flgSave;
basepaths = p.Results.basepaths;
fetSelect = p.Results.fetSelect;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get number of units
nUnits = height(fetTbl);
initType = fetTbl.initType;

% Identify good units (those that were classified)
idxGood = ~fetTbl.bad;

% Convert categorical initType to numeric codes 1=RS, 2=FS
unitType = zeros(size(initType));
unitType(initType == 'RS') = 1;
unitType(initType == 'FS') = 2;

% Define features if not provided
if isempty(fetSelect) && altClassify > 1
    switch altClassify
        case 2
            fetSelect = {'tp', 'lidor', 'mfr'};
        case 3
            fetSelect = {'asym', 'hpk', 'tp'};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% refine classification if requested
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Refine classification using GMM with threshold results as starting point
if ~isempty(fetSelect)

    % Prepare feature matrix for GMM
    fet = fetTbl{idxGood, fetSelect};

    % Get refined classification using GMM (on good units only)
    unitType(idxGood) = gm2units(fet, unitType(idxGood));
end

% Create clean matrix for units struct
clean = false(2, nUnits);
clean(1, idxGood) = unitType(idxGood) == 1;  % RS
clean(2, idxGood) = unitType(idxGood) == 2;  % FS

% Update unitType to span all units
unitType = zeros(nUnits, 1);
unitType(clean(1, :)) = 1;
unitType(clean(2, :)) = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save updated unit struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flgSave
    startIdx = 1;
    for iPath = 1 : length(basepaths)
        basepath = basepaths{iPath};
        [~, basename] = fileparts(basepath);
        uFile = fullfile(basepath, [basename, '.units.mat']);

        % Get number of units for this recording
        units = v(iPath).units;
        nUnits = size(units.clean, 2);
        uIdx = startIdx : startIdx + nUnits - 1;
        startIdx = uIdx(end) + 1;

        % Add new classification field
        fieldName = sprintf('clean%d', altClassify);
        units.(fieldName) = false(2, nUnits);
        units.(fieldName)(:, :) = clean(:, uIdx);

        % Place altClassify3 as default
        if altClassify == 3
            units.clean = clean(:, uIdx);
        end

        % Save updated units structure
        save(uFile, 'units');
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function typeOut = gm2units(fet, typeIn)
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
gm = fitgmdist(fet, 2, 'Start', typeIn,...
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
rsProbability = 0.8;
priors(rsIdx) = rsProbability;
priors(fsIdx) = 1 - rsProbability;
gm = gmdistribution(gm.mu, gm.Sigma, priors);

% Final Classification.
% Cluster assignment is now: Maximize (Likelihood * FixedPrior)
clusterLabels = cluster(gm, fet);

% Map internal GMM components back to output format (1=RS, 2=FS)
typeOut = ones(size(typeIn));          % Initialize all as RS
typeOut(clusterLabels ~= rsIdx) = 2;   % Set non-RS cluster to FS

end

% EOF



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
%  ========================================================================