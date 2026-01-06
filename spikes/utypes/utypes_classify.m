function tblUnit = utypes_classify(varargin)

% UTYPES_CLASSIFY classifies units into Regular Spiking (RS) and Fast
% Spiking (FS) types using a Gaussian Mixture Model (GMM). See notes at
% end.
%
% INPUT (Optional Key-Value Pairs):
%   basepaths    (cell array) Full paths to recording folders.
%   fetSelect    (numeric or cell) Features to use for GMM classification.
%   flgPlot      (logical) Flag to plot classification results.
%   rsPrior      (numeric) Prior probability of RS unit. {0.9}
%   regVal       (numeric) Regularization value. {0.01}
%
% OUTPUT:
%   tblUnit      (table) Table containing features, UnitType, etc.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addOptional(p, 'basepaths', {}, @(x) iscell(x));
addOptional(p, 'fetSelect', 1, @(x) iscell(x) || isnumeric(x));
addOptional(p, 'flgPlot', false, @islogical);
addOptional(p, 'rsPrior', 0.9, @isnumeric);
addOptional(p, 'regVal', 0.01, @isnumeric);

parse(p, varargin{:});
basepaths = p.Results.basepaths;
fetSelect = p.Results.fetSelect;
flgPlot = p.Results.flgPlot;
rsPrior = p.Results.rsPrior;
regVal = p.Results.regVal;

%% ========================================================================
%  CREATE TABLE
%  ========================================================================

if isempty(basepaths)
    basepaths = [mcu_basepaths('wt'), mcu_basepaths('mcu')];
end

vars = {'swv_metrics', 'st_metrics', 'fr', 'units'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

% Metadata
tagFiles = struct();
tagFiles.Mouse = get_mname(basepaths);
[~, fNames] = fileparts(basepaths);
tagFiles.File = fNames;

% Variable Map
varMap = struct();
varMap.FR = 'fr.mfr';
varMap.BRoy = 'st.royer';
varMap.BLidor = 'st.lidor';
varMap.BMiz = 'st.mizuseki';
varMap.TP = 'swv.tp';
varMap.Inverted = 'swv.inverted';
varMap.TPAmp = 'swv.tpAmp';
varMap.TPRatio = 'swv.tpRatio';
varMap.TPSlope = 'swv.tpSlope';
varMap.SpkW = 'swv.spkw';
varMap.Asym = 'swv.asym';
varMap.Hpk = 'swv.hpk';
varMap.TailSlope = 'swv.tailSlope';
varMap.TailAmp = 'swv.tailAmp';

tblUnit = v2tbl('v', v, 'varMap', varMap, 'tagFiles', tagFiles);
tblUnit = tbl_transform(tblUnit, 'logBase', 10);


%% ========================================================================
%  PRELIMINARY CLASSIFICATION
%  ========================================================================

% Bad units
tblUnit.bad = tblUnit.Inverted | isnan(tblUnit.TP);
idxGood = ~tblUnit.bad;

if sum(idxGood) == 0
    warning('No good units found.');
    return;
end

% Order: 1=RS, 2=FS, 3=Other
InitType = ones(height(tblUnit), 1) * 3;

% RS: TP >= 0.4
InitType(idxGood & (tblUnit.TP >= 0.4)) = 1;
% FS: TP < 0.4
InitType(idxGood & (tblUnit.TP < 0.4)) = 2;

tblUnit.InitType = categorical(InitType, [1, 2, 3], {'RS', 'FS', 'Other'});


%% ========================================================================
%  REFINE CLASSIFICATION (GMM)
%  ========================================================================

UnitType = double(tblUnit.InitType); % 1=RS, 2=FS, 3=Other

% Handle fetSelect presets
if isnumeric(fetSelect)
    switch fetSelect
        case 1
            fetSelect = {'Asym', 'Hpk', 'TP', 'BLidor', 'FR'};
        case 2
            fetSelect = {'Asym', 'Hpk', 'TP'};
        otherwise
            fetSelect = {};
    end
elseif iscell(fetSelect)
    fetSelect = fetSelect;
else
    fetSelect = {};
end

if ~isempty(fetSelect)
    % Prepare feature matrix
    fet = tblUnit{idxGood, fetSelect};

    % Only classified RS(1) or FS(2) units are refined
    typeIn = UnitType(idxGood);

    % GMM Refinement
    refinedType = gm2units(fet, typeIn, rsPrior, regVal);

    UnitType(idxGood) = refinedType;
end

tblUnit.UnitType = categorical(UnitType, [1, 2, 3], {'RS', 'FS', 'Other'});

%% ========================================================================
%  AUC CALCULATION
%  ========================================================================

idxNumeric = varfun(@isnumeric, tblUnit, 'OutputFormat', 'uniform');
nVar = length(idxNumeric);
varAuc = zeros(nVar, 1);
UnitType = double(tblUnit.UnitType(idxGood));

for iVar = 1:nVar
    if idxNumeric(iVar)
        [~,~,~,varAuc(iVar)] = perfcurve(UnitType, tblUnit{idxGood, iVar}, 1);
    end
end

%% ========================================================================
%  PLOT
%  ========================================================================

if flgPlot
    utypes_gui('tblUnit', tblUnit, 'basepaths', basepaths);
end

end

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function UnitType = gm2units(fet, typeIn, rsPrior, regVal)
% GM2UNITS Refines unit classification using Gaussian Mixture Model

% Learn the cluster shapes
gm = fitgmdist(fet, 2, 'Start', typeIn, 'RegularizationValue', regVal,...
    'CovarianceType', 'diagonal');

% Assign clusters
clusterLabels = cluster(gm, fet);
rsCounts = zeros(1, 2);
for iUnit = 1:2
    rsCounts(iUnit) = sum(clusterLabels == iUnit & typeIn == 1);
end
[~, rsIdx] = max(rsCounts);
fsIdx = 3 - rsIdx;

% Apply Bayesian Priors
priors = zeros(1, 2);
priors(rsIdx) = rsPrior;
priors(fsIdx) = 1 - rsPrior;
gm = gmdistribution(gm.mu, gm.Sigma, priors);

% Final Classification
clusterLabels = cluster(gm, fet);

% Map back
UnitType = ones(size(typeIn));          % RS
UnitType(clusterLabels ~= rsIdx) = 2;   % FS

end


%% ========================================================================
%  NOTE: FEATURE VALIDATION VIA AUC
%  ========================================================================
%  We use the Area Under the ROC Curve (AUC) to validate whether the
%  waveform-based classification predicts other physiological properties.
%
%  1. Ground Truth:
%     The classification labels (RS vs. FS) are defined
%     solely by the Trough-to-Peak (TP) waveform duration.
%
%  2. Interpretation of AUC:
%     The AUC quantifies how well *other* metrics (e.g., Firing Rate)
%     align with these waveform labels.
%     - AUC ≈ 1.0 or 0.0: Strong Correlation. The feature is an effective
%       predictor. (Note: Values < 0.5 simply indicate an inverse
%       relationship, e.g., Higher Rate = Lower Probability of RS).
%     - AUC ≈ 0.5: Random Chance. The feature provides no discriminative
%       information for this classification.
%  ========================================================================


%% ========================================================================
%  NOTE: REDUNDANCY ELIMINATION
%  ========================================================================
%  Feature selection is critical for Gaussian Mixture Model (GMM) performance.
%  GMMs rely on the Mahalanobis distance, which requires inverting the
%  covariance matrix. Highly correlated features (multicollinearity) cause
%  this matrix to become singular or unstable, leading to poor clustering.
%
%  1. Redundancy Criteria:
%     We assess redundancy using Spearman's Rank Correlation. Features with
%     |rho| > 0.85 are considered redundant.
%     - Visual Confirmation: Redundant pairs manifest as tight diagonal
%       lines in pairwise scatter plots.
%     - Independence: Effective feature pairs appear as "clouds" or
%       orthogonal distributions (rho < 0.3), providing unique information
%       to the classifier.
%
%  2. Selection Strategy (The "Survivor" Rule):
%     When a redundant pair is identified, we retain the feature that is:
%     - More Robust: Less sensitive to noise (e.g., Prefer 'FullWidth' over
%       'TailSlope').
%     - More Gaussian: Exhibits a distribution closer to Normality (or
%       Log-Normality), fitting the GMM assumptions better.
%     - Biologically Interpretable: Prioritize metrics with direct
%       physiologial correlates (e.g., Repolarization Time).
%  ========================================================================


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
