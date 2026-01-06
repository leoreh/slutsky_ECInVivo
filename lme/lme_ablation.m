function res = lme_ablation(tbl, frml, lmeCfg, varargin)
% LME_ABLATION Performs feature ablation analysis on a LME/GLME model.
%
%   RES = LME_ABLATION(TBL, FRML, LMECFG, ...)
%   Performs repeated stratified cross-validation to rank feature importance
%   by iteratively removing predictors from the model formula and measuring
%   the drop in performance (Balanced Accuracy and AUC).
%
% INPUTS:
%   tbl             (Required) Table: Data for the model.
%   frml            (Required) Char: Full model formula.
%   lmeCfg          (Optional) Struct: Configuration options.
%                     .dist - Distribution (default 'Binomial').
%   ...             (Optional) Name-Value Pairs:
%                     'nFolds'  - Number of CV folds (default 5).
%                     'nReps'   - Number of CV repetitions (default 5).
%                     'flgPlot' - Plot results (default true).
%
% OUTPUTS:
%   res             Struct containing:
%                     .acc       - Matrix (nTotal x nVars+1) of Bal. Acc.
%                     .auc       - Matrix (nTotal x nVars+1) of AUCs.
%                     .roc       - Cell array (nVars+1) of ROC curves.
%                     .varsAb    - List of ablated variables (1st is 'None').
%                     .impMat    - Matrix of Delta Balanced Accuracy.
%                     .impMatAUC - Matrix of Delta AUC.
%
% DEPENDENCIES:
%   lme_frml2vars, cvpartition, fitglme, predict, perfcurve.

%% ========================================================================
%  INPUT PARSING
%  ========================================================================
p = inputParser;
addRequired(p, 'tbl', @istable);
addRequired(p, 'frml', @ischar);
addOptional(p, 'lmeCfg', struct());
addParameter(p, 'nFolds', 5, @isnumeric);
addParameter(p, 'nReps', 5, @isnumeric);
addParameter(p, 'flgPlot', true, @islogical);
parse(p, tbl, frml, lmeCfg, varargin{:});

nFolds = p.Results.nFolds;
nReps = p.Results.nReps;
flgPlot = p.Results.flgPlot;

% Extract Dist
if isfield(lmeCfg, 'dist') && ~isempty(lmeCfg.dist)
    dist = lmeCfg.dist;
else
    dist = 'Binomial'; % Default context
end

%% ========================================================================
%  PREPARATION
%  ========================================================================

% Extract variables from formula
[varsFxd, varRsp, varsRand] = lme_frml2vars(frml);

% Extract Random Effects (Generic)
if isempty(varsRand)
    randStr = '';
else
    randStr = [' + ' strjoin(varsRand, ' + ')];
end

% Add 'None' to denote the Full Model (no ablation)
% varsIter will be our column mapping: [Full Model, No Var1, No Var2, ...]
varsIter = [{'None'}, varsFxd];

% Common Grid for ROC Interpolation
xGrid = linspace(0, 1, 100);

% Storage
nTotal = nReps * nFolds;
nVars  = length(varsFxd);       % Number of actual predictors
nCols  = length(varsIter);      % Columns = Full + nVars

res.vars  = varsIter;           % Column Labels
res.acc   = nan(nTotal, nCols); % Balanced Accuracy
res.auc   = nan(nTotal, nCols); % AUC
res.roc   = cell(nCols, 1);     % ROC Curves (interpolated)

for iCol = 1:nCols
    res.roc{iCol} = nan(nTotal, length(xGrid));
end

% Initialize Stratified Partition
cvp = cvpartition(tbl.(varRsp), 'KFold', nFolds, 'Stratify', true);

%% ========================================================================
%  VIF (COLLINEARITY)
%  ========================================================================

% Filter: Select only numeric variables from the predictor list
isNum = varfun(@isnumeric, tbl(:, varsFxd), 'OutputFormat', 'uniform');
varsVif = varsFxd(isNum);

if isempty(varsVif)
    res.vif = table();
else
    % Extract data matrix
    X = table2array(tbl(:, varsVif));

    % Correlation matrix
    R = corrcoef(X);

    % VIF = diagonal of inverse correlation matrix
    vifVal = diag(inv(R))';

    % Store
    res.vif = table(varsVif(:), vifVal(:), 'VariableNames', {'Predictor', 'VIF'});
end

%% ========================================================================
%  ABLATION LOOP
%  ========================================================================

fprintf('Starting Ablation Analysis (%d Reps, %d Folds)...\n', nReps, nFolds);

cnt = 0;
for iRep = 1 : nReps

    % Re-partition for new random splits (after first rep)
    if iRep > 1
        cvp = repartition(cvp);
    end

    for iFold = 1 : nFolds
        cnt = cnt + 1;

        tblTrn = tbl(training(cvp, iFold), :);
        tblTst = tbl(test(cvp, iFold), :);
        yTrue  = tblTst.(varRsp);

        % Iterate: Full Model first (iVar=1), then remove each variable
        for iVar = 1 : length(varsIter)

            varRem = varsIter{iVar}; % Variable to remove ('None' = Full)

            % Construct Formula
            % Remove current variable from list
            varsCurr = varsFxd;
            varsCurr(strcmp(varsCurr, varRem)) = [];
            % Rebuild formula: Response ~ Fixed + Random
            currFrml = sprintf('%s ~ %s%s', varRsp, strjoin(varsCurr, ' + '), randStr);

            % Fit Model
            mdl = fitglme(tblTrn, currFrml, 'Distribution', dist);

            % Predict
            yProb = predict(mdl, tblTst);
            
            % Metrics: AUC & ROC
            [x, y, T, aucCurr] = perfcurve(yTrue, yProb, 1);

            % Find Youden's Index (J = Sensitivity + Specificity - 1)
            J = y - x;
            [~, idxJ] = max(J);
            optThresh = T(idxJ);

            % Apply Optimal Threshold
            yPred = double(yProb >= optThresh);

            % Metrics: Balanced Accuracy (at Optimal Threshold)
            tp = sum(yTrue == 1 & yPred == 1);
            tn = sum(yTrue == 0 & yPred == 0);
            sens = tp / max(sum(yTrue == 1), eps);
            spec = tn / max(sum(yTrue == 0), eps);
            accCurr = mean([sens, spec]);

            % Robust Interpolation for ROC
            [uX, idxUni] = unique(x);
            rocCurr = interp1(uX, y(idxUni), xGrid, 'linear', 'extrap');

            % Store Results
            res.acc(cnt, iVar)   = accCurr;
            res.auc(cnt, iVar)   = aucCurr;
            res.roc{iVar}(cnt, :) = rocCurr;

        end
        fprintf('.');
    end
    fprintf('\n');
end

% Calculate Feature Importance (Delta)
% Delta = Full Model (Col 1) - Ablated Model (Col i)
% Matrix will have size (nTotal x nVars), corresponding to varsFxd
res.dacc      = res.acc(:,1) - res.acc(:, 2:end);
res.dauc      = res.auc(:,1) - res.auc(:, 2:end);

% Plotting
if flgPlot
    plot_ablation(res);
end

end     % EOF


%% ========================================================================
%  HELPER: PLOTTING
%  ========================================================================
function plot_ablation(res)

varsFxd = res.vars(2 : end);
nVars  = length(varsFxd);

figure('Color', 'w', 'Name', 'Feature Importance & ROC', 'Position', [50 100 1300 400]);
t = tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

% 1. Feature Importance (Balanced Accuracy)
nexttile;
muImp = mean(res.dacc, 1);
semImp = std(res.dacc, 0, 1) / sqrt(size(res.dacc, 1));
b = bar(muImp);
b.FaceColor = 'k'; b.FaceAlpha = 0.5;
hold on;
errorbar(1:nVars, muImp, semImp, 'k.', 'LineWidth', 1.5);
xticks(1:nVars); xticklabels(varsFxd);
ylabel('\Delta Balanced Accuracy');
title('Importance (BACC)');
grid on;

% 2. Feature Importance (AUC)
nexttile;
muImpAUC = mean(res.dauc, 1);
semImpAUC = std(res.dauc, 0, 1) / sqrt(size(res.dauc, 1));
b = bar(muImpAUC);
b.FaceColor = 'b'; b.FaceAlpha = 0.5;
hold on;
errorbar(1:nVars, muImpAUC, semImpAUC, 'k.', 'LineWidth', 1.5);
xticks(1:nVars); xticklabels(varsFxd);
ylabel('\Delta AUC');
title('Importance (AUC)');
grid on;

% 3. ROC Curves (Comparisons)
nexttile([1, 2]);
hold on;
clrs = lines(nVars);

% Plot Full Model First (Black)
xGrid = linspace(0, 1, 100);
yMeanFull = mean(res.roc{1}, 1);
plot(xGrid, yMeanFull, 'k-', 'LineWidth', 3, 'DisplayName', 'Full Model');

% Plot Reduced Models (Colored)
% Note: res.roc has columns [Full, No Var1, No Var2...]
% So index i+1 corresponds to removing varsAb{i}
for iVar = 1:nVars
    yMeanRed = mean(res.roc{iVar+1}, 1);
    plot(xGrid, yMeanRed, '-', 'Color', clrs(iVar,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('No %s', varsFxd{iVar}));
end

% Formatting
plot([0 1], [0 1], 'k--', 'HandleVisibility','off');
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title(sprintf('Model Comparison (Full AUC = %.2f)', mean(res.auc(:,1))));
legend('Location', 'best');
grid on;

end





%% ========================================================================
%  NOTE: ABLATION ANALYSIS
%  ========================================================================
% To perform an ablation analysis technically, you must iteratively compare
% the predictive performance of a "Full Model" containing all variables
% against several "Reduced Models," each missing exactly one predictor. The
% primary metric to track is **Balanced Accuracy**, which prevents
% performance inflation if your dataset has an unequal number of recovered
% versus non-recovered units. This process requires a cross-validation loop
% to ensure that the accuracy reflects the model's ability to generalize to
% unseen data, rather than just fitting the noise in your current table.
%
% Practically, you should write a loop that modifies your `frml` string for
% each iteration. In each pass, you remove one fixed effect (e.g., `fr` or
% `Group`), fit the model using your `lme_analyse` wrapper, and calculate
% the accuracy on a held-out test set. The resulting "Delta Accuracy" for a
% specific variable is the difference between the full model's accuracy and
% that reduced model's accuracy. This value directly quantifies the unique
% information that the specific predictor—whether unit-level like firing
% rate or network-level like dimensionality—contributes to the prediction
% of unit recovery.
%
%  1. Partition Data:
%     Use `cvpartition` to create training/test sets (e.g., 5-fold).
%
%  2. Full Model Baseline:
%     Fit `fitglme` using the complete formula.
%     Predict outcomes on the test set.
%     Calculate Balanced Accuracy (Bacc_Full).
%
%  3. Iterative Ablation:
%     Loop through each fixed effect in your formula.
%     Create a 'Reduced Formula' by removing one term.
%     Fit and calculate Bacc_Reduced.
%
%  4. Delta Calculation:
%     Delta_Acc = Bacc_Full - Bacc_Reduced.
%
%  BALANCED ACCURACY:
%  Standard accuracy is misleading if 80% of units recovered.
%  Balanced Accuracy = (Sensitivity + Specificity) / 2.
%
%  INTERPRETING THE OUTPUT:
%  - Positive Delta: The variable is essential.
%  - Zero Delta: The variable is redundant.
%  - Negative Delta: The variable may be inducing noise/overfitting.
%  ========================================================================

%% ========================================================================
%  NOTE: EVALUATING BINARY CLASSIFICATION PERFORMANCE
%  ========================================================================
% In physiological research, the evaluation of binary classification
% models—such as predicting whether a single unit will recover from a
% homeostatic perturbation—requires metrics that account for both the
% predictive power of the model and the underlying distribution of the
% data. While "Accuracy" is the most intuitive metric, it is often the
% least informative in biological systems characterized by class imbalance.
%
% Raw accuracy is defined as the total number of correct predictions
% divided by the total number of observations. While mathematically
% straightforward, this metric fails to distinguish between the model's
% ability to identify a biological phenomenon (e.g., unit recovery) and its
% tendency to predict the more frequent outcome. In a dataset where 90% of
% units naturally recover, a null model that predicts "recovery" for every
% unit will achieve 90% accuracy without capturing any physiological
% insight. Consequently, raw accuracy is rarely a sufficient benchmark for
% models involving skewed biological states.
%
% Balanced accuracy addresses the limitations of raw accuracy by
% calculating the arithmetic mean of sensitivity (the true positive rate)
% and specificity (the true negative rate). By normalizing correct
% predictions by the total number of samples within each specific class,
% this metric provides an unbiased assessment of performance. Even if one
% class is significantly larger than the other, balanced accuracy requires
% the model to perform well on both the "recovered" and "unrecovered"
% cohorts to achieve a high score.
% Balanced accuracy is a "point-in-time" metric; it requires you to convert
% a continuous probability into a discrete binary classification (0 or 1)
% by applying a specific threshold, typically 0.5.
%
% The ROC Curve: Threshold-Independent Performance
% The Receiver Operating Characteristic (ROC) curve is a graphical
% representation of a classifier's performance across all possible decision
% thresholds. In a logistic regression model, the output is a continuous
% probability between zero and one; a "threshold" is the point at which we
% decide a probability constitutes a positive prediction (e.g., 0.5). The
% ROC curve plots the True Positive Rate against the False Positive Rate as
% this threshold varies. This visualization reveals the trade-offs between
% sensitivity and specificity, showing how many "false alarms" (false
% positives) must be tolerated to achieve a specific level of "detection"
% (true positives).
%
% AUC: The Global Scalar of Separability The Area Under the Curve (AUC)
% provides a single scalar value that summarizes the entire ROC curve. the
% AUC is calculated by integrating the ROC curve. Conceptually, the AUC is
% equivalent to the Wilcoxon-Mann-Whitney U statistic; it represents the
% probability that the model will assign a higher recovery probability to a
% randomly selected "recovered" unit than to a randomly selected
% "unrecovered" unit. An AUC of 0.5 indicates performance no better than
% chance, while an AUC of 1.0 represents a perfect classifier. Unlike
% accuracy or balanced accuracy, which depend on a fixed decision
% threshold, the AUC evaluates the model's fundamental ability to separate
% the two physiological states regardless of where the cutoff is set.
%
% Technical Selection of Metrics
% The choice between these metrics depends on the goal of the analysis.
% Balanced accuracy is the preferred metric for reporting how well a
% specific, implemented classifier performs in a real-world predictive
% task. In contrast, the AUC is the gold standard for comparing the
% intrinsic quality of different models or features—such as comparing the
% predictive weight of firing rates versus network dimensionality—because
% it is immune to the biases introduced by class distribution and arbitrary
% threshold selection.
%% ========================================================================


%% ========================================================================
%  NOTE: COLLINEARITY
%  ========================================================================
% Collinearity, or multicollinearity, occurs when two or more independent
% predictors in your model are highly correlated with one another. In
% physiological datasets, this is a common challenge because many
% metrics—such as mean firing rate, burst fraction, and network
% dimensionality—often capture overlapping aspects of the same underlying
% neural state. Technically, collinearity means that one predictor can be
% linearly predicted from the others with a high degree of accuracy, which
% violates the assumption that each predictor provides unique information
% to the model.
%
% While collinearity does not necessarily reduce the overall predictive
% power of the model, it severely compromises the interpretability of
% individual coefficients. In a mixed-effects framework, high collinearity
% inflates the standard errors of the fixed effects, making it difficult to
% achieve statistical significance even for variables with strong
% biological effects. Furthermore, the estimated coefficients become highly
% unstable; small changes in your dataset can cause the sign of a
% coefficient to flip from positive to negative, leading to erroneous
% physiological conclusions.
%
% To assess whether you need to intervene, you should calculate the
% Variance Inflation Factor (VIF) for each continuous predictor. The VIF
% quantifies how much the variance of an estimated coefficient is increased
% due to collinearity. A general rule of thumb in the literature is that a
% VIF exceeding 5 to 10 indicates problematic redundancy. Additionally,
% generating a pairwise correlation matrix of your predictors before
% running `fitglme` provides a direct visual assessment of which variables
% are "tracking" each other too closely.
%
%  Impact:
%     - Unstable coefficient estimates (Beta).
%     - High sensitivity to outliers.
%     - Masked significance of predictors.
%
% See:
% https://www.mathworks.com/help/econ/time-series-regression-ii-collinearity-and-estimator-variance.html
%
%  If VIF analysis reveals severe redundancy between variables (e.g., fr
%  and dimensionality), consider the following technical steps:
%
%  1. Variable Removal:
%     Drop the less biologically relevant predictor.
%
%  2. Dimensionality Reduction (PCA):
%     Combine correlated metrics into one.
%     Use the first Principal Component (PC1).
%
%  3. Feature Centering:
%     Standardizing continuous variables reduces collinearity.
%     Particularly helpful for interaction terms.
%  ========================================================================

%% ========================================================================
%  NOTE: STEPWISE SELECTION
%  ========================================================================
%  MATLAB's built-in `stepwiseglm` does NOT support random effects.
%
% While MATLAB’s `stepwisefit` and `stepwiseglm` are standard tools for
% automated model selection, they are technically incompatible with your
% current analytical framework. These functions are designed for
% fixed-effects models and cannot accommodate the random-effects structure
% (`1|Name`) required for hierarchical physiological data. Using them would
% require you to ignore the non-independence of units within a culture,
% which significantly inflates Type I error rates and invalidates your
% statistical inferences.
%
% The `stepwiseglm` function operates by iteratively adding or removing
% predictors based on a specific criterion, such as -values, Akaike
% Information Criterion (AIC), or Bayesian Information Criterion (BIC).
% However, this "greedy" search algorithm is strictly limited to the
% `GeneralizedLinearModel` class. Because your recovery analysis requires a
% `GeneralizedLinearMixedModel` to account for mouse-specific or
% culture-specific variance, you cannot use these built-in automated
% routines. Furthermore, stepwise selection is frequently criticized in
% high-impact literature for -value bias, as it performs multiple hidden
% comparisons that are not reflected in the final model's output.
%
% There is a fundamental theoretical difference between the stepwise
% approach and the ablation analysis performed by Atsmon or Parks. Stepwise
% selection is a **search strategy** intended to find the most parsimonious
% model (the "best" subset of predictors). In contrast, ablation analysis
% is a **diagnostic tool** used to quantify the "predictive weight" of each
% variable within a pre-defined physiological hypothesis. While stepwise
% might tell you which variables to keep, ablation tells you how much
% information is lost if a variable is removed, which is the metric used to
% generate "Importance" or "Delta Accuracy" plots.
%
% To achieve an effect similar to stepwise selection while maintaining your
% random-effects structure, you must perform manual model comparison using
% the `compare` function in MATLAB. This allows you to test whether a model
% with a specific predictor is statistically superior to one without it
% using a Likelihood Ratio Test (LRT). However, for the specific "Delta
% Accuracy" goal you described, the manual loop strategy remains the
% superior technical path.
%
%  To technically determine if a predictor (e.g., 'Group') should be
%  included in your GLME, compare nested models:
%
%  - Model 1: uRcv ~ fr + bFrac + (1|Name)
%  - Model 2: uRcv ~ fr + bFrac + Group + (1|Name)
%
%  [~, p, stats] = compare(Mdl1, Mdl2);
%
%  If p < 0.05, the added predictor significantly improves the model fit
%  beyond what is expected by chance.
%
%  Ablation is a Machine Learning approach, not a Likelihood approach.
%  1. Purpose: To rank feature importance, not just "significance."
%  2. Performance Metric: Usually "Delta Balanced Accuracy" or "Delta AUC."
%  3. Implementation Requirement: Must be done via a custom loop to
%     preserve the `fitglme` random-effect structure in every iteration.
%  ========================================================================


%% ========================================================================
%  NOTE: FIXED VS. OPTIMAL THRESHOLDS (OPERATING POINTS)
%  =======================================================================
% When evaluating classification performance (e.g., uRcv), the choice of
% decision threshold - the "Operating Point" - directly dictates the
% calculated Balanced Accuracy (BACC).
%
% 1. The Standard Operating Point (p = 0.5):
%    By default, MATLAB's `predict` and `rocmetrics` utilize a threshold of
%    0.5. While standard for logistic regression, this fixed
%    cutoff often fails in physiological data characterized by class
%    imbalance or high variance in predictors like Firing Rate. If a
%    variable shows high AUC but negative Delta BACC, it typically indicates
%    that the 0.5 threshold is biologically suboptimal for that feature.
%
% 2. The Optimal Operating Point (Youden’s Index):
%    The optimal threshold is the point on the ROC curve that maximizes
%    separability (J = Sensitivity + Specificity - 1). This
%    point minimizes the "distance" to the top-left corner of the ROC space.
%    Using an adaptive threshold ensures that BACC reflects the true
%    discriminatory power of a predictor rather than its alignment with
%    an arbitrary 0.5 cutoff.
%
% 3. Technical Implementation:
%    To avoid "threshold-bias" in importance plots, calculate the best
%    threshold within each cross-validation fold using:
%
% 4. Interpretation:
%    - Fixed 0.5: Measures "Real-world" utility of a generic classifier.
%    - Optimal: Measures the "Theoretical Maximum" predictive weight of
%      biological features.
%  ========================================================================