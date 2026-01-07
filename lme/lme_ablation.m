function res = lme_ablation(tbl, frml, varargin)
% LME_ABLATION Performs feature ablation analysis on a LME/GLME model.
%
%   RES = LME_ABLATION(TBL, FRML, ...)
%   Performs repeated cross-validation to rank feature importance by
%   iteratively removing predictors from the model formula and measuring the
%   drop in performance.
%
%   MODES:
%       - Classification (Binomial): Uses Balanced Accuracy & AUC.
%       - Regression (Continuous): Uses RMSE & R-Squared.
%
% INPUTS:
%   tbl             (Required) Table: Raw data table.
%   frml            (Required) Char: Full model formula.
%   ...             (Optional) Name-Value Pairs:
%                     'dist'    - Distribution (default auto-select).
%                     'nFolds'  - Number of CV folds (default 5).
%                     'nReps'   - Number of CV repetitions (default 5).
%                     'flgPlot' - Plot results (default true).
%
% OUTPUTS:
%   res             Struct containing:
%                     .vars     - List of models (Full, No Var1, ...).
%                     .perfPrim - Matrix of primary metric (BACC or RMSE).
%                     .perfSec  - Matrix of secondary metric (AUC or R2).
%                     .impPrim  - Delta Primary Metric.
%                     .vif      - Variance Inflation Factors.
%
% DEPENDENCIES:
%   lme_fit, lme_frml2vars, cvpartition.

%% ========================================================================
%  INPUT PARSING
%  ========================================================================
p = inputParser;
addRequired(p, 'tbl', @istable);
addRequired(p, 'frml', @ischar);
addParameter(p, 'dist', '', @ischar);
addParameter(p, 'nFolds', 5, @isnumeric);
addParameter(p, 'nReps', 5, @isnumeric);
addParameter(p, 'flgPlot', true, @islogical);
parse(p, tbl, frml, varargin{:});

distInput = p.Results.dist;
nFolds = p.Results.nFolds;
nReps = p.Results.nReps;
flgPlot = p.Results.flgPlot;

% Run LME_ANALYSE to prepare data and select distribution
fprintf('[LME_ABLATION] Running initial analysis to prepare data...\n');
[~, ~, lmeInfo, lmeTbl] = lme_analyse(tbl, frml, 'dist', distInput);
dist = lmeInfo.distFinal;

% Determine Mode
isBinomial = strcmpi(dist, 'Binomial');
if isBinomial
    metricPrim = 'Balanced Accuracy';
    metricSec  = 'AUC';
    modeStr = 'Classification';
else
    metricPrim = 'RMSE';
    metricSec  = 'R-Squared';
    modeStr = 'Regression';
end

fprintf('[LME_ABLATION] Mode: %s (Dist: %s)\n', modeStr, dist);


%% ========================================================================
%  PREPARATION
%  ========================================================================

% Extract variables from formula
[varsFxd, varRsp, varsRand] = lme_frml2vars(frml);

% Extract Random Effects String
if isempty(varsRand)
    randStr = '';
else
    randStr = [' + ' strjoin(varsRand, ' + ')];
end

% varsIter: [Full Model, No Var1, No Var2, ...]
varsIter = [{'None'}, varsFxd];

% Storage
nTotal = nReps * nFolds;
nCols  = length(varsIter);

res.vars   = varsIter;
res.perfPrim = nan(nTotal, nCols); % BACC or RMSE
res.perfSec  = nan(nTotal, nCols); % AUC or R2
res.roc    = cell(nCols, 1);     % ROC (Binomial only)
xGrid      = linspace(0, 1, 100); % Common grid for ROC

for iCol = 1:nCols
    res.roc{iCol} = nan(nTotal, length(xGrid));
end

% Initialize CV Partition
if isBinomial
    cvp = cvpartition(lmeTbl.(varRsp), 'KFold', nFolds, 'Stratify', true);
else
    cvp = cvpartition(height(lmeTbl), 'KFold', nFolds);
end


%% ========================================================================
%  VIF (COLLINEARITY)
%  ========================================================================

% Calculate VIF on transformed predictors
isNum = varfun(@isnumeric, lmeTbl(:, varsFxd), 'OutputFormat', 'uniform');
varsVif = varsFxd(isNum);

if isempty(varsVif)
    res.vif = table();
else
    X = table2array(lmeTbl(:, varsVif));
    R = corrcoef(X, 'Rows', 'pairwise');

    % Handle singular matrix
    if rcond(R) < 1e-12
        vifVal = nan(1, length(varsVif));
        warning('Predictors are perfectly collinear. VIF calculation skipped.');
    else
        vifVal = diag(inv(R))';
    end

    res.vif = table(varsVif(:), vifVal(:), 'VariableNames', {'Predictor', 'VIF'});
end


%% ========================================================================
%  ABLATION LOOP
%  ========================================================================

fprintf('Starting Ablation (%d Reps, %d Folds)...\n', nReps, nFolds);

cnt = 0;
for iRep = 1 : nReps

    if iRep > 1
        cvp = repartition(cvp);
    end

    for iFold = 1 : nFolds
        cnt = cnt + 1;

        trainIdx = training(cvp, iFold);
        testIdx  = test(cvp, iFold);

        tblTrn = lmeTbl(trainIdx, :);
        tblTst = lmeTbl(testIdx, :);
        yTrue  = tblTst.(varRsp);

        % Iterate Models
        for iVar = 1 : length(varsIter)

            varRem = varsIter{iVar};

            % Construct Formula
            varsCurr = varsFxd;
            if ~strcmp(varRem, 'None')
                varsCurr(strcmp(varsCurr, varRem)) = [];
            end

            if isempty(varsCurr)
                currFrml = sprintf('%s ~ 1%s', varRsp, randStr);
            else
                currFrml = sprintf('%s ~ %s%s', varRsp, strjoin(varsCurr, ' + '), randStr);
            end

            % Fit Model (Use lme_fit)
            mdl = lme_fit(tblTrn, currFrml, 'dist', dist);
            yProb = predict(mdl, tblTst);

            if isBinomial
                % --- CLASSIFICATION METRICS ---
                [x, y, T, aucCurr] = perfcurve(yTrue, yProb, 1);

                % Optimal Threshold (Youden)
                J = y - x;
                [~, idxJ] = max(J);
                optThresh = T(idxJ);
                yPred = double(yProb >= optThresh);

                % Balanced Accuracy
                tp = sum(yTrue == 1 & yPred == 1);
                tn = sum(yTrue == 0 & yPred == 0);
                sens = tp / max(sum(yTrue == 1), eps);
                spec = tn / max(sum(yTrue == 0), eps);
                valPrim = mean([sens, spec]); % BACC
                valSec  = aucCurr;            % AUC

                % ROC Interp
                [uX, idxUni] = unique(x);
                rocCurr = interp1(uX, y(idxUni), xGrid, 'linear', 'extrap');
                res.roc{iVar}(cnt, :) = rocCurr;

            else
                % --- REGRESSION METRICS ---
                % RMSE
                rmse = sqrt(mean((yTrue - yProb).^2));

                % R-Squared
                sst = sum((yTrue - mean(yTrue)).^2);
                sse = sum((yTrue - yProb).^2);
                r2 = 1 - (sse / sst);

                valPrim = rmse;
                valSec  = r2;
            end

            % Store
            res.perfPrim(cnt, iVar) = valPrim;
            res.perfSec(cnt, iVar)  = valSec;
        end
        fprintf('.');
    end
    fprintf('\n');
end

% Calculate Deltas (Interpretation depends on metric)
% For BACC/AUC/R2: Higher is better -> Delta = Full - Reduced (Positive = Important)
% For RMSE: Lower is better -> Delta = Reduced - Full (Positive = Important, removing it increased error)

if strcmp(metricPrim, 'RMSE')
    % RMSE: Importance = Error(Reduced) - Error(Full)
    res.impPrim = res.perfPrim(:, 2:end) - res.perfPrim(:, 1);
else
    % BACC/AUC/R2: Importance = Score(Full) - Score(Reduced)
    res.impPrim = res.perfPrim(:, 1) - res.perfPrim(:, 2:end);
end

% Secondary Metric Deltas (AUC / R-Squared -> Higher is better)
res.impSec = res.perfSec(:, 1) - res.perfSec(:, 2:end);


%% ========================================================================
%  PLOTTING
%  ========================================================================
if flgPlot
    plot_ablation(res, metricPrim, metricSec, isBinomial);
end

end     % EOF


%% ========================================================================
%  HELPER: PLOTTING
%  ========================================================================
function plot_ablation(res, labPrim, labSec, isBinomial)

varsFxd = res.vars(2 : end);
nVars  = length(varsFxd);

figure('Color', 'w', 'Name', 'Feature Ablation', 'Position', [100 100 1200 400]);
t = tiledlayout(1, 3 + isBinomial, 'TileSpacing', 'compact', 'Padding', 'compact');

% 1. Primary Importance
nexttile;
muImp = mean(res.impPrim, 1, 'omitnan');
semImp = std(res.impPrim, 0, 1, 'omitnan') / sqrt(size(res.impPrim, 1));
bar(muImp, 'FaceColor', 'k', 'FaceAlpha', 0.5);
hold on;
errorbar(1:nVars, muImp, semImp, 'k.', 'LineWidth', 1.5);
xticks(1:nVars); xticklabels(varsFxd);
ylabel(['\Delta ' labPrim]);
title(['Importance (' labPrim ')']);
grid on;

% 2. Secondary Importance
nexttile;
muImpSec = mean(res.impSec, 1, 'omitnan');
semImpSec = std(res.impSec, 0, 1, 'omitnan') / sqrt(size(res.impSec, 1));
bar(muImpSec, 'FaceColor', 'b', 'FaceAlpha', 0.5);
hold on;
errorbar(1:nVars, muImpSec, semImpSec, 'k.', 'LineWidth', 1.5);
xticks(1:nVars); xticklabels(varsFxd);
ylabel(['\Delta ' labSec]);
title(['Importance (' labSec ')']);
grid on;

% 3. Model Comparison (Raw Values)
nexttile;
data = res.perfPrim;
labels = res.vars;
ylab = labPrim;
% Inline scatter_raw logic
nCols = size(data, 2);
hold on;
for iCol = 1:nCols
    x = iCol + randn(size(data,1), 1) * 0.05;
    scatter(x, data(:,iCol), 15, 'filled', 'MarkerFaceAlpha', 0.3);
    plot([iCol-0.2, iCol+0.2], [mean(data(:,iCol), 'omitnan'), mean(data(:,iCol), 'omitnan')], 'k-', 'LineWidth', 2);
end
xticks(1:nCols); xticklabels(labels);
ylabel(ylab);
title('Model Performance');
grid on;


% 4. ROC (Binomial Only)
if isBinomial
    nexttile;
    hold on;
    clrs = lines(nVars);
    xGrid = linspace(0, 1, 100);

    % Full
    yMeanFull = mean(res.roc{1}, 1, 'omitnan');
    plot(xGrid, yMeanFull, 'k-', 'LineWidth', 3, 'DisplayName', 'Full');

    for iCol = 1:nVars
        yMeanRed = mean(res.roc{iCol+1}, 1, 'omitnan');
        plot(xGrid, yMeanRed, '-', 'Color', clrs(iCol,:), 'LineWidth', 1.5, ...
            'DisplayName', ['w/o ' varsFxd{iCol}]);
    end
    plot([0 1], [0 1], 'k:', 'HandleVisibility', 'off');
    xlabel('FPR'); ylabel('TPR');
    title('ROC Curves');
    legend('Location', 'best');
end

end     % EOF


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