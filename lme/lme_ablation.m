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
addParameter(p, 'nReps', 5, @isnumeric);
addParameter(p, 'flgPlot', true, @islogical);
parse(p, tbl, frml, varargin{:});

dist = p.Results.dist;
nReps = p.Results.nReps;
flgPlot = p.Results.flgPlot;

%% ========================================================================
%  PREPARATION
%  ========================================================================

% Run LME_ANALYSE to prepare data and select distribution (Global Check)
fprintf('[LME_ABLATION] Running initial analysis to prepare data...\n');
[~, ~, lmeInfo, tblLme] = lme_analyse(tbl, frml, 'dist', dist);

% Retrieve VIF from lme_analyse output
res.vif = lmeInfo.vif;

% Update dist
dist = lmeInfo.distSelected;
distFinal = lmeInfo.distFinal;

% Determine Mode
flgBinomial = strcmpi(dist, 'Binomial');
if flgBinomial
    metricPrim = 'Balanced Accuracy';
    metricSec  = 'AUC';
    modeStr = 'Classification';
else
    metricPrim = 'RMSE';
    metricSec  = 'R-Squared';
    modeStr = 'Regression';
end
fprintf('[LME_ABLATION] Mode: %s (Dist: %s)\n', modeStr, dist);

% Extract variables from formula and add full model as first iteration
[varsFxd, varRsp, ~, varsIntr] = lme_frml2vars(frml);
varsIter = [{'None'}, varsFxd, varsIntr];


%% ========================================================================
%  MODEL COMPARISON (LIKELIHOOD RATIO TEST)
%  ========================================================================
%  Compare Full Model vs. Reduced Models on the ENTIRE dataset to get
%  statistical significance (p-values) for each predictor.

fprintf('[LME_ABLATION] Running statistical model comparison...\n');

% Determine Fit Method for Comparison
%    - LME: Use 'ML' (Maximum Likelihood) to compare fixed effects.
%    - GLME: Use 'Laplace'
if flgBinomial || ~strcmpi(distFinal, 'Normal')
    % Generalized: Use Laplace for accurate Likelihood approximation
    fitMet = 'Laplace';
else
    % Linear: Use ML to compare models with different Fixed Effects
    fitMet = 'ML';
end

% Initialize Stats Storage
nVars = length(varsIter);
statNames = {'Predictor', 'AIC', 'BIC', 'LogLik', 'LRStat', 'pValue', 'DeltaAIC'};
statsData = cell(nVars, 7);

% Fit Full Model
mdlFull = lme_fit(tblLme, frml, 'dist', distFinal, 'fitMethod', fitMet);

% Iterate Ablations
for iVar = 1:nVars
    varRmv = varsIter{iVar};

    if strcmp(varRmv, 'None')
        % Full Model Baseline
        statsData(iVar, :) = {varRmv, mdlFull.ModelCriterion.AIC, ...
            mdlFull.ModelCriterion.BIC, mdlFull.LogLikelihood, ...
            NaN, NaN, 0};
    else
        % Reduced Model
        currFrml = lme_frml2rmv(frml, varRmv);
        mdlRed = lme_fit(tblLme, currFrml, 'dist', distFinal, 'fitMethod', fitMet);

        % Compare: Reduced (H0) vs Full (H1)
        % Returns table, we extract pValue and LRStat
        cmpTbl = compare(mdlRed, mdlFull);
        lrStat = cmpTbl.LRStat(2);
        pVal   = cmpTbl.pValue(2);
        dAIC   = mdlRed.ModelCriterion.AIC - mdlFull.ModelCriterion.AIC;

        statsData(iVar, :) = {varRmv, mdlRed.ModelCriterion.AIC, ...
            mdlRed.ModelCriterion.BIC, mdlRed.LogLikelihood, ...
            lrStat, pVal, dAIC};
    end
end

% Convert to Table
res.stats = cell2table(statsData, 'VariableNames', statNames);


%% ========================================================================
%  PARTITIONING
%  ========================================================================
% Generate training/testing indices.
% If 'Name' and 'Group' exist, use Stratified Grouped Random Subsampling.
% Otherwise, use standard Stratified K-Fold.

fprintf('[LME_ABLATION] Generating partitions...\n');
matTrn = get_partitions(tblLme, nReps, flgBinomial);

% Storage
nTotal = size(matTrn, 2);
nCols  = length(varsIter);

res.vars        = varsIter;
res.perfPrim    = nan(nTotal, nCols);        % BACC or RMSE
res.perfSec     = nan(nTotal, nCols);        % AUC or R2
res.roc         = cell(nCols, 1);            % ROC (Binomial only)
xGrid           = linspace(0, 1, 100);       % Common grid for ROC

for iCol = 1:nCols
    res.roc{iCol} = nan(nTotal, length(xGrid));
end


%% ========================================================================
%  ABLATION LOOP
%  ========================================================================

fprintf('Starting Ablation (%d Iterations)...\n', nTotal);

for iRep = 1 : nTotal

    % Get indices for this iteration
    trainIdx = matTrn(:, iRep);
    testIdx  = ~trainIdx;

    % Train Data (Raw)
    tblTrn = tbl(trainIdx, :);

    % Iterate Models (Variables to remove)
    for iVar = 1 : length(varsIter)

        varRmv = varsIter{iVar};
        currFrml = lme_frml2rmv(frml, varRmv);

        % Train and Get Params
        try
            [mdl, ~, infoTrn] = lme_analyse(tblTrn, currFrml, ...
                'dist', dist, 'verbose', false);

            % Test: Apply EXACT same transformations as Training
            tblTst = tbl(testIdx, :); % Start with RAW data
            tblTst = tbl_trans(tblTst, 'template', infoTrn.transParams);

            yTrue = tblTst.(varRsp);
            yProb = predict(mdl, tblTst);

        catch ME
            warning('Fit failed for %s (Iter %d): %s', varRmv, iRep, ME.message);
            yProb = nan(size(yTrue)); % Ensure yProb is defined for metrics
        end

        % --- METRICS CALCULATION ---
        if flgBinomial      % CLASSIFICATION

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
            res.roc{iVar}(iRep, :) = rocCurr;

        else                % REGRESSION

            % Back transform response values
            pVar = infoTrn.transParams.varsTrans.(varRsp);
            if ~isempty(pVar.logBase)
                c = pVar.offset;
                yTrue = exp(yTrue) - c;
                yProb = exp(yProb) - c;
            end

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
        res.perfPrim(iRep, iVar) = valPrim;
        res.perfSec(iRep, iVar)  = valSec;
    end

    % Progress dot
    if mod(iRep, 2) == 0
        fprintf('.');
    end
end
fprintf('\n');

% Calculate Deltas (Percentage Change relative to Full Model)
% Index 1 is always 'None' (Full Model)
basePrim = res.perfPrim(:, 1);
baseSec  = res.perfSec(:, 1);

if strcmp(metricPrim, 'RMSE')
    % RMSE: importance = % Increase in Error (Lower is better)
    % (Reduced - Full) / Full * 100
    res.impPrim = (res.perfPrim(:, 2:end) - basePrim) ./ basePrim * 100;
else
    % BACC/AUC/R2: importance = % Drop in Performance (Higher is better)
    % (Full - Reduced) / Full * 100
    res.impPrim = (basePrim - res.perfPrim(:, 2:end)) ./ basePrim * 100;
end

% Secondary Metric Deltas (AUC / R2 -> Higher is better)
res.impSec = (baseSec - res.perfSec(:, 2:end)) ./ baseSec * 100;


%% ========================================================================
%  PLOTTING
%  ========================================================================
if flgPlot
    plot_ablation(res, metricPrim, metricSec, flgBinomial);
end

end     % EOF


%  ========================================================================
%  HELPER: PARTITIONS
%  ========================================================================
function matTrn = get_partitions(tbl, nReps, flgBinomial)

% Detect Grouped Mode
hasName  = ismember('Name', tbl.Properties.VariableNames);
hasGroup = ismember('Group', tbl.Properties.VariableNames);

if hasName && hasGroup
    % --- Exhaustive Leave-One-Pair-Out (Deterministic) ---
    % 1. Identify all unique subjects and their groups.
    % 2. Generate every possible combination of ONE subject from EACH group.
    %    (e.g., 5 Ctrl * 5 KO = 25 unique test sets).
    % 3. Each combination becomes a Test Fold.
    % Note: Only works for two groups currently

    [uNames, idxFirst] = unique(tbl.Name, 'stable');
    nameGrps = tbl.Group(idxFirst);
    uGrps = unique(nameGrps);

    % Collect subjects by group
    grpSubjects = cell(numel(uGrps), 1);
    for iGrp = 1:numel(uGrps)
        grpSubjects{iGrp} = uNames(nameGrps == uGrps(iGrp));
    end

    % Generate Cartesian Product of subjects (1 from each group)
    for iGrp = 2:numel(uGrps)
        currSub = grpSubjects{iGrp};
        nOld = size(grpSubjects{1}, 1);
        nNew = numel(currSub);

        tmpOld = repmat(grpSubjects{1}, nNew, 1);
        tmpNew = repelem(currSub, nOld, 1);
        combos = [tmpOld, tmpNew];
    end

    nTotal = size(combos, 1);
    fprintf('[LME_ABLATION] Found %d unique test combinations.\n', nTotal);

    matTrn = true(height(tbl), nTotal);

    for iRep = 1 : nTotal
        testNames = combos(iRep, :); % 1xGroups

        % Defines TEST set (Leave-These-Out)
        matTrn(:, iRep) = ~ismember(tbl.Name, testNames);
    end

else
    % --- Standard Repeated K-Fold ---
    nFolds = 5;
    nTotal = nReps * nFolds;
    matTrn = true(height(tbl), nTotal);

    cnt = 0;
    for iRep = 1 : nReps
        if flgBinomial
            cvp = cvpartition(tbl{:, 1}, 'KFold', nFolds);
        else
            cvp = cvpartition(height(tbl), 'KFold', nFolds);
        end

        for iFold = 1 : nFolds
            cnt = cnt + 1;
            matTrn(:, cnt) = training(cvp, iFold);
        end
    end
end

end

%% ========================================================================
%  HELPER: PLOTTING
%  ========================================================================
function plot_ablation(res, labPrim, labSec, flgBinomial)

varsFxd = res.vars(2 : end);
nVars  = length(varsFxd);

figure('Color', 'w', 'Name', 'Feature Ablation', 'Position', [100 100 1200 400]);
t = tiledlayout(1, 3 + flgBinomial, 'TileSpacing', 'compact', 'Padding', 'compact');

% 1. Primary Importance
nexttile;
muImp = mean(res.impPrim, 1, 'omitnan');
semImp = std(res.impPrim, 0, 1, 'omitnan') / sqrt(size(res.impPrim, 1));
bar(muImp, 'FaceColor', 'k', 'FaceAlpha', 0.5);
hold on;
errorbar(1:nVars, muImp, semImp, 'k.', 'LineWidth', 1.5);
xticks(1:nVars); xticklabels(varsFxd);
ylabel(['% \Delta ' labPrim]);
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
ylabel(['% \Delta ' labSec]);
title(['Importance (' labSec ')']);
grid on;

% 3. Model Comparison (Raw Values)
nexttile;
data = res.perfPrim;
labels = res.vars;
ylab = labPrim;
% Connecting Lines for paired observations
plot(1:size(data,2), data', '-', 'Color', [0.7 0.7 0.7 0.3], 'LineWidth', 0.5);
hold on;

% Inline scatter_raw logic (Means + Scatters)
nCols = size(data, 2);
for iCol = 1:nCols
    x = iCol + randn(size(data,1), 1) * 0.05;
    scatter(x, data(:,iCol), 15, 'filled', 'MarkerFaceAlpha', 0.8);
    plot([iCol-0.2, iCol+0.2], [mean(data(:,iCol), 'omitnan'), mean(data(:,iCol), 'omitnan')], 'k-', 'LineWidth', 2);
end
xticks(1:nCols); xticklabels(labels);
ylabel(ylab);
title('Model Performance');
grid on;

% 4. ROC (Binomial Only)
if flgBinomial
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
%  NOTE: METRICS FOR CONTINUOUS SCALES
%  ========================================================================
%  Transitioning from binary classification (e.g., Recovery T/F) to
%  continuous regression (e.g., Recovery Time, Firing Rate) requires
%  different metrics. "Accuracy" and "ROC" are invalid because there is no
%  categorical threshold.
%
%  1. Primary Metric: Coefficient of Determination (R-Squared)
%     Represents the proportion of variance explained by the model.
%     - Marginal R2: Explained by Fixed Effects only (e.g., predictors).
%     - Conditional R2: Explained by Fixed + Random Effects (e.g., animal).
%
%  2. Error Metrics (RMSE / MAE)
%     Quantifies the absolute error in physiological units (e.g., Hz).
%     - RMSE (Root Mean Square): Penalizes large outliers heavily.
%     - MAE (Mean Absolute): Linear penalty, more robust to outliers.
%
%  3. Ablation Analysis (Delta R2)
%     To rank feature importance, calculate "Delta R2":
%        Delta = R2_Full_Model - R2_Reduced_Model
%     - Large Positive Delta: Variable is a critical driver.
%     - Near Zero Delta: Variable is redundant.
%
%  4. Model Selection (AIC / BIC)
%     Used to balance goodness-of-fit with complexity.
%     - Penalizes overfitting (adding variables that don't help much).
%     - Lower is better.
%
%  5. Standardized Coefficients (Beta)
%     If predictors are Z-scored, Beta represents the "Effect Size" in
%     standard deviations.
%     - Beta = 0.5 -> 1 SD increase in X leads to 0.5 SD increase in Y.
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


%% ========================================================================
%  NOTE: STRATIFIED GROUPED K-FOLD
%  ========================================================================
%  In hierarchical physiological datasets (e.g., units nested within
%  cultures/animals), standard cross-validation strategies often lead to
%  two critical failures that invalidate ablation results:
%
%  1. THE PROBLEM: DATA LEAKAGE VS. MATRIX SINGULARITY
%     - Standard CV: Splitting units from the same 'Name' between training
%       and testing violates the independence assumption of the (1|Name)
%       random effect. This results in "Data Leakage,"
%       where the model overestimates accuracy by 'memorizing' specific
%       culture baselines.
%     - Grouped CV: Simply grouping by 'Name' ensures units stay together
%       but often creates "Empty Folds". If a training fold
%       randomly lacks any representatives of a specific 'Group' (e.g., all
%       MCU-KO names end up in the test set), the model matrix becomes
%       singular. This causes catastrophic numerical
%       failures (RMSE > 10^75).
%
%  2. THE SOLUTION: MANUAL STRATIFIED GROUPED PARTITIONING
%     The most rigorous solution is to partition the 'Biological Units'
%     (the Names/Cultures) rather than the table rows, while simultaneously
%     stratifying by the experimental condition (Group). This ensures that
%     every training fold contains a balanced ratio of Control and KO
%     cultures, guaranteeing mathematical stability for fixed-effect
%     estimation.
%
%  3. COMBINATORICS AND REPETITIONS
%     With 5 Control and 5 MCU-KO cultures (N=10), there are exactly 25
%     unique combinations of (1 Ctrl + 1 KO) test pairs. Because the
%     training set is defined as "everyone else," each test pair is linked
%     to a unique, deterministic training set.
%
%     - Exhaustive Search: Instead of random shuffling, we generate all 25
%       combinations (Cartesian product) to serve as the test sets.
%     - Completeness: This removes sampling noise and ensures the results
%       represent the true, complete distribution of the model's predictive
%       power across the cohort.
%
%  4. TECHNICAL IMPLEMENTATION
%     - Step A: Extract unique Names and their Group assignments.
%     - Step B: specific 'cvpartition' is NOT used. Instead, we generate
%       indices for every possible combination of one-subject-per-group.
%     - Step C: The number of iterations is automatically set to the total
%       number of combinations (e.g., 25), ignoring 'nReps'.
%
%  This approach is optimal for testing if predictors (pBspk, fr, etc.)
%  can generalize to a completely new, unseen animal or culture.
%  ========================================================================

%% ========================================================================
%  NOTE: NUMERICAL ROBUSTNESS VS. MODEL FIT
%  ========================================================================
%  In hierarchical modeling, a conflict often arises between the "Best
%  Fitting" model (lowest AIC) and the "Most Robust" model for
%  cross-validation. While AIC measures relative information loss on a
%  static dataset, it does not guarantee the numerical stability of the
%  Iteratively Reweighted Least Squares (IRLS) algorithm when data is
%  partitioned.
%
%  1. The Gamma/Log-Link Sensitivity:
%     The Gamma distribution with a Log link is frequently selected by AIC
%     for firing rate data. However, this combination is highly sensitive
%     to data perturbations. In Leave-One-Experiment-Out (LOEO) schemes,
%     removing a single culture can leave the training set with an
%     unconstrained slope. This causes the pseudo-likelihood weights to
%     explode, leading to numerical model collapse (RMSE > 10^70).
%
%  2. The Log-Normal Stability:
%     Modeling log-transformed data using a Gaussian distribution (Log-Normal)
%     is often the superior choice for ablation analysis.
%     - Uses closed-form mathematical solutions.
%     - Avoids iterative weight scaling issues.
%     - Provides stable out-of-sample predictions.
%     - Aligns with neurophysiological theory.
%
%  3. Justification for Overriding AIC:
%     An AIC-driven selection should be overridden when the primary
%     analytical goal is feature ranking rather than population
%     description. Numerical stability is a prerequisite for interpretable
%     importance metrics. If a Gamma model fails to converge in any
%     cross-validation fold, its global AIC advantage is invalidated for
%     that specific task.
%
%  4. Practical Recommendations:
%     - Report p-values using Gamma.
%     - Use Log-Normal for importance.
%     - Verify convergence in every fold.
%     - Z-score inside the CV loop.
%     - Monitor for badly scaled weights.
%
%  By prioritizing numerical robustness during ablation, you ensure that
%  "Importance" reflects biological signal rather than the failure
%  of a specific link function's convergence algorithm.
%  ========================================================================


%% ========================================================================
%  NOTE: LOG VS. LINEAR PERFORMANCE METRICS
%  ========================================================================
%  When performing ablation analysis on variables that require
%  log-transformation (e.g., Firing Rate, which is often Log-Normal), the
%  choice of scale for metric calculation (RMSE, R-Squared) fundamentally
%  alters the "importance" ranking of predictors.
%
%  1. LOG-SPACE METRICS (RELATIVE ERROR):
%     Calculating RMSE/R2 on log-transformed values treats errors as
%     proportional changes (fold-changes).
%     - Impact: An error in predicting 0.1 Hz vs 0.2 Hz is penalized as
%       heavily as an error in predicting 10 Hz vs 20 Hz.
%     - Interpretation: This "democratizes" the data, giving equal weight
%       to low-firing and high-firing units. Importance reflects the
%       model's ability to capture the underlying biological mechanism
%       across the entire population range.
%
%  2. LINEAR-SPACE METRICS (ABSOLUTE ERROR):
%     Back-transforming (exp(y)) before calculating metrics shifts the
%     focus to absolute magnitudes.
%     - Impact: An error of 10 Hz on a high-firing unit becomes
%       catastrophic to the RMSE, while a 100% error on a 0.1 Hz unit
%       (Error = 0.1) becomes mathematically negligible.
%     - Interpretation: This "oligarchic" approach allows high-magnitude
%       outliers to dominate the error sum. Importance reflects which
%       variables are best at predicting the "peaks" or high-value states.
%
%  3. WHY FEATURE IMPORTANCE SHIFTS:
%     If a variable (e.g., 'pBspk') loses importance in linear-space but
%     was critical in log-space, it means that variable was primarily
%     useful for fine-tuning the predictions of low-magnitude units.
%     Conversely, variables that gain importance in linear-space are the
%     primary drivers of the highest-magnitude observations.
%
%  RECOMMENDATION:
%  - Use LOG-metrics to characterize biological "States" and "Rules."
%  - Use LINEAR-metrics to characterize "Real-world" predictive accuracy
%    in physiological units (Hz).
% ========================================================================


%% ========================================================================
%  NOTE: WALD TEST VS. LIKELIHOOD RATIO TEST
%  ========================================================================
%  The p-values from `lme_analyse` (Wald Test) and `lme_ablation`
%  (Likelihood Ratio Test) are nearly identical. This is a fundamental
%  property of statistical theory, not a redundancy error.
%
%  1. ASYMPTOTIC EQUIVALENCE:
%     - The Wald Test estimates the curvature of the Likelihood surface at
%       the peak (Maximum Likelihood Estimate) assuming a perfect quadratic
%       shape.
%     - The LRT measures the actual drop in height of the Likelihood
%       surface when a variable is removed.
%     - For large sample sizes and "well-behaved" data, the Likelihood
%       surface approximates a perfect parabola. In this scenario, the
%       curvature at the peak (Wald) perfectly predicts the drop in height
%       (LRT), causing the p-values to converge.
%
%  2. WHY USE ABLATION (LRT) IF THEY ARE EQUAL?
%     Despite the numerical similarity, the ablation loop provides distinct
%     advantages required for rigorous model selection:
%
%     - Metric of Information (AIC/BIC):
%       Wald p-values only indicate non-zero significance. They do not
%       quantify *how much* information a variable adds. Ablation allows
%       you to rank variables by Delta-AIC (e.g., Variable A adds 50 bits
%       of information vs. Variable B's 2 bits).
%
%     - Categorical Variables (ANOVA):
%       `lme_analyse` outputs separate p-values for every level of a
%       Group (e.g., Level 2 vs 1, Level 3 vs 1). Ablation removes the
%       entire factor to provide a single, holistic p-value for "Group"
%       existence.
%
%     - Robustness & Stability:
%       Identical results are a diagnostic "Health Check." If Wald and LRT
%       diverged significantly, it would indicate a non-quadratic
%       likelihood surface, separation issues (in Binomial models), or
%       small-sample bias. Their convergence confirms your model is stable.
%
%  3. SUMMARY:
%     "Identical" results confirm the model is mathematically robust. Use
%     Wald for quick significance checks; use Ablation for information
%     ranking (AIC) and variable importance.
%  ========================================================================