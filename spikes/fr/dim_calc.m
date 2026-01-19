function dim = dim_calc(frMat, varargin)
% DIM_CALC Computes the effective dimensionality of neural activity.
%
%   dim = DIM_CALC(frMat, ...) calculates the dimensionality of the population
%   activity matrix frMat (Neurons x Time).
%
%   INPUTS:
%       frMat       - (matrix) Activity Matrix (Neurons x Time).
%       varargin    - (param/value) Optional parameters:
%                     'method'     : (char) 'PR' (Participation Ratio) or
%                                    'thr' (Variance Explained) {'PR'}
%                     'thrVal'     : (num) Var explained cutoff (0-1) {0.8}
%                     'flgFrac'    : (log) Normalize by N neurons {false}
%
%   OUTPUTS:
%       dim         - (scalar) Effective dimensionality.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'frMat', @isnumeric);
addParameter(p, 'method', 'PR', @ischar);
addParameter(p, 'thrVal', 0.8, @isnumeric);
addParameter(p, 'flgFrac', false, @islogical);
addParameter(p, 'flgPlot', false, @islogical);

parse(p, frMat, varargin{:});
method = p.Results.method;
thrVal = p.Results.thrVal;
flgFrac = p.Results.flgFrac;
flgPlot = p.Results.flgPlot;


%% ========================================================================
%  PRE-PROCESSING
%  ========================================================================

% Remove silent neurons (rows with no activity)
mY = mean(frMat, 2);
silentIdx = mY == 0;
frMat(silentIdx, :) = [];

% Assert frMat orientation. Heurestic; more observations then variables
[nTime, nUnits] = size(frMat);
if nUnits > nTime
    frMat = frMat';
    [nTime, nUnits] = size(frMat);
end

% Safety check for low dimensionality
if nUnits < 2
    dim = 0;
    return;
end


%% ========================================================================
%  CALC
%  ========================================================================

[~, ~, eigVals] = pca(frMat);

switch lower(method)
    case 'pr'       % Participation Ratio (PR)

        dim = sum(eigVals)^2 / sum(eigVals.^2);

    case 'thr'      % Variance Explained Threshold

        pSpec = eigVals / sum(eigVals);
        cumVar = cumsum(pSpec);
        dim = find(cumVar > thrVal, 1, 'first');
end

if flgFrac
    dim = dim / nUnits;
end


%% ========================================================================
%  VISUALIZATION
%  ========================================================================

if flgPlot

    % Re-calculate PCA for visualization (scores needed)
    [~, score, ~, ~, explained] = pca(frMat);

    % Prepare figure
    figure('Color', 'w', 'Name', 'Dimensionality Analysis', ...
        'NumberTitle', 'off', 'Position', [100, 100, 1200, 900]);
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % ---------------------------------------------------------------------
    % Scree Plot (Pareto)
    % ---------------------------------------------------------------------
    nexttile;
    hold on;
    yyaxis left
    bar(1:min(nUnits, 20), explained(1:min(nUnits, 20)), ...
        'FaceColor', [0.2 0.2 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    ylabel('Variance Explained (%)');

    yyaxis right
    plot(cumsum(explained), 'k-o', 'LineWidth', 2, 'MarkerFaceColor', 'w');
    yline(thrVal*100, 'r--', sprintf('Threshold (%.0f%%)', thrVal*100), ...
        'LabelHorizontalAlignment', 'left', 'LineWidth', 1.5);
    ylabel('Cumulative Variance (%)');

    xlabel('Principal Component');
    title({'Scree Plot', sprintf('Estimated Dim: %.2f', dim)});
    grid on;
    box on;
    xlim([0.5, min(nUnits, 20) + 0.5]);

    % ---------------------------------------------------------------------
    % Manifold Projection (3D)
    % ---------------------------------------------------------------------
    nexttile;
    % Create time vector for coloring
    timeColors = parula(nTime);

    plot3(score(:,1), score(:,2), score(:,3), 'Color', [0.8 0.8 0.8 0.5]); % Trace
    hold on;
    scatter3(score(:,1), score(:,2), score(:,3), 20, timeColors, 'filled');

    xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
    ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
    zlabel(sprintf('PC3 (%.1f%%)', explained(3)));
    title('Neural Manifold (First 3 PCs)');
    grid on;
    view(3);
    axis tight;

    % ---------------------------------------------------------------------
    % Saturation Curve (Subsampling)
    % ---------------------------------------------------------------------
    nexttile;

    nSteps = 10;
    subSampleSizes = unique(round(linspace(2, nUnits, nSteps)));
    nBoots = 5;
    dimCurve = zeros(length(subSampleSizes), 2); % Mean, Std

    for iSub = 1:length(subSampleSizes)
        k = subSampleSizes(iSub);
        currDims = zeros(nBoots, 1);

        for iBoot = 1:nBoots
            % Randomly select k neurons
            randIdx = randperm(nUnits, k);
            subMat = frMat(randIdx, :);

            % Recalculate dim for this subset (using same method)
            % Pass flgPlot=false to avoid recursion/plotting
            currDims(iBoot) = dim_calc(subMat, 'method', method, 'thrVal', thrVal, ...
                'flgFrac', flgFrac, 'flgPlot', false);
        end
        dimCurve(iSub, 1) = mean(currDims);
        dimCurve(iSub, 2) = std(currDims);
    end

    errorbar(subSampleSizes, dimCurve(:,1), dimCurve(:,2), 'k-o', ...
        'LineWidth', 2, 'MarkerFaceColor', 'w', 'CapSize', 0);
    xlabel('Number of Neurons');
    ylabel('Estimated Dimensionality');
    title('Saturation Curve (Subsampling)');
    grid on;
    axis tight;

    % ---------------------------------------------------------------------
    % Raw Data Correlation
    % ---------------------------------------------------------------------
    nexttile;

    % Calculate correlation matrix
    C = corr(frMat);
    imagesc(C);
    colormap(gca, jet);
    colorbar;
    clim([-1 1]);

    title('Neuron Correlation Matrix');
    xlabel('Neuron ID');
    ylabel('Neuron ID');
    axis square;

    sgtitle('Dimensionality Analysis Summary', 'FontSize', 14, 'FontWeight', 'bold');
end

end     % EOF





%% ========================================================================
%  NOTE: NEURAL DIMENSIONALITY
%  ========================================================================
% This function estimates the effective dimensionality of neural population
% activity using two primary metrics: the Participation Ratio (PR) and the
% Variance Explained Threshold. These methods help distinguish between the
% "structural" dimensionality (number of neurons recorded) and the
% "functional" dimensionality (the subspace containing the majority of
% signal variance).
%
% THEORETICAL BASIS: PCA AND EIGENVALUES
% Given an activity matrix frMat (Neurons x Time), the function computes the
% eigenvalues (lambda) of the covariance matrix. These eigenvalues
% represent the variance explained by each principal component. If neural
% activity is highly correlated, a few eigenvalues will be large while the
% rest are near zero, indicating low effective dimensionality. Conversely,
% if neurons are independent, the eigenvalues will be more uniform,
% indicating high dimensionality.
%
% PARTICIPATION RATIO (PR)
% The 'PR' method provides a continuous estimate of dimensionality that is
% less sensitive to arbitrary thresholds. It is defined mathematically as:
%                       sum(lambda)^2 / sum(lambda^2)
%
% If all variance is concentrated in a single dimension, PR = 1. If
% variance is spread equally across all  dimensions, PR = N. In neural
% data, PR often yields a non-integer value that represents the "effective"
% number of dimensions. It is particularly useful because it accounts for
% the entire distribution of the eigenvalue spectrum (the "spectral
% probability density") rather than just a subset.
%
% VARIANCE EXPLAINED THRESHOLD
% The 'thr' method calculates the minimum number of principal
% components required to account for a user-defined fraction of the total
% variance (defaulting to 80%). This is a discrete measure:
% 1. Eigenvalues are normalized to sum to 1, creating a cumulative variance
%    distribution.
% 2. The function finds the smallest integer k such that the sum of the
%    top k eigenvalues exceeds the thr.
% While intuitive, this method is sensitive to the chosen thr and
% may overemphasize dimensions with very low signal-to-noise ratios if the
% thr is set too high.
%
% DATA PRE-PROCESSING: SILENT NEURONS
% The function automatically identifies and removes "silent neurons" (rows
% where the mean activity is zero). This is a critical step because
% non-responsive neurons contribute zero variance; leaving them in the
% matrix would artificially inflate the neuron count (N) used for
% normalization (flgFrac) without contributing to the underlying manifold,
% leading to inaccurate "fractional dimensionality" scores.
%
% FRACTIONAL DIMENSIONALITY (flgFrac)
% When 'flgFrac' is set to true, the output is normalized by the number
% of active neurons (N). This produces a value between 0 and 1,
% representing the "compression" of the neural code. A low fractional
% dimensionality suggests that the population activity is highly
% constrained to a low-dimensional manifold, a common observation in
% motor cortex and structured decision-making tasks.
%  ========================================================================

