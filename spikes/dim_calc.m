function dim = dim_calc(Y, varargin)
% DIM_CALC Computes the effective dimensionality of neural activity.
%
%   dim = DIM_CALC(Y, ...) calculates the dimensionality of the population
%   activity matrix Y (Neurons x Time).
%
%   INPUTS:
%       Y           - (matrix) Activity Matrix (N x T). See Note.
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
addRequired(p, 'Y', @isnumeric);
addParameter(p, 'method', 'PR', @ischar);
addParameter(p, 'thrVal', 0.8, @isnumeric);
addParameter(p, 'flgFrac', false, @islogical);

parse(p, Y, varargin{:});
method = p.Results.method;
thrVal = p.Results.thrVal;
flgFrac = p.Results.flgFrac;


%% ========================================================================
%  PRE-PROCESSING
%  ========================================================================

% Remove silent neurons (rows with no activity)
mY = mean(Y, 2);
silentIdx = mY == 0;
Y(silentIdx, :) = [];

[nUnits, ~] = size(Y);

% Safety check for low dimensionality
if nUnits < 2
    dim = 0;
    return;
end


%% ========================================================================
%  CALC
%  ========================================================================

[~, ~, eigVals] = pca(Y);

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
% Given an activity matrix Y (Neurons x Time), the function computes the
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

