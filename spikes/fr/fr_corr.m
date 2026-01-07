function [mcc, cc] = fr_corr(Y, varargin)
% FR_CORR Computes pairwise correlations of single-unit FR trajectories.
%
%   [mcc, cc] = FR_CORR(Y) calculates the Pearson correlation
%   coefficients between all pairs of single-unit FR trajectories.
%
%   The methodology is as follows:
%   1. Omit cells that did not fire (silent neurons).
%   2. Z-score each unit (subtract time-average FR and divide by temporal STD).
%   3. Compute temporal covariance matrix on this data.
%   4. Subtract the diagonal (since it is all ones).
%   5. Compute average of all elements of this matrix.
%
%   INPUTS:
%       Y           - (matrix) Activity Matrix (Neurons x Time).
%       varargin    - (param/value) Optional parameters:
%                     'flgPlot' : (log) Plot correlation matrix {false}
%
%   OUTPUTS:
%       mcc         - (scalar) Average correlation (of all elements).
%       cc          - (matrix) Pairwise correlation matrix (N x N) with
%                              diagonal subtracted (zeros).
%
%   See also DIM_CALC


%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'Y', @isnumeric);
addParameter(p, 'flgPlot', false, @islogical);

parse(p, Y, varargin{:});
flgPlot = p.Results.flgPlot;


%% ========================================================================
%  PRE-PROCESSING
%  ========================================================================

% Remove silent neurons (rows with no activity) for calculation
% But keep track of indices to maintain output size
nUnits = size(Y, 1);
mY = mean(Y, 2);
validIdx = mY > 0;

% Safety check
if sum(validIdx) < 2
    mcc = NaN;
    cc = nan(nUnits, nUnits);
    return;
end

% zscore only valid units
Z = zscore(Y(validIdx, :), 0, 2);


%% ========================================================================
%  CALC
%  ========================================================================

% cov expects observations in rows and variables in columns.
% We want correlation between neurons (variables), so transpose.
ccValid = cov(Z');
ccValid = ccValid - diag(diag(ccValid));

% Calculate mean correlation from the valid sub-matrix
mcc = mean(ccValid(:));

% Reconstruct full matrix
cc = nan(nUnits, nUnits);
cc(validIdx, validIdx) = ccValid;


%% ========================================================================
%  PLOT
%  ========================================================================

if flgPlot
    figure;
    imagesc(cc);
    colorbar;
    axis square;
    title(['Mean Corr: ' num2str(mcc, '%.3f')]);
    xlabel('Neuron ID');
    ylabel('Neuron ID');
    set(gca, 'TickDir', 'out');
end

end     % EOF



%% ========================================================================
%  NOTE: CORRELATION METHODOLOGY
%  ========================================================================
% This function implements a specific pairwise correlation methodology for
% analyzing neural population coupling.
%
% METHODOLOGY STEPS:
% 1. Silent Cell Removal: Neurons with zero firing rate in the window are
%    excluded to prevent division by zero during z-scoring and to focus
%    analysis on active participants.
% 2. Z-Scoring: Normalizes the firing rate trajectories of each neuron.
%    This is crucial because the raw firing rates may differ significantly
%    between neurons. Z-scoring converts activity into units of standard
%    deviation from the mean, emphasizing temporal fluctuations rather than
%    absolute rate.
% 3. Covariance Calculation: Computing the covariance of z-scored variables
%    is mathematically equivalent to computing the Pearson correlation
%    coefficient of the original variables.
%    Cov(Z_x, Z_y) = E[(X - mu_x)/sigma_x * (Y - mu_y)/sigma_y] = rho_xy
% 4. Diagonal Subtraction: The auto-correlation (diagonal) is always 1.
%    Removing it ensures that the "average correlation" metric reflects
%    cross-neuronal coupling, not self-correlation.
% 5. Averaging: The final metric collapses the entire unit-unit interaction
%    space into a single scalar representing the global synchrony or coupling
%    state of the network.
%  ========================================================================
