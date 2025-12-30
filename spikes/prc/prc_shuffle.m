function shuffledMat = prc_shuffle(rasterMat, nSwaps)
% PRC_SHUFFLE Shuffles binary spike matrix while preserving marginal sums.
% REPLACED WITH MEX VERSION
%
%   shuffledMat = PRC_SHUFFLE(RASTERMAT, NSWAPS) performs a controlled shuffle
%   of the binary spike matrix using the "raster marginals model" (Okun et al., 2015).
%
%   INPUTS:
%       rasterMat   - (bin) Binary matrix (neurons x time_bins).
%       nSwaps      - (num) Number of swap attempts {nBins * nUnits}.
%
%   OUTPUTS:
%       shuffledMat - (bin) Shuffled matrix with preserved marginal sums.
%
%   See also: PRC_CALC
%
%   NOTE: METHODOLOGY
%   This function swaps 2x2 submatrices to preserve both row (neuron rate)
%   and column (population activity) sums:
%       [[1, 0]; [0, 1]]  <-->  [[0, 1]; [1, 0]]

%% ========================================================================
%  ARGUMENTS
%  ========================================================================
[nUnits, nBins] = size(rasterMat);
if nargin < 2 || isempty(nSwaps)
    nSwaps = nBins * nUnits;
end

%% ========================================================================
%  INITIALIZE
%  ========================================================================
shuffledMat = rasterMat;

% Pre-generate all random indices for rows and columns
r1_all = randi(nUnits, nSwaps, 1);
r2_all = randi(nUnits-1, nSwaps, 1);
r2_all = r2_all + (r2_all >= r1_all);  % Ensure f2 != f1
c1_all = randi(nBins-1, nSwaps, 1);
c2_all = randi(nBins-1, nSwaps, 1);
c2_all = c2_all + (c2_all >= c1_all);  % Ensure t2 != t1

%% ========================================================================
%  SWAP LOOP
%  ========================================================================
for iSwap = 1:nSwaps
    % Get current indices
    r1 = r1_all(iSwap);
    r2 = r2_all(iSwap);
    c1 = c1_all(iSwap);
    c2 = c2_all(iSwap);

    % Get values directly
    v11 = shuffledMat(r1, c1);
    v12 = shuffledMat(r1, c2);
    v21 = shuffledMat(r2, c1);
    v22 = shuffledMat(r2, c2);

    % Swap if pattern matches [1,0;0,1] or [0,1;1,0]
    if (v11 == 1 && v12 == 0 && v21 == 0 && v22 == 1)
        shuffledMat(r1, c1) = 0;
        shuffledMat(r1, c2) = 1;
        shuffledMat(r2, c1) = 1;
        shuffledMat(r2, c2) = 0;
    elseif (v11 == 0 && v12 == 1 && v21 == 1 && v22 == 0)
        shuffledMat(r1, c1) = 1;
        shuffledMat(r1, c2) = 0;
        shuffledMat(r2, c1) = 0;
        shuffledMat(r2, c2) = 1;
    end
end

end     % EOF
