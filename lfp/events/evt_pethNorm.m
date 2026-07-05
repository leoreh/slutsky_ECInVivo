function pethZ = evt_pethNorm(pethRaw, refData, tstamps, sigma)
% EVT_PETHNORM Smooth and z-score a PETH against a reference distribution.
%
%   pethZ = EVT_PETHNORM(pethRaw, refData, tstamps, sigma)
%
%   SUMMARY:
%       Builds a Gaussian smoothing kernel from the PETH time base, convolves
%       (with edge-effect correction), then z-scores against a reference.
%       When refData has the same number of rows as pethRaw (>1), each row is
%       normalized by its own mean/SD (per-unit). Otherwise refData is pooled
%       into a single distribution (per-population).
%
%   INPUTS:
%       pethRaw - (Mat) [M x nBins] Raw PETH (units or events x bins).
%       refData - (Mat) [M x nBins] or [1 x nBins] Reference PETH data.
%       tstamps - (Vec) [1 x nBins] PETH time base [s] (sets the kernel dt).
%       sigma   - (Num) Gaussian SD [s]. (Default: 0.001).
%
%   OUTPUTS:
%       pethZ   - (Mat) [M x nBins] Smoothed, z-scored PETH.
%
%   DEPENDENCIES:
%       None.
%
%   HISTORY:
%       Created: 05 Jul 2026 (extracted from ripp_wrapper's local peth_norm;
%                the kernel construction, previously inline in the wrapper, now
%                lives here so per-unit and population PETHs share one code path).

if nargin < 4 || isempty(sigma)
    sigma = 0.001;
end

% Gaussian kernel, built from the PETH time base
dt = mode(diff(tstamps));
nSteps = ceil(3 * sigma / dt);
kRng = (-nSteps : nSteps) * dt;
kd = normpdf(kRng, 0, sigma);
kd = kd / sum(kd);

% Edge-effect correction: convolving a unity vector reveals boundary loss
nBins = size(pethRaw, 2);
kCorr = conv(ones(1, nBins), kd, 'same');

% Smooth target
pethSmooth = conv2(pethRaw, kd, 'same') ./ kCorr;

% Reference statistics
if size(refData, 1) == size(pethRaw, 1) && size(refData, 1) > 1
    % Per-unit normalization (reference is Units x Bins)
    refMean = mean(refData, 2, 'omitnan');
    refSd = std(refData, [], 2, 'omitnan');
else
    % Pooled reference (e.g. population mean vector 1 x Bins)
    refMean = mean(refData, 'all', 'omitnan');
    refSd = std(refData, [], 'all', 'omitnan');
end

if any(refSd == 0)
    refSd(refSd == 0) = 1;
end

% Z-score
pethZ = (pethSmooth - refMean) ./ refSd;

end
