function [thrCal, nf] = ripp_noiseFloor(lfp, fs, method, varargin)
% RIPP_NOISEFLOOR Calibrate a detection threshold to the 1/f noise floor.
%
%   [thrCal, nf] = RIPP_NOISEFLOOR(lfp, fs, method, varargin)
%
%   SUMMARY:
%       Sets a ripple detection threshold from the recording's own
%       aperiodic (1/f) background, so that background alone yields at most
%       a target false-positive rate. The steps mirror van Schalkwijk &
%       Helfrich 2026: estimate the spectral exponent by a straight-line
%       fit of the power spectrum in log-log space below the ripple band;
%       synthesize a matched 1/f surrogate with no oscillation; run the
%       identical detector (ripp_sigPrep + ripp_times) on the surrogate;
%       and return the peak threshold (in SD) at which the surrogate's
%       event rate first drops to the target. A recording whose noise sits
%       in the detector's sensitive range is assigned a higher threshold,
%       which equalizes the false-positive rate across recordings rather
%       than the raw SD.
%
%   INPUTS:
%       lfp      - <vec>    raw detection-channel LFP (microvolts).
%       fs       - <num>    sampling frequency (Hz).
%       method   - <struct> one ripp_methods element; uses .passband,
%                           .detectMet, .thr (the gap between start and peak
%                           thresholds is preserved).
%       varargin - Parameter/Value:
%           'targetFP'  - <num> target surrogate rate (events/s). {0.05}
%           'fitRange'  - <vec> [lo hi] Hz for the log-log slope. {[20 45]}
%           'nremTimes' - <mat> [N x 2] NREM bouts (s); restricts the exponent
%                               fit to background epochs. {[]}
%           'durSurr'   - <num> surrogate duration (s). {600}
%           'seed'      - <num> rng seed for a reproducible surrogate. {0}
%
%   OUTPUTS:
%       thrCal   - <num>    calibrated peak threshold (SD). Falls back to
%                           method.thr(2) if the fit or sweep fails.
%       nf       - <struct> .chi .targetFP .fitRange .rate (per grid point)
%                           .thrGrid .thrCal .rateCal .r2 (slope fit quality).
%
%   DEPENDENCIES:
%       ripp_sigPrep, ripp_times.
%
%   HISTORY:
%       260717 native 1/f noise-floor threshold calibration.

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'lfp', @isnumeric);
addRequired(p, 'fs', @isnumeric);
addRequired(p, 'method', @isstruct);
addParameter(p, 'targetFP', 0.05, @isnumeric);
addParameter(p, 'fitRange', [20 45], @isnumeric);
addParameter(p, 'nremTimes', [], @isnumeric);
addParameter(p, 'durSurr', 600, @isnumeric);
addParameter(p, 'seed', 0, @(x) isempty(x) || isscalar(x));
parse(p, lfp, fs, method, varargin{:});
targetFP  = p.Results.targetFP;
fitRange  = p.Results.fitRange;
nremTimes = p.Results.nremTimes;
durSurr   = p.Results.durSurr;
seed      = p.Results.seed;

lfp    = double(lfp(:));
thrPk  = method.thr(2);
thrGap = method.thr(2) - method.thr(1);     % keep start below peak by this

nf = struct('chi', NaN, 'targetFP', targetFP, 'fitRange', fitRange, ...
    'rate', [], 'thrGrid', [], 'thrCal', thrPk, 'rateCal', NaN, 'r2', NaN);
thrCal = thrPk;

%% ========================================================================
%  APERIODIC EXPONENT (log-log slope below the ripple band)
%  ========================================================================

% restrict the spectrum estimate to background (NREM) epochs when given
sigBg = lfp;
if ~isempty(nremTimes)
    idx = bouts2idx(nremTimes, fs, numel(lfp));
    if numel(idx) > 4 * fs
        sigBg = lfp(idx);
    end
end

% Welch spectrum, then a straight line to log10(PSD) over fitRange
win = round(2 * fs);
if numel(sigBg) < 2 * win
    return                          % too little data; keep the fixed threshold
end
[pxx, faxis] = pwelch(sigBg, hann(win), round(win / 2), [], fs);
mask = faxis >= fitRange(1) & faxis <= fitRange(2) & pxx > 0;
if nnz(mask) < 5
    return
end

xLog = log10(faxis(mask));
yLog = log10(pxx(mask));
coef = polyfit(xLog, yLog, 1);
nf.chi = -coef(1);                  % PSD ~ f^coef(1); chi is the positive form

yHat  = polyval(coef, xLog);
ssRes = sum((yLog - yHat) .^ 2);
ssTot = sum((yLog - mean(yLog)) .^ 2);
nf.r2 = 1 - ssRes / max(ssTot, eps);

%% ========================================================================
%  MATCHED SURROGATE + THRESHOLD SWEEP
%  ========================================================================

nSurr = round(durSurr * fs);
surr  = gen1of(nSurr, nf.chi, seed);

% identical detection front-end; global z (whole surrogate as one NREM bout)
sigSurr = ripp_sigPrep(surr, fs, 'detectMet', method.detectMet, ...
    'passband', method.passband, 'zMet', 'nrem', 'nremTimes', [0 durSurr]);

% sweep the peak threshold; surrogate rate falls monotonically with it
thrGrid = 2 : 0.25 : 8;
rate    = nan(size(thrGrid));
for iThr = 1:numel(thrGrid)
    thr = method.thr;
    thr(2) = thrGrid(iThr);
    thr(1) = max(0.5, thrGrid(iThr) - thrGap);
    rSurr = ripp_times(sigSurr, fs, 'thr', thr, 'limDur', method.limDur);
    rate(iThr) = size(rSurr.times, 1) / durSurr;
end
nf.rate    = rate;
nf.thrGrid = thrGrid;

% first grid point at or below target; keep the fixed threshold if none reach it
iHit = find(rate <= targetFP, 1);
if ~isempty(iHit)
    thrCal      = thrGrid(iHit);
    nf.thrCal   = thrCal;
    nf.rateCal  = rate(iHit);
end

end     % EOF


% =========================================================================
%  LOCALS
% =========================================================================
function y = gen1of(n, chi, seed)
% Gaussian noise whose power spectrum follows 1/f^chi (no oscillation), via
% spectral shaping of white noise. Returned with zero mean and unit variance;
% amplitude is irrelevant since the detector z-scores.
if ~isempty(seed)
    rng(seed);
end
x = randn(n, 1);
X = fft(x);
nHalf = floor(n / 2);
f = (1:nHalf)';
amp = f .^ (-chi / 2);              % power ~ f^-chi -> amplitude ~ f^(-chi/2)

H = ones(n, 1);
H(2:nHalf + 1) = amp;
H(n:-1:n - nHalf + 2) = amp(1:nHalf - 1);   % mirror so ifft is real

y = real(ifft(X .* H));
y = y - mean(y);
sd = std(y);
if sd > 0
    y = y / sd;
end
end

% -------------------------------------------------------------------------
function idx = bouts2idx(bouts, fs, nMax)
% linear sample indices covered by [start end] bouts (s), clipped to nMax
idx = [];
smp = round(bouts * fs) + 1;
for iBout = 1:size(smp, 1)
    i0 = max(1, smp(iBout, 1));
    i1 = min(nMax, smp(iBout, 2));
    if i0 <= i1
        idx = [idx, i0:i1];  %#ok<AGROW>
    end
end
idx = idx(:);
end
