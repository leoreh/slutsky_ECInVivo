function frFit = mea_frFit(frUnit, recovOnset, pertOnset, varargin)
% MEAFRFIT Fits double-exponential recovery model to firing rate data.
%
% SUMMARY:
% This function fits a double-exponential model to capture firing rate
% recovery dynamics after a perturbation. The model captures three phases:
%   1. Baseline period: Robust linear fit to pre-perturbation data
%   2. Perturbation period: Linear decline from baseline to trough
%   3. Recovery period: Double-exponential model with primary recovery
%      and optional overshoot components
%
% INPUT (Required):
%   frUnit      - Vector of firing rate values [Hz] for a single unit.
%   recovOnset  - Index of recovery onset (trough of population activity).
%   pertOnset   - Index of perturbation onset.
%
% INPUT (Optional):
%   opts        - Optimization options for lsqcurvefit {default options}.
%
% OUTPUT:
%   frFit       - Structure containing fit results:
%     .troughFr       - Firing rate at recovery onset (trough) [Hz].
%     .p0             - Initial parameter guesses [frSS, kRecov, aOver, kOver].
%     .pFit           - Fitted parameters [frSS, kRecov, aOver, kOver].
%     .rsquare        - R-squared goodness of fit.
%     .overshootTime  - Time to peak overshoot [s] (if applicable).
%     .overShootFr    - Peak overshoot firing rate [Hz] (if applicable).
%     .fitCurve       - Complete fitted curve [baseline; perturbation; recovery].
%     .bslFr          - Baseline firing rate from robust linear fit [Hz].
%     .bslFit         - Robust linear fit parameters [intercept, slope].
%
% DEPENDENCIES:
%   Optimization Toolbox (for lsqcurvefit)
%   Statistics and Machine Learning Toolbox (for robustfit)
%
% HISTORY:
%   Sep 2024 - Extracted from mea_frRecovery.m as standalone function.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING & INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addRequired(p, 'frUnit', @isnumeric);
addRequired(p, 'recovOnset', @(x) isnumeric(x) && isscalar(x) && x > 0);
addRequired(p, 'pertOnset', @(x) isnumeric(x) && isscalar(x) && x > 0);
addOptional(p, 'opts', [], @(x) isempty(x) || isstruct(x));

parse(p, frUnit, recovOnset, pertOnset, varargin{:});
opts = p.Results.opts;

% Set default optimization options if not provided
if isempty(opts)
    opts = optimoptions('lsqcurvefit', 'Display', 'off', 'TolFun', 1e-6,...
        'MaxIter', 1000, 'FunctionTolerance', 1e-6);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE OUTPUT STRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

frFit = struct('troughFr', NaN, 'p0', NaN(1, 4), 'pFit', NaN(1, 4), ...
    'rsquare', NaN, 'overshootTime', NaN, 'overShootFr', NaN, ...
    'fitCurve', NaN, 'bslFr', NaN, 'bslFit', NaN);

% Hard limit on number of bins with spikes
if sum(frUnit > 0) < 10
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL FITTING PREPARATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert indices to relative time
tIdx = (1:length(frUnit))' - recovOnset;

% Select data for fitting
frPost = frUnit(recovOnset:end);
tPost = tIdx(recovOnset:end);
frFit.troughFr = frUnit(recovOnset);

% Steady-state rate: Mean of last 10% of data
frSS = mean(frPost(round(0.9 * length(frPost)):end), 'omitnan');

% Recovery rate: Based on time to 50% recovery, or last bin if no clear recovery
halfRecovFr = frFit.troughFr + 0.5 * (frSS - frFit.troughFr);
halfRecovIdx = find(frPost >= halfRecovFr, 1, 'first');
if ~isempty(halfRecovIdx)
    % Clear recovery: use time to 50% recovery
    tHalf = tPost(halfRecovIdx);
    kRecovGuess = log(2) / tHalf;
else
    % No clear recovery: use rate from last bin. This will be slow if the
    % unit hasn't recovered.
    lastBinFr = frPost(end);
    kRecovGuess = (lastBinFr - frFit.troughFr) / tPost(end);
end

% Overshoot amplitude: Based on max post-pert FR
[maxFr, maxIdx] = max(frPost);
if maxFr > frSS
    aOver = (maxFr - frSS) / (tPost(maxIdx) * exp(-1));
else
    aOver = 0;
end

% Overshoot rate: Based on time to peak
if aOver > 1e-6
    kOver = 1 / tPost(maxIdx);
else
    kOver = 1/1800;  % Default: characteristic time of 30 min
end

% Store initial parameters
frFit.p0 = [frSS, kRecovGuess, aOver, kOver];

% Define parameter bounds
lb = [0,   1e-6, 0,    1e-6];
ub = [3*max(frUnit), 1, 3*max(frUnit), 1];

% Override overshoot and exclude it from the model
% frFit.p0(3) = 0;
% ub(3) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PERFORM MODEL FIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define model function
modelFunc = @(p, tData) p(1) - (p(1) - frFit.troughFr) .* exp(-p(2) * tData) + ...
    p(3) * tData .* exp(-p(4) * tData);

% Perform the fit
[pFit, resnorm, ~, exitflag] = lsqcurvefit(modelFunc, frFit.p0, tPost, frPost(:), lb, ub, opts);
if exitflag == 0
    return
end

frFit.pFit = pFit;

% Calculate R-squared
ssTot = sum((frPost - mean(frPost)).^2);  % Total variance in the data
ssRes = resnorm;  % Sum of squared residuals (model vs data)
frFit.rsquare = 1 - (ssRes / ssTot);

% Calculate peak overshoot
if pFit(3) > 1e-6 && pFit(4) > 1e-6
    tPeak = 1 / pFit(4);
    frPeak = modelFunc(pFit, tPeak);
    frFit.overshootTime = tPeak;
    frFit.overShootFr = frPeak;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE BASELINE WITH ROBUST LINEAR FIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define baseline window (pre-perturbation onset) with margin
mrgnBins = 10;
bslWin = 1:max(1, pertOnset - mrgnBins);

% Fit robust linear regression to baseline period
[bslFit, ~] = robustfit(tIdx(bslWin), frUnit(bslWin));

% Generate baseline curve using the linear fit
bslCurve = bslFit(1) + bslFit(2) * tIdx(1:(pertOnset-1));

% Use fitted line for calculating baseline FR
frFit.bslFr = mean(bslCurve);
frFit.bslFit = bslFit;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE COMPLETE FITTED CURVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Perturbation period: exponential decay
frTrough = max([frFit.troughFr, eps]);
tPert = tIdx(pertOnset:(recovOnset-1));
pertCurve = linspace(bslCurve(end), frTrough, length(tPert))';

% Recovery period: fitted model
recovCurve = modelFunc(pFit, tPost);

% Combine segments
frFit.fitCurve = [bslCurve; pertCurve; recovCurve];

end 