function frr = mea_frRecovery(spktimes, varargin)
% MEAFRRECOVERY Calculates firing rate recovery parameters for all neurons.
%
% SUMMARY:
% This function analyzes firing rate recovery after perturbation:
%   1. Calculates firing rates using times2rate
%   2. Detects perturbation time from population MFR
%   3. Calculates recovery metrics for each unit
%
% INPUT (Required):
%   spktimes      - Cell array. spktimes{i} contains spike times (s) for neuron i.
%
% INPUT (Optional Key-Value Pairs):
%   basepath      - Path to recording session directory {pwd}.
%   flgSave       - Logical flag to save results to .frr.mat file {true}.
%   flgPlot       - Logical flag to plot denoising results {false}.
%   winLim        - [start end] time window to analyze [s] {[0 Inf]}.
%   spkThr        - Minimum number of spikes required per unit {300}.
%   frThr         - Minimum baseline firing rate [Hz] {0.1}.
%   bslThr  - Tolerance threshold for recovery (fraction of baseline) {0.1}.
%   wavelet       - Wavelet type for denoising {'sym4'}.
%   level         - Decomposition level for denoising {3}.
%
% OUTPUT:
%   frr           - Structure containing recovery analysis results:
%     .recovTime   - Recovery time for each unit [s].
%     .recovSlope  - Recovery slope (Hz/min) for each unit.
%     .steadyState - Steady state firing rate for each unit [Hz].
%     .pertTime    - Detected perturbation time [s].
%     .bslFr  - Baseline firing rates [Hz].
%     .fr          - Full firing rate curves for each neuron.
%     .t           - Time vector for rate curves [s].
%     .info        - Analysis parameters and metadata.
%
% DEPENDENCIES:
%   times2rate.m
%
% HISTORY:
%   Aug 2024 LH - Initial version.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING & INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define input parameters and their defaults/validation functions
p = inputParser;
addRequired(p, 'spktimes', @iscell);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgSave', true, @islogical);
addParameter(p, 'flgPlot', false, @islogical);
addParameter(p, 'flgDenoise', true, @islogical);
addParameter(p, 'winLim', [0 Inf], @(x) isnumeric(x) && numel(x)==2);
addParameter(p, 'spkThr', 300, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'frThr', 0.0001, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'bslThr', 0.1, @(x) isnumeric(x) && isscalar(x) && x>0);

% Parse input arguments
parse(p, spktimes, varargin{:});
spktimes = p.Results.spktimes;
basepath = p.Results.basepath;
flgSave = p.Results.flgSave;
flgPlot = p.Results.flgPlot;
flgDenoise = p.Results.flgDenoise;
winLim = p.Results.winLim;
spkThr = p.Results.spkThr;
frThr = p.Results.frThr;
bslThr = p.Results.bslThr;

% Set basepath and get basename
cd(basepath);
[~, basename] = fileparts(basepath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME AXIS & INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create time axis for binning
if isinf(winLim(2))
    winLim(2) = max(cellfun(@(x) max(x), spktimes, 'uni', true));
end

% Window spikes to analysis window
spktimes = cellfun(@(x) x(x >= winLim(1) & x <= winLim(2)), spktimes, 'UniformOutput', false);
nSpks = cellfun(@(x) length(x), spktimes, 'UniformOutput', true);
uGood = find(nSpks > spkThr);
nUnits = length(spktimes);

% Initialize output arrays
frBsl = nan(nUnits, 1);
recovTime = nan(nUnits, 1);
recovSlope = nan(nUnits, 1);
frRecov = nan(nUnits, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIRING RATE CALCULATION & DENOISING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate firing rates for all units
[frOrig, ~, t] = times2rate(spktimes(uGood), 'binsize', 30, 'c2r', true);
frOrig = frOrig';

% First smooth with Gaussian kernel to handle sparse spikes
kw = 3;
gk = gausswin(kw) / sum(gausswin(kw));

% Smooth each unit
for iUnit = 1:size(frOrig, 2)
    frOrig(:, iUnit) = conv(frOrig(:, iUnit), gk, 'same');
end

% Calculate population mean firing rate after smoothing
mfrOrig = mean(frOrig, 2, 'omitnan');

% Apply denoising if requested
if flgDenoise
    [fr, mfr] = denoise_fr(frOrig, mfrOrig, t);
else
    fr = frOrig;
    mfr = mfrOrig;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PERTURBATION DETECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find peaks in negative MFR (which gives us minima in original MFR)
[peakProms, pertIdx] = findpeaks(-mfr, ...
    'MinPeakProminence', std(mfr)/2, ... % Minimum prominence relative to std
    'MinPeakDistance', 30, ...           % Increased minimum distance between peaks (bins)
    'MinPeakHeight', -mean(mfr));        % Must be below mean

% If multiple peaks found, take the most prominent one
if length(pertIdx) > 1
    [~, maxPromIdx] = max(peakProms);
    pertIdx = pertIdx(maxPromIdx);
end

% If no significant peak found, use the global minimum
if isempty(pertIdx)
    [~, pertIdx] = min(mfr);
end

% Get the perturbation time
pertTime = t(pertIdx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RECOVERY ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define margin bins for more robust calculations
mrgnBins = 5;  % Number of bins to add/remove for calculations

% Analyze recovery for each unit
for iGood = 1:length(uGood)
    iUnit = uGood(iGood);
    
    % Get unit's firing rate
    uFr = fr(:, iGood);
    
    % Calculate baseline (pre-perturbation) firing rate with margin
    baseIdx = t < t(pertIdx - mrgnBins);
    frBsl(iUnit) = mean(uFr(baseIdx), 'omitnan');
    
    % Skip if baseline firing rate is too low
    if frBsl(iUnit) < frThr
        continue
    end
    
    % Find recovery time and calculate metrics
    [recovTime(iUnit), recovIdx(iUnit), recovSlope(iUnit), frRecov(iUnit)] = ...
        fit_recovery(uFr, t, pertTime, frBsl(iUnit), bslThr, mrgnBins);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT RECOVERY RELATIONSHIPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get valid units (those with recovery)
validUnits = find(~isnan(recovTime) & ~isnan(recovSlope));

% Create figure for recovery relationships
figure('Name', 'Recovery Relationships', 'NumberTitle', 'off');

% Plot baseline vs recovery firing rate
subplot(1,2,1);
plot_scatterCorr(frBsl(validUnits), fr(validUnits), ...
    'xLbl', 'Baseline Firing Rate (Hz)', 'yLbl', 'Recovery Firing Rate (Hz)');

% Plot baseline vs recovery slope
subplot(1,2,2);
plot_scatterCorr(frBsl(validUnits), recovSlope(validUnits), ...
    'xLbl', 'Baseline Firing Rate (Hz)', 'yLbl', 'Recovery Slope (Hz/min)');

% Adjust figure
set(gcf, 'Position', [100 100 1200 500]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT RECOVERY SLOPES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get valid units (those with recovery)
nSmpl = 10;
validUnits = find(~isnan(recovTime) & ~isnan(recovSlope));
% Select random sample
rng(1); % For reproducibility
smplUnits = validUnits(randperm(length(validUnits), nSmpl));

% Create figure
figure('Name', 'Recovery Slopes', 'NumberTitle', 'off');

% Plot MFR for context
subplot(2,1,1);
plot(t/60, mfr, 'k', 'LineWidth', 1);
hold on;
xline(pertTime/60, 'r--', 'Perturbation');
xlabel('Time (min)');
ylabel('MFR (Hz)');
title('Population Mean Firing Rate');

% Plot recovery slopes
subplot(2,1,2);
colors = lines(length(smplUnits));
for iUnit = 1:length(smplUnits)
    uIdx = smplUnits(iUnit);
    % Plot firing rate
    plot(t/60, fr(:, uIdx), 'Color', [colors(iUnit,:), 0.3], 'LineWidth', 1);
    hold on;
    % Plot recovery slope
    tRecov = [t(pertIdx), t(recovIdx(uIdx))];
    ySlope = [fr(pertIdx, uIdx), fr(recovIdx(uIdx), uIdx)];
    plot(tRecov/60, ySlope, '--', 'Color', colors(iUnit,:), 'LineWidth', 2);
end
xline(pertTime/60, 'r--', 'Perturbation');
xlabel('Time (min)');
ylabel('Firing Rate (Hz)');
title(sprintf('Recovery Slopes (n=%d units)', length(smplUnits)));

% Add legend with recovery times
legStr = arrayfun(@(x) sprintf('Unit %d (%.1f min)', x, recovTime(x)/60), ...
    smplUnits, 'UniformOutput', false);
legend(legStr, 'Location', 'eastoutside');

% Adjust figure
set(gcf, 'Position', [100 100 1200 800]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORGANIZE OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Store results
frr.bslFr = frBsl;              % Baseline firing rates
frr.frRecov = frRecov;          % Steady state firing rate
frr.recovTime = recovTime;      % Recovery time for each unit
frr.recovSlope = recovSlope;    % Recovery slope (Hz/min)
frr.pertTime = pertTime;        % Detected perturbation time
frr.fr = fr;                    % Full firing rate curves
frr.t = t;                      % Time vector

% Create info struct with parameters and metadata
frr.info = struct(...
    'nUnits', nUnits, ...
    'winLim', winLim, ...
    'spkThr', spkThr, ...
    'frThr', frThr, ...
    'bslThr', bslThr);

% Save results if requested
if flgSave
    save(fullfile(basepath, [basename, '.frr.mat']), 'frr', '-v7.3');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION: DENOISE FIRING RATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fr, mfr] = denoise_fr(frOrig, mfrOrig, t)
% DENOISE_FR Denoises firing rate data using wavelet denoising.
%
% Hardcoded parameters for consistent denoising across all analyses.
%
% INPUT:
%   frOrig      - Original firing rate matrix [time x units]
%   mfrOrig     - Original mean firing rate vector
%   t           - Time vector
%
% OUTPUT:
%   fr          - Denoised firing rate matrix
%   mfr         - Denoised mean firing rate vector

% Hardcoded denoising parameters
wavelet = 'sym8';
level = 10;
thresholdRule = 'Hard';

% Apply wavelet denoising
fr = wdenoise(frOrig, level, 'Wavelet', wavelet, 'ThresholdRule', thresholdRule);
mfr = wdenoise(mfrOrig, level, 'Wavelet', wavelet, 'ThresholdRule', thresholdRule);

% Ensure non-negative values
fr = max(0, fr);
mfr = max(0, mfr);

% Always plot results
figure('Name', 'Firing Rate Denoising', 'NumberTitle', 'off');

% Plot MFR
subplot(2,1,1); hold on
h1 = plot(t/60, mfrOrig, 'k', 'LineWidth', 1);
h2 = plot(t/60, mfr, 'r', 'LineWidth', 1.5);
set(h1, 'Color', [0 0 0 0.3]); % Make original signal semi-transparent
xlabel('Time (min)');
ylabel('MFR (Hz)');
title('Population Mean Firing Rate');
legend([h1 h2], {'Smoothed', 'Denoised'}, 'Location', 'best');

% Plot example units
subplot(2,1,2); hold on
nSmpl = min(5, size(fr,2));
rng(1); % For reproducibility
smplUnits = randperm(size(fr,2), nSmpl);
colors = lines(nSmpl);

for i = 1:nSmpl
    uIdx = smplUnits(i);
    h1 = plot(t/60, frOrig(:, uIdx), 'Color', [colors(i,:) 0.3], 'LineWidth', 1);
    h2 = plot(t/60, fr(:, uIdx), '--', 'Color', colors(i,:), 'LineWidth', 1.5);
    set(h1, 'Color', [colors(i,:) 0.3]); % Make original signal semi-transparent
end
xlabel('Time (min)');
ylabel('Firing Rate (Hz)');
title(sprintf('Example Units (n=%d)', nSmpl));

% Create legend entries
legEntries = cell(1, nSmpl*2);
for i = 1:nSmpl
    legEntries{i*2-1} = sprintf('Unit %d (smoothed)', smplUnits(i));
    legEntries{i*2} = sprintf('Unit %d (denoised)', smplUnits(i));
end
legend(legEntries, 'Location', 'eastoutside');

set(gcf, 'Position', [100 100 1200 800]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION: FIT RECOVERY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [recovTime, recovIdx, recovSlope, frRecov] = fit_recovery(fr, t, pertTime, bslFr, bslThr, mrgnBins)
% FIT_RECOVERY Finds recovery time and calculates recovery metrics.
%
% INPUT:
%   fr           - Firing rate vector.
%   t            - Time vector.
%   pertTime     - Perturbation time.
%   bslFr        - Baseline firing rate.
%   bslThr       - Tolerance threshold for recovery (fraction of baseline).
%   mrgnBins     - Number of margin bins to add/remove.
%
% OUTPUT:
%   recovTime    - Time to recovery [s].
%   recovIdx     - Index of recovery time.
%   recovSlope   - Recovery slope [Hz/min].
%   frRecov      - Recovery firing rate [Hz].

% Find perturbation index
[~, pertIdx] = min(abs(t - pertTime));

% Get post-perturbation data with margin
postIdx = t >= t(pertIdx + mrgnBins);
recovFr = fr(postIdx);
recovT = t(postIdx);

% Skip if not enough data points
if sum(~isnan(recovFr)) < 10
    recovTime = NaN;
    recovIdx = NaN;
    recovSlope = NaN;
    frRecov = NaN;
    return
end

% Calculate tolerance window
tolLow = bslFr * (1 - bslThr);
tolHigh = bslFr * (1 + bslThr);

% Find first time firing rate returns to baseline tolerance
inTol = recovFr >= tolLow & recovFr <= tolHigh;

% Find first sustained return to baseline (at least 3 consecutive bins)
sustainedDur = 3;
recovFound = false;
for i = 1:(length(inTol) - sustainedDur + 1)
    if all(inTol(i:i + sustainedDur - 1))
        recovIdx = find(t == recovT(i));
        recovTime = recovT(i) - pertTime;
        recovFound = true;
        break
    end
end

% If no sustained recovery found, use end of recording
if ~recovFound
    recovTime = t(end);
    recovIdx = length(t);
end

% Calculate recovery slope
pertIdxMarg = pertIdx + mrgnBins;
recovIdxMarg = recovIdx - mrgnBins;

% Ensure we have enough points for slope calculation
if recovIdxMarg - pertIdxMarg > 2
    x = t(pertIdxMarg:recovIdxMarg) / 60;  % Convert to minutes
    y = fr(pertIdxMarg:recovIdxMarg);
    % Use robust linear fit if enough points
    [p, ~] = robustfit(x, y);
    recovSlope = p(2);  % Slope in Hz/min
else
    % Use simple two-point slope if not enough points
    recovSlope = (fr(recovIdxMarg) - fr(pertIdxMarg)) / ...
        ((t(recovIdxMarg) - t(pertIdxMarg)) / 60);
end

% Calculate recovery firing rate
if recovFound
    % Use mean after recovery point
    frRecov = mean(fr(recovIdx - mrgnBins : end), 'omitnan');
else
    % Use mean of last 20 bins if no recovery found
    lastBins = 20;
    frRecov = mean(fr(end - lastBins + 1 : end), 'omitnan');
end

end

