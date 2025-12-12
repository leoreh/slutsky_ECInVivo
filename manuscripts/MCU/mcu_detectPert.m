function [pertIdx, tAxis, debugData] = mcu_detectPert(frMat, varargin)
% MCU_DETECTPERT Detects perturbation onset using a robust Hybrid Strategy.
%
% This function implements a robust 3-stage strategy to detect network
% perturbations (e.g., Baclofen onset) from population firing rates avoiding
% false positives from transient fluctuations.
%
% Strategy:
%   1. Identify the top steep drops in the population median firing rate
%      and select the one leading to the lowest state.
%   2. Refine the onset time using a Geometric Knee Detection method.
%   3. Validate that the drop is a global network event by checking the
%      participation rate of individual units.
%
% INPUTS:
%   frMat        (double) Firing rate matrix [nUnits x nSamples].
%
% OPTIONAL (Key-Value Pairs):
%   binSize      (double) Bin size of each FR sample in seconds. {60}
%   flgPlot      (logical) Whether to plot the detection results. {false}
%
% OUTPUTS:
%   pertIdx      (double) Index of the perturbation onset (1-based).
%   tAxis        (double) Time axis vector (in hours), zero-centered at
%                the perturbation onset.
%   debugData    (struct) Intermediate data for debugging/plotting.
%
% See also: MCU_FRTBL, FR_DENOISE, MEA_FRR

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addRequired(p, 'frMat', @isnumeric);
addParameter(p, 'binSize', 60, @isnumeric); % Default 60 seconds
addParameter(p, 'flgPlot', false, @islogical);

parse(p, frMat, varargin{:});
binSize = p.Results.binSize;
flgPlot = p.Results.flgPlot;

% Constants
srchStart = 20;
srchEnd = 96;
winDropCheck = 30;      % Minutes to check for lowest point after drop
bslDur = 4;             % Hours to calculate baseline prior to searchStart
nCandidates = 10;       % Number of steep drops to evaluate

% Convert time constants to bins
srchStart = max(1, round(srchStart * 3600 / binSize));
winDrop   = round(winDropCheck * 60 / binSize);
winBsl    = round(bslDur * 3600 / binSize);

%% ========================================================================
%  COARSE IDENTIFICATION
%  ========================================================================

% Calculate Robust Population State. Use median to be resistant to
% single-unit outliers.
mfr = median(frMat, 1, 'omitnan');

% Denoise the trace
mfrSm = fr_denoise(mfr, 'flgPlot', false, 'frameLen', 60);

% Calculate Derivative and smooth
smoothWin = round(300 / binSize); % 5 min smoothing
mfrD = movmean(diff(mfrSm), smoothWin, 'Endpoints', 'fill');

% Find Top N Steepest Drops
% Restrict search to valid window
if isinf(srchEnd)
    srchEnd = length(mfrD) - winDrop;
else
    srchEnd = min(length(mfrD), round(srchEnd * 3600 / binSize)) - winDrop;
end

% Ensure indices are valid
srchStart = min(srchStart, length(mfrD)-1);
srchEnd = max(srchStart, min(srchEnd, length(mfrD)-1));
srchWin = srchStart:srchEnd;

% Find valleys (most negative peaks) in derivative
% We invert dPop to use findpeaks for local minima
[~, locs] = findpeaks(-mfrD(srchWin), 'SortStr', 'descend', 'NPeaks', nCandidates);
candidateIdx = srchWin(locs); % Map back to original indices

if isempty(candidateIdx)
    error('mcu:detectPert', 'No drops detected');
end

% Select the candidate that results in the *lowest absolute firing rate*
% immediately following the drop. This distinguishes major transitions
% from transient fluctuations.
minPostVal = inf;
anchorIdx = candidateIdx(1);

for iLoc = 1:length(candidateIdx)
    idx = candidateIdx(iLoc);

    % Check window after drop
    winEnd = min(length(mfrSm), idx + winDrop);
    postWin = idx : winEnd;
    minValInWin = min(mfrSm(postWin));

    if minValInWin < minPostVal
        minPostVal = minValInWin;
        anchorIdx = idx;
    end
end

%% ========================================================================
%  FINE REFINEMENT
%  ========================================================================
%  Geometric Knee Detection: Intersection of tangent and baseline.

% Define Pre-Event Baseline
% Median of the 4 hours prior to the *search window* (global baseline)
blEnd   = srchStart;
blStart = max(1, blEnd - winBsl);

% Fallback if recording is short
if blEnd <= 10 % Just arbitrary small number
    blEnd = max(1, anchorIdx - round(600/binSize));
    blStart = max(1, blEnd - winBsl);
end

baselineVal = median(mfr(blStart:blEnd), 'omitnan');

% Get Tangent at Anchor
% Slope at anchorIdx from the smoothed derivative
slope = mfrD(anchorIdx);
yAnchor = mfrSm(anchorIdx);
xAnchor = anchorIdx;

% Project Backward to Baseline
% Line Eq: y - y1 = m(x - x1)
% Intersection where y = baselineVal:
% xPert = xAnchor + (baselineVal - yAnchor) / slope

if abs(slope) < 1e-6
    % Slope too flat, fallback to anchor
    pertIdx = anchorIdx;
else
    dx = (baselineVal - yAnchor) / slope;
    pertIdxRaw = xAnchor + dx;
    pertIdx = round(pertIdxRaw);
end

% Bound check
pertIdx = max(1, min(length(mfr), pertIdx));

%% ========================================================================
%  CONSISTENCY CHECK
%  ========================================================================

% Check % of units decreasing
% Pre: 1 hour before pert
% Post: 1 hour after pert
winChk = round(3600 / binSize);
preIdx = max(1, pertIdx - winChk) : max(1, pertIdx);
postIdx = min(pertIdx + 1, size(frMat, 2)) : min(pertIdx + winChk + 1, size(frMat, 2));

muPre = mean(frMat(:, preIdx), 2, 'omitnan');
muPost = mean(frMat(:, postIdx), 2, 'omitnan');

% Count Decreasers
isDecr = muPost < muPre;
% Avoid division by zero
validUnits = ~isnan(muPre + muPost);
if sum(validUnits) > 0
    pctDecr = sum(isDecr & validUnits) / sum(validUnits) * 100;
else
    pctDecr = 0;
end

if pctDecr < 50
    warning('mcu:detectPert:Consistency', ...
        'Low consistency drop: only %.1f%% of units decreased FR.', pctDecr);
end

%% ========================================================================
%  CONSTRUCT OUTPUT
%  ========================================================================

% Time axis
nBins = size(frMat, 2);
dt_min = binSize / 60;
tAxis = ((1:nBins) - pertIdx) * dt_min / 60;

% Detection data
debugData.popMed = mfrSm;
debugData.candidates = candidateIdx;
debugData.anchorIdx = anchorIdx;
debugData.baseline = baselineVal;
debugData.consistency = pctDecr;
debugData.tangentSlope = slope;
debugData.tangentPoint = [xAnchor, yAnchor];

%% ========================================================================
%  PLOT
%  ========================================================================

if flgPlot

    hFig = figure('Name', 'Perturbation Detection', ...
        'Color', 'w', 'Position', [100 100 1000 600]);
    hold on;

    % Population Median
    plot(tAxis, mfrSm, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Pop. Median');

    % Baseline
    yline(baselineVal, 'b--', 'LineWidth', 1, 'DisplayName', 'Baseline');

    % Candidates
    plot(tAxis(candidateIdx), mfrSm(candidateIdx), 'ro', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 4, 'DisplayName', 'Candidates');

    % Anchor
    plot(tAxis(anchorIdx), mfrSm(anchorIdx), 'ks', ...
        'MarkerFaceColor', 'y', 'MarkerSize', 8, 'DisplayName', 'Anchor');

    % Tangent Line
    % Determine range for line plotting (e.g. +/- 30 mins)
    winTan = round(30 * 60 / binSize);
    xRng = [pertIdx - winTan, anchorIdx + winTan];
    yC = yAnchor + slope * (xRng - xAnchor);
    plot((xRng - pertIdx) * dt_min, yC, 'g-', 'LineWidth', 1, 'DisplayName', 'Tangent');

    % Detected Onset
    xline(tAxis(pertIdx), 'r--', 'LineWidth', 2, 'DisplayName', 'Detected Onset');

    title(sprintf('Perturbation Detection (Consistency=%.1f%%)', pctDecr));
    xlabel('Time (min)');
    ylabel('Median FR (Hz)');

    % Derivative on right axis
    yyaxis right
    plot(tAxis(1:length(mfrD)), mfrD, '-', 'Color', [0.5 0.5 0.5, 0.3], ...
        'LineWidth', 1, 'DisplayName', 'Derivative');
    ylabel('Derivative');
    hAx = gca;
    hAx.YAxis(2).Color = [0.5 0.5 0.5];

    legend('Location', 'best');
    grid on;
end

end     % EOF

