function fr = fr_denoise(frOrig, varargin)
% FR_DENOISE Smooths firing rate data using normalized Gaussian convolution.
%
% SUMMARY:
% This function smooths firing rate data using normalized Gaussian
% convolution, which naturally handles NaN gaps and segment edges.
%
% NaN positions contribute zero weight in both the numerator (signal)
% and denominator (normalization), so they are automatically excluded
% from the local weighted average. This eliminates edge artifacts by
% construction — a weighted average cannot produce values outside the
% range of its inputs (no polynomial overshoot).
%
% INPUT (Required):
%   frOrig       - Matrix of raw firing rate values. Units are rows.
%
% INPUT (Optional Key-Value Pairs):
%   flgPlot      - Logical flag to generate smoothing visualization {false}.
%   frameLen     - Frame length for Gaussian kernel in samples {60}.
%
% OUTPUT:
%   fr           - Matrix of smoothed firing rate values [Hz]. Same size as frOrig.
%
% DEPENDENCIES:
%   Signal Processing Toolbox (for gausswin)
%
% HISTORY:
%   Sep 2024 - Extracted from mea_frRecovery.m as standalone function.
%   Dec 2024 - Renamed to fr_denoise and added gap handling logic.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING & INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addRequired(p, 'frOrig', @isnumeric);
addParameter(p, 'flgPlot', false, @islogical);
addParameter(p, 'frameLen', 60, @(x) isnumeric(x) && isscalar(x) && x > 0);

parse(p, frOrig, varargin{:});
flgPlot = p.Results.flgPlot;
frameLen = p.Results.frameLen;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KERNEL SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ensure odd kernel length
if mod(frameLen, 2) == 0
    frameLen = frameLen + 1;
end

% gaussian kernel
gk = gausswin(frameLen)';
gk = gk / sum(gk);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% APPLY FILTERING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize
fr = nan(size(frOrig));
nUnits = size(frOrig, 1);

for iUnit = 1:nUnits

    frUnit = frOrig(iUnit, :);

    % grab nan indices
    nanIdx = isnan(frUnit);

    if all(nanIdx)
        continue
    end

    % fill small gaps (< 5 samples)
    frUnit = fillmissing(frUnit, 'linear', 'MaxGap', 4);

    % normalized gaussian convolution. nan positions are zeroed in both
    % signal and weight vectors, so they are excluded from the local
    % weighted average. division by the convolved weights corrects for
    % partial kernel overlap at edges and near gaps.
    frZero = frUnit;
    frZero(isnan(frZero)) = 0;
    wt = double(~isnan(frUnit));
    frUnit = conv(frZero, gk, 'same') ./ conv(wt, gk, 'same');

    % enforce non-negativity and restore original nan positions
    frUnit(frUnit < 0) = 0;
    frUnit(nanIdx) = NaN;
    fr(iUnit, :) = frUnit;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flgPlot
    figure('Name', 'Firing Rate Smoothing (Savitzky-Golay)', ...
        'NumberTitle', 'off', 'Position', [100 100 1200 800]);

    % --- Plot Example Units ---
    ax = subplot(1,1,1);
    hold on;
    nUnits = size(frOrig, 1);
    nSmpl = min(5, nUnits);
    rng(1); % for reproducibility
    smplIndices = randperm(nUnits, nSmpl);
    colors = lines(nSmpl);

    hLgd = gobjects(nSmpl, 1);
    txtLgd = cell(nSmpl, 1);

    t = 1:size(frOrig, 2);

    for iSmpl = 1 : nSmpl
        idx = smplIndices(iSmpl);
        % Plot original trace (semi-transparent)
        plot(t, frOrig(idx, :), 'Color', [colors(iSmpl,:), 0.4], 'LineWidth', 1);
        % Plot smoothed trace
        hLgd(iSmpl) = plot(t, fr(idx, :), 'Color', colors(iSmpl,:), 'LineWidth', 2);
        txtLgd{iSmpl} = sprintf('Unit %d', idx);
    end

    xlabel('Time (Samples)');
    ylabel('Firing Rate');
    title(sprintf('Example Smoothed Units (n=%d)\\nGaussian Filter: Frame %d Smpls', ...
        nSmpl, frameLen));
    legend(hLgd, txtLgd, 'Location', 'eastoutside');
    grid on;
    box on;
    xlim([t(1), t(end)]);
end

end
