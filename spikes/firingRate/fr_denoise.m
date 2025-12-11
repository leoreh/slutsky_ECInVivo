function fr = fr_denoise(frOrig, t, varargin)
% FR_DENOISE Applies Savitzky-Golay filter to smooth firing rate data.
%
% SUMMARY:
% This function applies a Savitzky-Golay filter to smooth firing rate data
% while dealing with gaps (NaNs).
%
% The smoothing process follows these steps:
%   1. Small gaps (< 5 samples) are interpolated using linear filling.
%   2. The filter is applied to each resulting non-NaN segment separately.
%   3. Original NaNs are restored in the final output (masking).
%
% INPUT (Required):
%   frOrig       - Matrix of raw firing rate values [Hz]. Units are rows.
%   t            - Time vector corresponding to firing rate data [s].
%
% INPUT (Optional Key-Value Pairs):
%   flgPlot      - Logical flag to generate smoothing visualization {false}.
%   polyOrder    - Polynomial order for Savitzky-Golay filter {3}.
%   frameLenSec  - Frame length for Savitzky-Golay filter [s] {600}.
%
% OUTPUT:
%   fr           - Matrix of smoothed firing rate values [Hz]. Same size as frOrig.
%
% DEPENDENCIES:
%   Signal Processing Toolbox (for sgolayfilt)
%
% HISTORY:
%   Sep 2024 - Extracted from mea_frRecovery.m as standalone function.
%   Dec 2024 - Renamed to fr_denoise and added gap handling logic.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING & INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addRequired(p, 'frOrig', @isnumeric);
addRequired(p, 't', @isnumeric);
addParameter(p, 'flgPlot', false, @islogical);
addParameter(p, 'polyOrder', 3, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'frameLenSec', 600, @(x) isnumeric(x) && isscalar(x) && x > 0);

parse(p, frOrig, t, varargin{:});
flgPlot = p.Results.flgPlot;
polyOrder = p.Results.polyOrder;
frameLenSec = p.Results.frameLenSec;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILTER PARAMETER SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the filter frame length in number of bins (samples)
dt = mean(diff(t)); % Average sampling interval
frameLen = round(frameLenSec / dt);

% Savitzky-Golay filter requires an odd frame length.
if mod(frameLen, 2) == 0
    frameLen = frameLen + 1;
end

% The polynomial order must be less than the frame length.
if polyOrder >= frameLen
    polyOrder = frameLen - 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% APPLY FILTERING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize
fr = nan(size(frOrig)); % Start with NaNs
nUnits = size(frOrig, 1);

for iUnit = 1:nUnits
    
    frUnit = frOrig(iUnit, :);
    
    % Grab nan indices
    nanIdx = isnan(frUnit);

    % Fill small gaps (< 5 samples)
    frUnit = fillmissing(frUnit, 'linear', 'MaxGap', 4);

    % Identify non-NaN segments
    bouts = binary2bouts('vec', ~nanIdx);
    bouts(:, 2) = bouts(:, 2) - 1;
    
    for iBout = 1 : size(bouts, 1)
        boutIdx = bouts(iBout, 1) : bouts(iBout, 2);
        if length(boutIdx) >= frameLen
            frUnit(boutIdx) = sgolayfilt(frUnit(boutIdx), polyOrder, frameLen);
        end
    end
    
    % Enforce non-negativity, restore nan, and fill
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

    for iSmpl = 1 : nSmpl
        idx = smplIndices(iSmpl);
        % Plot original trace (semi-transparent)
        plot(t/60, frOrig(idx, :), 'Color', [colors(iSmpl,:), 0.4], 'LineWidth', 1);
        % Plot smoothed trace
        hLgd(iSmpl) = plot(t/60, fr(idx, :), 'Color', colors(iSmpl,:), 'LineWidth', 2);
        txtLgd{iSmpl} = sprintf('Unit %d', idx);
    end

    xlabel('Time (min)');
    ylabel('Firing Rate (Hz)');
    title(sprintf('Example Smoothed Units (n=%d)\\nSav-Gol Filter: Order %d, Frame %.1f min', ...
        nSmpl, polyOrder, frameLenSec/60));
    legend(hLgd, txtLgd, 'Location', 'eastoutside');
    grid on;
    box on;
    xlim([t(1)/60, t(end)/60]);
end

end
