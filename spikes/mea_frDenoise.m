function fr = mea_frDenoise(frOrig, t, varargin)
% MEAFRDENOISE Applies Savitzky-Golay filter to smooth firing rate data.
%
% SUMMARY:
% This function applies a Savitzky-Golay filter to smooth firing rate data
% while preserving important temporal features. This method is chosen for
% its ability to smooth data while preserving the shape of important
% features like sharp onsets/offsets, peaks, and non-monotonic trends
% (e.g., recovery with overshoot).
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
fr = zeros(size(frOrig));
nUnits = size(frOrig, 1);

for iUnit = 1:nUnits
    % Apply filter
    frClean = sgolayfilt(frOrig(iUnit, :), polyOrder, frameLen);

    % Enforce non-negativity, as polynomial fits can sometimes dip below zero
    frClean(frClean < 0) = 0;

    fr(iUnit, :) = frClean;
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