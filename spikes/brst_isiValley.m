function [isiVal, stats] = brst_isiValley(spktimes, varargin)
% BRST_ISIVALLEY Detects the ISI valley separating bursts from non-bursts.
%
%   isiValley = BRST_ISIVALLEY(SPKTIMES, ...) calculates the inter-spike
%   intervals (ISIs) for all cells, aggregates them, and finds the valley
%   in the log-log histogram (or KDE) that typically separates burst ISIs
%   from tonic ISIs.
%
%   INPUTS:
%       spktimes    - (cell) Spike times per unit (e.g., {unit1, unit2}).
%                     Times should be in seconds.
%       varargin    - (param/value) Optional parameters:
%                     'rngVal'       : (vec) [min, max] ISI (s) to search for valley {[0.005, 0.25]}
%                     'bandwidth'    : (num) Bandwidth for KDE (log10 scale) {0.05}
%                     'flgPlot'      : (log) Plot the PDF and detected valley {true}
%
%   OUTPUTS:
%       isiVal      - (num) The ISI value (in seconds) at the valley.
%       stats       - (struct) Additional stats (PDF, grid, peaks).
%
%   EXAMPLE:
%       val = brst_isiValley(spktimes, 'flgPlot', true);
%

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

p = inputParser;
addParameter(p, 'rngVal', [0.005, 0.25], @isnumeric);
addParameter(p, 'bandwidth', 0.05, @isnumeric);
addParameter(p, 'flgPlot', true, @islogical);

parse(p, varargin{:});
rngVal       = p.Results.rngVal;
bw           = p.Results.bandwidth;
flgPlot      = p.Results.flgPlot;


%% ========================================================================
%  PREPARATIONS
%  ========================================================================

% Flatten spktimes if needed and handle single vector input
if isnumeric(spktimes)
    spktimes = {spktimes};
end

% Collect all ISIs
isiVec = [];
for iUnit = 1:length(spktimes)
    st = spktimes{iUnit};
    if isempty(st) || length(st) < 2
        continue;
    end
    st = st(:); % Ensure column
    isiVec = [isiVec; diff(st)]; %#ok<AGROW>
end

% Remove invalid ISIs
isiVec(isiVec <= 0) = [];

if isempty(isiVec)
    warning('brst_isiValley:NoISIs', 'No valid ISIs found.');
    isiVal = nan;
    stats = struct();
    return;
end


%% ========================================================================
%  CALCULATION
%  ========================================================================

% Log transform
isiLog = log10(isiVec);

% Kernel Density Estimation
% Define evaluation grid: from min to max, with sufficient resolution
gridMin = min(isiLog);
gridMax = max(isiLog);
nPoints = 1000;
xi = linspace(gridMin, gridMax, nPoints);

% Use ksdensity for robust smoothing
[isiPdf, xi] = ksdensity(isiLog, xi, 'Bandwidth', bw);

% Find local minima (valleys) and maxima (peaks)
% We effectively invert the PDF to use findpeaks for valleys
[~, locsMin] = findpeaks(-isiPdf);
[~, locsPeaks] = findpeaks(isiPdf);

% Filter for valleys within the search range
% Convert search range to log10
rngLog = log10(rngVal);
valid_min_indices = locsMin(xi(locsMin) >= rngLog(1) & xi(locsMin) <= rngLog(2));

isiVal = nan;

if isempty(valid_min_indices)
    % Fallback: Find the absolute minimum in the search range
    % (Useful if the distribution is bimodal but has no clear local minimum
    % peak-to-peak, or is just a slope in that region)
    idx_in_range = xi >= rngLog(1) & xi <= rngLog(2);
    if any(idx_in_range)
        [~, min_idx_local] = min(isiPdf(idx_in_range));
        % Map back to global index
        idxs = find(idx_in_range);
        best_idx = idxs(min_idx_local);
        isiVal = 10^xi(best_idx);
    end
else
    % If multiple valleys, pick the one with the lowest density
    [~, best_min_idx_in_locs] = min(isiPdf(valid_min_indices));
    best_idx = valid_min_indices(best_min_idx_in_locs);
    isiVal = 10^xi(best_idx);
end

% Prepare stats output
stats.pdf = isiPdf;
stats.xi  = xi;
stats.peaks_log = xi(locsPeaks);
stats.all_isis_log = isiLog;


%% ========================================================================
%  PLOTTING
%  ========================================================================

if flgPlot
    [hFig, hAx] = plot_axSize('szOnly', false, 'flgFullscreen', true, ...
        'flgPos', true);

    % Plot Histogram
    hHist = histogram(isiLog, 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    hold on;

    % Plot KDE
    hPlt = plot(xi, isiPdf, 'k-', 'LineWidth', 2);

    % Plot Valley
    if ~isnan(isiVal)
        xline(log10(isiVal), 'r--', 'LineWidth', 1.5, 'Label', sprintf('Valley: %.1f ms', isiVal*1000));
        plot(log10(isiVal), interp1(xi, isiPdf, log10(isiVal)), 'ro', 'MarkerFaceColor', 'r');
    end

    % Formatting
    xlabel('Log_{10}(ISI) [sec]');
    ylabel('PDF');
    title('ISI Distribution & Burst Valley');
    grid on;

    % Add secondary x-axis labels for time
    xticks_log = ceil(min(xi)):floor(max(xi));
    xticks(xticks_log);
    xticklabels(arrayfun(@(x) sprintf('10^{%d}', x), xticks_log, 'UniformOutput', false));
end

end
