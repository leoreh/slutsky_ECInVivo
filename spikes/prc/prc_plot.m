function prc_plot(prc, varargin)
% PRC_PLOT Generates summary plots for population coupling analysis.
%
%   PRC_PLOT(PRC, ...) creates a summary figure visualizing the population
%   coupling analysis results, including STPR curves, shuffled distributions,
%   and coupling strength metrics.
%
%   INPUTS:
%       prc         - (struct) Output from prc_calc.m containing:
%                     .prc0_norm : Z-scored population coupling
%                     .prc0      : Raw STPR at zero lag
%                     .stpr      : Full STPR curves
%                     .stpr_shfl : Shuffled STPR curves
%                     .prc0_shfl : Shuffled STPR values at zero lag
%                     .t         : Time vector for STPR curves
%                     .info      : Analysis parameters
%       varargin    - (param/value) Optional parameters:
%                     'basepath' : (char) Path to recording session directory {pwd}
%                     'flgSave'  : (log) Save figure to disk {true}
%
%   OUTPUTS:
%       None. Generates and optionally saves a figure.
%
%   See also: PRC_CALC, PLOT_STDSHADE, PLOT_SCATTERCORR
%
%   HISTORY:
%   Aug 2024 LH - Created and updated

%% ========================================================================
%  ARGUMENTS
%  ========================================================================

% Define input parameters and their defaults
p = inputParser;
addRequired(p, 'prc', @isstruct);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgSave', true, @islogical);

% Parse input arguments
parse(p, prc, varargin{:});
prc = p.Results.prc;
basepath = p.Results.basepath;
flgSave = p.Results.flgSave;

% Set basepath and get basename
cd(basepath);
[~, basename] = fileparts(basepath);

%% ========================================================================
%  PREPARATIONS
%  ========================================================================

% Extract key parameters
nUnits = length(prc.prc0_norm);
winStpr = prc.info.winStpr;
nShuffles = prc.info.nShuffles;
lagBins = round((winStpr / 2) / prc.info.binSize);
tStpr = (-lagBins:lagBins) * prc.info.binSize;

% Sort units by coupling strength for visualization
[~, sortIdx] = sort(prc.prc0_norm, 'descend');

% Select example units for detailed visualization
nExUnits = min(5, nUnits);
exUnits = round(linspace(1, nUnits, nExUnits));

%% ========================================================================
%  GRAPHICS GENERATION
%  ========================================================================

% Set up figure with consistent graphics settings
fh = figure('Name', [basename '_populationCoupling'], 'NumberTitle', 'off');

% --- Raw STPR curves for all units ---
subplot(3, 3, [1, 2]);
imagesc(tStpr, 1:nUnits, prc.stpr(sortIdx,:));
colormap('jet');
colorbar;
xlabel('Time [s]');
ylabel('Unit (sorted by coupling)');
title('STPR Curves');
box off;

% --- Example unit STPR curves ---
subplot(3, 3, [4, 5]);
hold on;
colors = lines(nExUnits);
for iUnit = 1:nExUnits
    unit = exUnits(iUnit);
    plot(tStpr, prc.stpr(unit,:), 'Color', colors(iUnit,:), 'LineWidth', 1);
    plot(tStpr, prc.stpr_shfl(unit,:), '--', 'Color', colors(iUnit,:), 'LineWidth', 0.5);
end
xline(0, '--k');
xlabel('Time [s]');
ylabel('Population Rate [Hz]');
title('Example Unit STPR');
legend(arrayfun(@(x) sprintf('Unit %d', x), exUnits, 'UniformOutput', false),...
    'Location', 'best', 'Box', 'off');
box off;

% --- Population-averaged STPR curves ---
subplot(3, 3, [7, 8]);
hold on;
plot_stdShade(prc.stpr, tStpr, gca, [0 0 0], 0.3);  % Actual data
plot_stdShade(prc.stpr_shfl, tStpr, gca, [1 0 0], 0.3);  % Shuffled data
xline(0, '--k');
xlabel('Time [s]');
ylabel('Population Rate [Hz]');
title('Population-Averaged STPR');
legend('Actual', 'Shuffled', 'Location', 'best');
box off;

% --- Correlation between raw STPR and z-score ---
subplot(3, 3, 3);
plot_scatterCorr(prc.prc0, prc.prc0_norm, 'xLbl', 'Raw STPR at Zero Lag', ...
    'yLbl', 'Population Coupling (Z-score)', 'flgOtl', true);
title('STPR vs. Coupling Strength');

% --- Population Shuffle Overlay ---
subplot(3, 3, 6);
hold on;
% Histogram of Data (Actual Population Coupling)
h1 = histogram(prc.prc0, 30, 'Normalization', 'pdf', ...
    'FaceColor', [0.85 0.325 0.098], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
% Histogram of Null Distribution (All Shuffles)
h2 = histogram(prc.prc0_shfl(:), 30, 'Normalization', 'pdf', ...
    'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.4);

xlabel('STPR at Zero Lag (Hz)');
ylabel('Probability Density');
title('Population Shuffle Overlay');
legend([h1, h2], {'Data', 'Null'}, 'Location', 'best', 'Box', 'off');
box off;

% --- MCMC Chain Stability ---
subplot(3, 3, 9);
hold on;
% Calculate mean metric across units for each shuffle iteration
% prc0_shfl is (nUnits x nShuffles)
mean_shfl_trace = mean(prc.prc0_shfl, 1, 'omitnan');
mean_data_val = mean(prc.prc0, 'omitnan');

plot(mean_shfl_trace, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
yline(mean_data_val, '--r', 'Data Mean', 'LabelHorizontalAlignment', 'left', 'FontSize', 8);

xlabel('Shuffle Iteration');
ylabel('Mean STPR(0) [Hz]');
title('MCMC Stability');
grid on; box off;

% Add stats as annotation in subplot 3 (Scatter) or overlay on stability plot
% Let's overlay compact stats on the stability plot
statsText = {
    sprintf('N=%d', nUnits),
    sprintf('Iter=%d', nShuffles)
    };
text(0.95, 0.95, statsText, 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'BackgroundColor', 'none', 'FontSize', 8);

% Adjust figure layout
set(gcf, 'Position', [100, 100, 1200, 1000]);

% Save figure if requested
if flgSave
    figFile = fullfile(basepath, [basename, '.prc.fig']);
    savefig(fh, figFile);
    fprintf('Saved figure to %s\n', figFile);
end

end