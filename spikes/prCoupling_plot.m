function prCoupling_plot(prc, varargin)
% PRCOUPLING_PLOT Generates summary plots for population coupling analysis.
%
% SUMMARY:
% This function takes the output structure from prCoupling.m and creates a
% figure summarizing key aspects of population coupling analysis. It visualizes:
%   - Distribution of z-scored population coupling (prc0_norm) across units
%   - Raw STPR curves for all units, sorted by coupling strength
%   - Comparison of actual vs shuffled STPR values at zero lag
%   - Relationship between raw STPR and z-scored coupling
%   - Population-averaged STPR and shuffled STPR curves
%   - Example unit STPR curves with shuffled distributions
%   - Summary statistics and analysis parameters
%
% INPUT:
%   prc             Structure containing population coupling results from
%                   prCoupling.m. Must include fields:
%                     .prc0_norm: Z-scored population coupling
%                     .prc0: Raw STPR at zero lag
%                     .stpr: Full STPR curves
%                     .stpr_shfl: Shuffled STPR curves
%                     .prc0_shfl: Shuffled STPR values at zero lag
%                     .t: Time vector for STPR curves
%                     .info: Analysis parameters
%   basepath        (Optional) Path to recording session directory {pwd}
%   flgSaveFig      (Optional) Logical flag to save the figure {true}
%
% OUTPUT:
%   None. Generates and optionally saves a figure.
%
% DEPENDENCIES:
%   setMatlabGraphics (custom)
%   plot_stdshade (custom)
%   plot_scatterCorr (custom)
%
% HISTORY:
%   Aug 2024 LH - Created and updated

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARGUMENT PARSING & INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define input parameters and their defaults
p = inputParser;
addRequired(p, 'prc', @isstruct);
addParameter(p, 'basepath', pwd, @ischar);
addParameter(p, 'flgSaveFig', true, @islogical);

% Parse input arguments
parse(p, prc, varargin{:});
prc = p.Results.prc;
basepath = p.Results.basepath;
flgSaveFig = p.Results.flgSaveFig;

% Set basepath and get basename
cd(basepath);
[~, basename] = fileparts(basepath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARATIONS & PARAMETER EXTRACTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAPHICS GENERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up figure with consistent graphics settings
setMatlabGraphics(true);
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
for i = 1:nExUnits
    unit = exUnits(i);
    plot(tStpr, prc.stpr(unit,:), 'Color', colors(i,:), 'LineWidth', 1);
    plot(tStpr, prc.stpr_shfl(unit,:), '--', 'Color', colors(i,:), 'LineWidth', 0.5);
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

% --- Shuffled distribution for example units ---
subplot(3, 3, 6);
hold on;
for i = 1:nExUnits
    unit = exUnits(i);
    histogram(prc.prc0_shfl(unit,:), 20, 'Normalization', 'probability',...
        'FaceColor', colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    xline(prc.prc0(unit), '--', 'Color', colors(i,:), 'LineWidth', 1);
end
xlabel('STPR at Zero Lag');
ylabel('Probability');
title('Shuffled Distributions');
box off;

% --- Summary statistics and parameters ---
subplot(3, 3, 9);
hold on;
% Create text box with summary statistics
stats = {
    sprintf('Total Units: %d', nUnits),
    sprintf('Mean Coupling: %.2f', mean(prc.prc0_norm, 'omitnan')),
    sprintf('Median Coupling: %.2f', median(prc.prc0_norm, 'omitnan')),
    sprintf('STPR Window: %.1f s', winStpr),
    sprintf('Shuffles: %d', nShuffles),
    sprintf('Bin Size: %.3f s', prc.info.binSize),
    sprintf('Gaussian HW: %.3f s', prc.info.gkHw)
};
text(0.1, 0.5, stats, 'Units', 'normalized', 'FontName', 'Consolas');
axis off;

% Adjust figure layout
set(gcf, 'Position', [100, 100, 1200, 1000]);

% Save figure if requested
if flgSaveFig
    figFile = fullfile(basepath, [basename, '.prc.fig']);
    savefig(fh, figFile);
    fprintf('Saved figure to %s\n', figFile);
end

end 