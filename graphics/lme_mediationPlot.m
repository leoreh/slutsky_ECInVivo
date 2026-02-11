function lme_mediationPlot(res, varargin)
% LME_MEDIATIONPLOT Visualizes the results of lme_mediation.
%
%   LME_MEDIATIONPLOT(RES) creates a 2x2 tile layout visualizing the
%   mediation paths (A, B, C, C') using the data and statistics from the
%   results structure.
%
%   Plots:
%       1. Path C (Total): X vs Y.
%       2. Path A (Mediation): X vs M.
%       3. Path B (Partial): M vs Y_partial (Effect of M on Y | X).
%       4. Path C' (Direct): X vs Y_partial (Effect of X on Y | M).
%
%   INPUTS:
%       res - (struct) Output from lme_mediation.
%
%   OPTIONAL:
%       'Parent' - (handle) Target figure or layout.
%
%   See also: LME_MEDIATION, PLOT_SCAT

%% ========================================================================
%  INPUT
%  ========================================================================
p = inputParser;
addRequired(p, 'res', @isstruct);
addParameter(p, 'Parent', [], @(x) isempty(x) || isgraphics(x));
parse(p, res, varargin{:});

hPar = p.Results.Parent;
if isempty(hPar)
    figure('Color', 'w', 'Position', [100 100 1000 800]);
    hPar = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
else
    if isa(hPar, 'matlab.ui.Figure')
        figure(hPar);
        hPar = tiledlayout(hPar, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    end
end

% Check if plot data exists
if ~isfield(res, 'plot') || isempty(res.plot)
    error('Res structure does not contain plot data. Please re-run lme_mediation.');
end

% Determine Fit Types
fitType = 'linear'; 
if ~isnumeric(res.plot.X)
    fitType = 'None'; 
end

% Extract Variable Names from Paths table for labeling
rNames = res.paths.Properties.RowNames;
pathA_Row = find(contains(rNames, 'Path A'));
pathB_Row = find(contains(rNames, 'Path B'));
pathC_Row = find(contains(rNames, 'Path C ('));
pathCp_Row = find(contains(rNames, 'Path C'''));

% Helper to get stats string
getStats = @(rowIdx) sprintf('%c = %.3f, p = %.4f', ...
    946, res.paths.Estimate(rowIdx), res.paths.PValue(rowIdx)); % Beta symbol


%% ========================================================================
%  PATH A (MEDIATION: X -> M)
%  ========================================================================
nexttile(hPar);
plot_scat([], res.plot.X, res.plot.M, ...
    'fitType', fitType, 'flgStats', true, 'alpha', 0.5);
title({'Path A', getStats(pathA_Row)}, 'FontWeight', 'normal');
xlabel('X');
ylabel('M');

%% ========================================================================
%  PATH B (PARTIAL: M -> Y | X)
%  ========================================================================
% X-Axis: Mediator (M)
% Y-Axis: Partial Residuals of Y (slope should correspond to beta_B)
nexttile(hPar);
plot_scat([], res.plot.M, res.plot.Y_part_M, ...
    'fitType', fitType, 'flgStats', true, 'alpha', 0.5);
title({'Path B', getStats(pathB_Row)}, 'FontWeight', 'normal');
xlabel('M');
ylabel('Y | X');

%% ========================================================================
%  PATH C' (DIRECT: X -> Y | M)
%  ========================================================================
% X-Axis: Treatment (X)
% Y-Axis: Partial Residuals of Y (slope should correspond to beta_C')
nexttile(hPar);
plot_scat([], res.plot.X, res.plot.Y_part_X, ...
    'fitType', fitType, 'flgStats', true, 'alpha', 0.5);
title({'Path C'' (Direct)', getStats(pathCp_Row)}, 'FontWeight', 'normal');
xlabel('X');
ylabel('Y | M');

%% ========================================================================
%  PATH C (TOTAL EFFECT: X -> Y)
%  ========================================================================
nexttile(hPar);
plot_scat([], res.plot.X, res.plot.Y, ...
    'fitType', fitType, 'flgStats', true, 'alpha', 0.5);
title({'Path C (Total)', getStats(pathC_Row)}, 'FontWeight', 'normal');
xlabel('X');
ylabel('Y');

end
