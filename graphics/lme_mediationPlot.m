function lme_mediationPlot(res, mdlA, mdlBC, tbl, varargin)
% LME_MEDIATIONPLOT Visualizes mediation paths from lme_mediation results.
%
%   LME_MEDIATIONPLOT(RES, MDLA, MDLBC, TBL, ...) creates a 2x2 scatter
%   plot of the four mediation paths (A, B, C, C'). When RES is a struct
%   array (from stratified analysis), groups are overlaid on the same axes.
%
%   Plots:
%       1. Path A (X -> M)
%       2. Path B (M -> Y | X)         partial residuals
%       3. Path C' (Direct X -> Y | M) partial residuals
%       4. Path C (Total X -> Y)
%
%   INPUTS:
%       res     - (struct array) Output from lme_mediation.
%       mdlA    - (LME/GLME) Pre-fitted Model A (for fitted values).
%       mdlBC   - (LME/GLME) Pre-fitted Model BC (for residuals).
%       tbl     - (table) Data table (all groups, before any plot transforms).
%       varargin - (param/value) Optional parameters:
%                  'xVar'   : (char) X variable name. Default: inferred
%                             from res, but can override for transformed
%                             columns (e.g., 'pBspk_trans').
%                  'mVar'   : (char) M variable name. Default: mdlA
%                             ResponseName, but can override.
%                  'grpVar' : (char) Grouping variable for colors {''}.
%                  'Parent' : (handle) Target figure or layout {[]}.
%
%   EXAMPLE:
%       res = lme_mediation(mdlA, mdlC, mdlBC, 'pBspk', 'ss_frBspk', ...
%           'grpVar', 'Group');
%       lme_mediationPlot(res, mdlA, mdlBC, tblMea, 'grpVar', 'Group')
%
%   See also: LME_MEDIATION, PLOT_SCAT

%% ========================================================================
%  INPUT
%  ========================================================================

p = inputParser;
addRequired(p, 'res', @isstruct);
addRequired(p, 'mdlA');
addRequired(p, 'mdlBC');
addRequired(p, 'tbl', @istable);
addParameter(p, 'xVar', '', @ischar);
addParameter(p, 'mVar', '', @ischar);
addParameter(p, 'grpVar', '', @ischar);
addParameter(p, 'Parent', [], @(x) isempty(x) || isgraphics(x));
parse(p, res, mdlA, mdlBC, tbl, varargin{:});

xVar   = p.Results.xVar;
mVar   = p.Results.mVar;
grpVar = p.Results.grpVar;
hPar   = p.Results.Parent;

% Default variable names from models
if isempty(mVar)
    mVar = mdlA.ResponseName;
end

% Determine xVar from the paths table Description if not provided
if isempty(xVar)
    % Fall back to first predictor that isn't the mediator
    predNames = mdlA.PredictorNames;
    predNames(strcmp(predNames, mVar)) = [];
    if ~isempty(grpVar)
        predNames(strcmp(predNames, grpVar)) = [];
    end
    xVar = predNames{1};
end

yVar = mdlBC.ResponseName;


%% ========================================================================
%  FIGURE SETUP
%  ========================================================================

if isempty(hPar)
    figure('Color', 'w', 'Position', [100 100 1000 800]);
    hPar = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
else
    if isa(hPar, 'matlab.ui.Figure')
        figure(hPar);
        hPar = tiledlayout(hPar, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    end
end


%% ========================================================================
%  COMPUTE PLOT DATA
%  ========================================================================
% Extract raw data and partial residuals from models. For grouped analyses,
% subset by group membership in the model data.

nGrp = numel(res);
X_all = cell(nGrp, 1);
M_all = cell(nGrp, 1);
Y_all = cell(nGrp, 1);
Y_partM_all = cell(nGrp, 1);
Y_partX_all = cell(nGrp, 1);
G_all = cell(nGrp, 1);

fitA_raw  = fitted(mdlA);
resBC_raw = residuals(mdlBC, 'ResidualType', 'Raw');
yBC_raw   = mdlBC.Variables.(yVar);

for iGrp = 1:nGrp

    grpLevel = res(iGrp).grpLevel;

    % Determine row indices for this group in the model data
    if ~isempty(grpVar) && ~isempty(grpLevel)
        idxA  = mdlA.Variables.(grpVar)  == grpLevel;
        idxBC = mdlBC.Variables.(grpVar) == grpLevel;
        idxTbl = tbl.(grpVar) == grpLevel;
    else
        idxA  = true(mdlA.NumObservations, 1);
        idxBC = true(mdlBC.NumObservations, 1);
        idxTbl = true(height(tbl), 1);
    end

    % Raw data from tbl (may be pre-transformed for plotting)
    X_all{iGrp} = tbl.(xVar)(idxTbl);
    M_all{iGrp} = tbl.(mVar)(idxTbl);

    % Y from model (transformed scale, e.g. log for log-normal)
    Y_all{iGrp} = yBC_raw(idxBC);

    % Extract betaB for partial residuals
    betaB = res(iGrp).paths.Estimate(2);

    % Partial residuals
    rawResid = resBC_raw(idxBC);
    M_mdl    = mdlBC.Variables.(mdlA.ResponseName)(idxBC);
    Y_partM_all{iGrp} = rawResid + betaB * M_mdl;
    Y_partX_all{iGrp} = Y_all{iGrp} - betaB * M_mdl;

    % Group label
    if ~isempty(grpLevel)
        G_all{iGrp} = repmat(string(grpLevel), sum(idxTbl), 1);
    else
        G_all{iGrp} = repmat("", sum(idxTbl), 1);
    end
end

% Concatenate across groups
X = vertcat(X_all{:});
M = vertcat(M_all{:});
Y = vertcat(Y_all{:});
Y_partM = vertcat(Y_partM_all{:});
Y_partX = vertcat(Y_partX_all{:});
G = vertcat(G_all{:});

% Fit type for regression line
fitType = 'ortho';
if ~isnumeric(X)
    fitType = 'None';
end

% Stats annotation helper: combine groups or show per-group
getStats = @(pathRow) arrayfun(@(r) sprintf('%c=%.3f, p=%.4f', ...
    946, r.paths.Estimate(pathRow), r.paths.pValue(pathRow)), ...
    res, 'UniformOutput', false);

fmtTitle = @(pathName, pathRow) [{pathName}, getStats(pathRow)];


%% ========================================================================
%  PATH A (X -> M)
%  ========================================================================
hAx = nexttile(hPar);
plot_scat([], X, M, 'g', G, ...
    'fitType', fitType, 'flgStats', false, 'alpha', 0.5, 'hAx', hAx);
title(fmtTitle('Path A', 1), 'FontWeight', 'normal');
xlabel(xVar, 'Interpreter', 'none');
ylabel(mVar, 'Interpreter', 'none');


%% ========================================================================
%  PATH B (M -> Y | X)
%  ========================================================================
hAx = nexttile(hPar);
plot_scat([], M, Y_partM, 'g', G, ...
    'fitType', fitType, 'flgStats', false, 'alpha', 0.5, 'hAx', hAx);
title(fmtTitle('Path B', 2), 'FontWeight', 'normal');
xlabel(mVar, 'Interpreter', 'none');
ylabel([yVar ' | X'], 'Interpreter', 'none');


%% ========================================================================
%  PATH C' (DIRECT: X -> Y | M)
%  ========================================================================
hAx = nexttile(hPar);
plot_scat([], X, Y_partX, 'g', G, ...
    'fitType', fitType, 'flgStats', false, 'alpha', 0.5, 'hAx', hAx);
title(fmtTitle('Path C'' (Direct)', 4), 'FontWeight', 'normal');
xlabel(xVar, 'Interpreter', 'none');
ylabel([yVar ' | M'], 'Interpreter', 'none');


%% ========================================================================
%  PATH C (TOTAL: X -> Y)
%  ========================================================================
hAx = nexttile(hPar);
plot_scat([], X, Y, 'g', G, ...
    'fitType', fitType, 'flgStats', false, 'alpha', 0.5, 'hAx', hAx);
title(fmtTitle('Path C (Total)', 3), 'FontWeight', 'normal');
xlabel(xVar, 'Interpreter', 'none');
ylabel(yVar, 'Interpreter', 'none');

end
