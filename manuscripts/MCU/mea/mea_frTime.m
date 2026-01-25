%% ========================================================================
%  MEA FIRING RATE OVER TIME (WRAPPER)
%  ========================================================================
%  Script to load MEA data, cluster units by baseline firing rate, normalize
%  temporal dynamics, and visualize using tblGUI_xy.

% Load
% [tbl, xVec, basepaths, v] = mcu_tblMea('presets', {'time'});
tblPlot = tbl;


%% ========================================================================
%  CLUSTERS
%  ========================================================================
%  Cluster units into percentiles based on a specific variable

varClu = 'pBspk';
nClu   = 1;         % Number of clusters (percentiles)
alpha  = 2;         % Scaling factor for percentile spacing

% Initialize Cluster Label Column
tblPlot.cluLbl = strings(height(tblPlot), 1);

% Get Unique Groups
grps = unique(tblPlot.Group);

for iGrp = 1:length(grps)

    idxGrp = tblPlot.Group == grps(iGrp);

    % Extract Data for Clustering
    grpData = tblPlot.(varClu)(idxGrp);

    % Calculate Percentile Edges
    p = linspace(0, 1, nClu + 1) .^ alpha;
    percEdges = prctile(grpData, 100 * (1 - p));
    percEdges = sort(percEdges);

    % Assign Clusters
    for iClu = 1:nClu
        edgeLo = percEdges(iClu);
        edgeHi = percEdges(iClu + 1);

        if iClu == 1
            idxClu = grpData <= edgeHi;
        else
            idxClu = grpData > edgeLo & grpData <= edgeHi;
        end

        % Create Label
        lbl = sprintf('P%d (%.1f-%.1f)', iClu, edgeLo, edgeHi);

        % Map back to full table
        idxGlobal = find(idxGrp);
        tblPlot.cluLbl(idxGlobal(idxClu)) = lbl;
    end
end

% Convert to categorical for GUI grouping
tblPlot.cluLbl = categorical(tblPlot.cluLbl);

%% ========================================================================
%  NORMALIZE
%  ========================================================================
%  Normalize traces to baseline percentage using tbl_tNorm

% Define Baseline Window (Indices where Time < 0)
winNorm = [1, find(xVec >= 0, 1) - 1];

tblVars = tbl.Properties.VariableNames;
tVars = tblVars(contains(tblVars, 't_'));
tblPlot = tbl_tNorm(tblPlot, 'varsInc', tVars, 'winNorm', winNorm, ...
    'Method', 'percentage', 'varsGrp', {});

% Convert Ratio to Percentage
% tblPlot.t_fr = tblPlot.t_fr * 100;


%% ========================================================================
%  CLEAN OUTLIERS
%  ========================================================================
%  Remove units with anomalous temporal dynamics using PCA.
%  Units are identified as outliers if their trajectory in PC-space
%  deviates significantly from the population (Mahalanobis distance).

otl = false(height(tblPlot), 1);
for iVar = 1:length(tVars)
    vName = tVars{iVar};
    data = tblPlot.(vName);

    % Use nearest neighbor to fill gaps so PCA can run
    dataClean = fillmissing(data, 'nearest', 2);

    % centering is automatic in pca()
    [~, score] = pca(dataClean);

    % Detect Outliers in Score Space (Top 3 PCs)
    % Use Mahalanobis distance to find units far from the cluster center
    otlVar = all(isoutlier(score(:, 1:3), 'median', 'ThresholdFactor', 100), 2);
    otl = otl | otlVar;
end
fprintf('Removed %d units.\n', sum(otl));
tblPlot(otl, :) = [];




%% ========================================================================
%  PLOT
%  ========================================================================

tblGUI_xy(xVec, tblPlot, ...
    'yVar', 't_fr', ...
    'grpVar', 'cluLbl', ...    % Group lines by Cluster
    'tileVar', 'Group', ...    % Separate tiles by Group (Control vs KO)
    'tileFlow', 'vertical', ...
    'xLbl', 'Time (Hours)');
