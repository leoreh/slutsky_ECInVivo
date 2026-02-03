%% ========================================================================
%  MEA FIRING RATE OVER TIME (WRAPPER)
%  ========================================================================
%  Script to load MEA data, cluster units by baseline firing rate, normalize
%  temporal dynamics, and visualize using tblGUI_xy.

% Load
% [tbl, xVec, basepaths, v] = mcu_tblMea('presets', {'time', 'rcv'});
tblPlot = tbl;

% Add logit pBspk
tblTrans = tbl_trans(tblPlot, 'varsInc', {'pBspk'}, 'logBase', 'logit');
tblPlot.pBspk_trans = tblTrans.pBspk;

%% ========================================================================
%  CLUSTERS
%  ========================================================================
%  Cluster units into percentiles based on a specific variable

varClu = 'pBspk';
nClu   = 1;           % Number of clusters (percentiles)
alpha  = 2.5;         % Scaling factor for percentile spacing

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
winNorm = [0, find(xVec >= 0, 1) - 1];
floorVal = 1 / (median(diff(xVec)) * 3600);

tblVars = tblPlot.Properties.VariableNames;
tVars = tblVars(contains(tblVars, 't_'));
tblPlot = tbl_tNorm(tblPlot, 'varsInc', tVars, 'winNorm', winNorm, ...
    'Method', 'percentage', 'flgGeom', true, 'floorVal', floorVal, 'varsGrp', {});


%% ========================================================================
%  PLOT
%  ========================================================================

tblGUI_xy(xVec, tblPlot, ...
    'yVar', 't_fr', ...
    'grpVar', 'cluLbl', ...    % Group lines by Cluster
    'tileVar', 'Group', ...    % Separate tiles by Group (Control vs KO)
    'tileFlow', 'vertical', ...
    'xLbl', 'Time (Hours)');
