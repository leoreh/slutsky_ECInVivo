%% ========================================================================
%  MEA FIRING RATE OVER TIME (WRAPPER)
%  ========================================================================
%  Script to load MEA data, cluster units by baseline firing rate, normalize
%  temporal dynamics, and visualize using tblGUI_xy.

%% ========================================================================
%  LOAD
%  ========================================================================
clear; clc;

% Load table with temporal dynamics ('time' preset)
% [tbl, xVec, basepaths, v] = mcu_tblMea('presets', {'time'});

% Correct time vector (Max data skips every other hour after perturbation)
% Reference: manuscripts/MCU/mea/obsolete/mcu_meaFRt.m
xVec(xVec > 0) = xVec(xVec > 0) * 2;


%% ========================================================================
%  CLUSTERS
%  ========================================================================
%  Cluster units into percentiles based on a specific variable (fr)

nClu   = 3;         % Number of clusters (percentiles)
varClu = 'fr';      % Variable to cluster by ('fr' = Baseline Firing Rate)
alpha  = 1.5;       % Scaling factor for percentile spacing

% Initialize Cluster Label Column
tbl.cluLbl = strings(height(tbl), 1);

% Get Unique Groups
grps = unique(tbl.Group);

for iGrp = 1:length(grps)

    idxGrp = tbl.Group == grps(iGrp);

    % Extract Data for Clustering
    grpData = tbl.(varClu)(idxGrp);

    % Calculate Percentile Edges
    % Logic adapted from mcu_meaFRt.m
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
        lbl = sprintf('P%d (%.1f-%.1f Hz)', iClu, edgeLo, edgeHi);

        % Map back to full table
        % Find indices in full table that match group AND cluster
        % (Subset of idxGrp)
        idxGlobal = find(idxGrp);
        tbl.cluLbl(idxGlobal(idxClu)) = lbl;
    end
end

% Convert to categorical for GUI grouping
tbl.cluLbl = categorical(tbl.cluLbl);

%% ========================================================================
%  NORMALIZE
%  ========================================================================

%  Normalize traces to baseline percentage

% Define Baseline Window (Indices where Time < 0)
winNorm = [1, find(xVec >= 0, 1) - 1];

% Normalize using tbl_tNorm
tbl = tbl_tNorm(tbl, 'varsInc', {'t_fr'}, 'winNorm', winNorm, ...
    'Method', 'percentage', 'varsGrp', 'Name');

% Convert Ratio to Percentage
% tbl.t_fr = tbl.t_fr * 100;

%% ========================================================================
%  PLOT
%  ========================================================================

tblGUI_xy(xVec, tbl, ...
    'yVar', 't_fr', ...
    'grpVar', 'cluLbl', ...    % Group lines by Cluster
    'tileVar', 'Group', ...    % Separate tiles by Group (Control vs KO)
    'tileFlow', 'flow', ...
    'xLbl', 'Time (Hours)');
