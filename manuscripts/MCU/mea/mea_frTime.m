%% ========================================================================
%  MEA FIRING RATE OVER TIME (WRAPPER)
%  ========================================================================
%  Script to load MEA data, cluster units by baseline firing rate, normalize
%  temporal dynamics, and visualize using tblGUI_xy.

% Load
[tbl, xVec, basepaths, v] = mcu_tblMea('presets', {'time', 'rcv', 'steadyState'});
tblPlot = tbl;

% Add logit pBurst
tblTrans = tbl_trans(tblPlot, 'varsInc', {'pBurst'}, 'logBase', 'logit');
tblPlot.pBurst_trans = tblTrans.pBurst;

%% ========================================================================
%  CLUSTERS
%  ========================================================================
%  Cluster units into percentiles based on a specific variable

varClu = 'pBurst';

% --- Manual boundaries (leave empty [] to use percentile mode) ---
% Defines fixed edges applied identically to both groups.
% Example: [0.1, 0.25] creates 3 clusters: <0.1, 0.1-0.25, >0.25
manualEdges = [0.2];
% manualEdges = [];

% --- Percentile mode (used when manualEdges is empty) ---
nClu   = 3;           % Number of clusters (percentiles)
alpha  = 2;           % Scaling factor for percentile spacing

% Initialize Cluster Label Column
tblPlot.cluLbl = strings(height(tblPlot), 1);

% Build edges
if ~isempty(manualEdges)
    % Manual mode: same edges for all groups
    edges = [-Inf, sort(manualEdges(:)'), Inf];
    nClu  = length(edges) - 1;
    useManual = true;
else
    useManual = false;
end

% Get Unique Groups
grps = unique(tblPlot.genotype);

for iGrp = 1:length(grps)

    idxGrp  = tblPlot.genotype == grps(iGrp);
    grpData = tblPlot.(varClu)(idxGrp);

    if useManual
        percEdges = edges;
    else
        % Percentile mode: edges computed per group
        p = linspace(0, 1, nClu + 1) .^ alpha;
        percEdges = [-Inf, prctile(grpData, 100 * (1 - p(2:end-1))), Inf];
        percEdges = sort(percEdges);
    end

    % Assign Clusters
    for iClu = 1:nClu
        edgeLo = percEdges(iClu);
        edgeHi = percEdges(iClu + 1);

        idxClu = grpData > edgeLo & grpData <= edgeHi;

        % Create Label
        if isinf(edgeLo)
            lbl = sprintf('P%d (<%.2g)', iClu, edgeHi);
        elseif isinf(edgeHi)
            lbl = sprintf('P%d (>%.2g)', iClu, edgeLo);
        else
            lbl = sprintf('P%d (%.2g-%.2g)', iClu, edgeLo, edgeHi);
        end

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

% Floor Value: 1 event per Total Recording Duration (approx 0.0001 Hz)
% Prevents clipping of valid low-rate dynamics when using Geometric Mean
floorVal = 1 / (max(xVec) * 3600);

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
    'tileVar', 'genotype', ...    % Separate tiles by Group (Control vs KO)
    'tileFlow', 'vertical', ...
    'xLbl', 'Time (Hours)');


%% ========================================================================
%  PRISM
%  ========================================================================

% Loop over clusters and calculate geometric stats for each
idxGrp = tblPlot.genotype == 'Control';

% Grab raw matrix
tblPlot.t_frTot(idxGrp, :)';

prismMat = [];
for iClu = 1:nClu

    % Extract Data for specific Cluster within Group
    % Note: cluLbl format is "P<d> (<num>-<num>)"
    patLbl  = sprintf('P%d (', iClu);
    idxClu  = contains(string(tblPlot.cluLbl), patLbl);
    finalIdx = idxGrp & idxClu;

    prismData = tblPlot.t_frTot(finalIdx, :)';

    % Log-transform (floorVal already applied in tbl_tNorm)
    logData = log(prismData);

    n    = sum(~isnan(logData), 2);
    mLog = mean(logData, 2, 'omitnan');
    sLog = std(logData, 0, 2, 'omitnan') ./ sqrt(n);

    geoMean = exp(mLog);
    lower   = exp(mLog - sLog);
    upper   = exp(mLog + sLog);

    % Append to matrix: [Mean, Upper, Lower]
    prismMat = [prismMat, geoMean, upper, lower];
end

%% ========================================================================
%  AGGREGATE BY NAME
%  ========================================================================
%  Create a summary table where each unit (row) is the average of all units
%  belonging to the same 'Name' (Animal).

% % Define Grouping Variables
% grpVars = {'sbjID', 'genotype'};
% 
% % Identify Numeric Variables to Average (Time-Series columns & others)
% % We specifically target time-series variables starting with 't_'
% tVars = tblPlot.Properties.VariableNames(contains(tblPlot.Properties.VariableNames, 't_'));
% 
% % Aggregate
% tblPlot = groupsummary(tblPlot, grpVars, 'mean', tVars);
% 
% % Cleanup: Remove 'mean_' prefix from variable names
% for iVar = 1:length(tVars)
%     oldName = ['mean_' tVars{iVar}];
%     if ismember(oldName, tblPlot.Properties.VariableNames)
%         tblPlot.Properties.VariableNames{oldName} = tVars{iVar};
%     end
% end
% 
% % Check
% disp('Table aggregated by Name. Rows:');
% disp(height(tblPlot));
% 
% tblGUI_xy(xVec, tblPlot, ...
%     'yVar', 't_fr', ...
%     'grpVar', 'cluLbl', ...    % Group lines by Cluster
%     'tileVar', 'genotype', ...    % Separate tiles by Group (Control vs KO)
%     'tileFlow', 'vertical', ...
%     'xLbl', 'Time (Hours)');
% 
% idxGrp = tblPlot.genotype == 'MCU-KO';
% prismMat = tblPlot.t_frTot(idxGrp, :)';