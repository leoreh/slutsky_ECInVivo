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

varClu = 'fr';

% Set manualEdges to fixed inner boundaries (same for all groups).
% Example: [0.1, 0.25] → 3 clusters: <0.1, 0.1-0.25, >0.25
% Leave empty [] to use percentile-based clustering (nClu, alpha).
manualEdges = [];

nClu  = 1;    % number of clusters (percentile mode only)
alpha = 2;    % percentile spacing exponent (percentile mode only)

% Initialize
tblPlot.cluLbl = strings(height(tblPlot), 1);
grps = unique(tblPlot.genotype);

for iGrp = 1:length(grps)
    idxGrp    = tblPlot.genotype == grps(iGrp);
    grpData   = tblPlot.(varClu)(idxGrp);
    idxGlobal = find(idxGrp);

    if ~isempty(manualEdges)
        edges = [-Inf, sort(manualEdges(:)'), Inf];
    else
        % Percentile edges computed per group; inner boundaries only
        p     = linspace(0, 1, nClu + 1) .^ alpha;
        pcts  = sort(prctile(grpData, 100 * p));
        edges = [-Inf, pcts(2:end-1), Inf];
    end

    for iClu = 1:length(edges) - 1
        edgeLo = edges(iClu);
        edgeHi = edges(iClu + 1);
        idxClu = grpData > edgeLo & grpData <= edgeHi;

        if isinf(edgeLo)
            lbl = sprintf('P%d (<%.2g)',       iClu, edgeHi);
        elseif isinf(edgeHi)
            lbl = sprintf('P%d (>%.2g)',       iClu, edgeLo);
        else
            lbl = sprintf('P%d (%.2g-%.2g)',   iClu, edgeLo, edgeHi);
        end

        tblPlot.cluLbl(idxGlobal(idxClu)) = lbl;
    end
end

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
tVars = tblVars(startsWith(tblVars, 't_'));
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

meanType = 'geometric';     % 'geometric' or 'arithmetic'
idxGrp   = tblPlot.genotype == 'MCU-KO';

prismMat = [];
for iClu = 1:length(edges) - 1

    patLbl   = sprintf('P%d (', iClu);
    idxClu   = contains(string(tblPlot.cluLbl), patLbl);
    prismData = tblPlot.t_frTot(idxGrp & idxClu, :)';

    if strcmp(meanType, 'geometric')
        logData  = log(prismData);
        n        = sum(~isnan(logData), 2);
        mLog     = mean(logData, 2, 'omitnan');
        semLog   = std(logData, 0, 2, 'omitnan') ./ sqrt(n);
        mu       = exp(mLog);
        lo       = exp(mLog - semLog);
        hi       = exp(mLog + semLog);
    else
        n        = sum(~isnan(prismData), 2);
        mu       = mean(prismData, 2, 'omitnan');
        sem      = std(prismData, 0, 2, 'omitnan') ./ sqrt(n);
        lo       = mu - sem;
        hi       = mu + sem;
    end

    % Append columns: [Mean, Upper, Lower]
    prismMat = [prismMat, mu, hi, lo];
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