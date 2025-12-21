



%% ========================================================================
%  BOUT ANALYSIS (Baseline)
%  ========================================================================

% Files
basepaths = [mcu_basepaths('wt_bsl'), mcu_basepaths('wt_bsl_ripp'),...
    mcu_basepaths('mcu_bsl')];
basepaths = natsort(unique(basepaths));

% Load
vars = {'sleep_states'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

% Config
cfg = mcu_cfg();
cfg.lbl.states = v(1).ss.info.names;
clear varMap
varMap.BoutLen = 'bouts.blen';
sStates = [1, 4, 5];

% Pre-process: extract bout length per state
for iState = sStates
    for iFile = 1:length(basepaths)
        v(iFile).ss.bouts.blen = v(iFile).ss.bouts.boutLen{iState};
    end

    % TABLE

    % Metadata
    tagFiles = struct();
    tagFiles.Name = get_mname(basepaths);
    [~, fileNames] = fileparts(basepaths);
    tagFiles.File = fileNames;
    tagAll.State = cfg.lbl.states{iState};

    % Table
    tbl = v2tbl('v', [v(:).ss], 'varMap', varMap, 'tagAll',...
        tagAll, 'tagFiles', tagFiles, 'idxCol', []);

    % Group
    tbl.Group = ones(height(tbl), 1) * 1;
    tbl.Group(ismember(tbl.Name, cfg.miceMCU), :) = 2;
    tbl.Group = categorical(tbl.Group, [1, 2], cfg.lbl.grp);

    tblCell{iState} = tbl;
end

tblSs = vertcat(tblCell{:});

% Reorder columns
varOrder = {'Group', 'Name', 'File', 'UnitID', 'State'};
tblSs = movevars(tblSs, varOrder, 'Before', 1);

% Plot
hFig = tblGUI_bar(tblSs, 'yVar', 'BoutLen', 'xVar', 'State', 'GrpVar', 'Group');



% -------------------------------------------------------------------------
% LME

% Formula
frml = 'BoutLen ~ Group * State + (1|Name)';

% Check best model
% statsPark = lme_parkTest(tblSs, frml)
% statsDist = lme_compareDists(tblSs, frml)

% Run LME
% Comparison reveals InverseGaussian
cfgLme.dist = 'InverseGaussian';
cfgLme.contrasts = [1 : 9];
[lmeStats, lmeMdl] = lme_analyse(tblSs, frml, cfgLme);


% Prism
idxRow = tblSs.State == 'WAKE';
prismMat = tbl2prism(tblSs(idxRow, :), 'yVar', 'BoutLen');






%% ========================================================================
%  STATE DURATION (Baseline)
%  ========================================================================

% Files
basepaths = [mcu_basepaths('wt_bsl'), mcu_basepaths('wt_bsl_ripp'),...
    mcu_basepaths('mcu_bsl')];
basepaths = natsort(unique(basepaths));

% Load
vars = {'sleep_states'};
v = basepaths2vars('basepaths', basepaths, 'vars', vars);

% Config
cfg = mcu_cfg();
cfg.lbl.states = v(1).ss.info.names;
clear varMap
varMap.StatePrct = 'bouts.StatePrct';
sStates = [1, 4, 5];

% Pre-process
for iState = sStates
    
    for iFile = 1:length(basepaths)
        v(iFile).ss.bouts.StatePrct = v(iFile).ss.bouts.prctDur(1, iState);
    end

    % Metadata
    tagFiles = struct();
    tagFiles.Name = get_mname(basepaths);
    [~, fileNames] = fileparts(basepaths);
    tagFiles.File = fileNames;
    tagAll.State = cfg.lbl.states{iState};

    % Table
    tbl = v2tbl('v', [v(:).ss], 'varMap', varMap, 'tagAll',...
        tagAll, 'tagFiles', tagFiles, 'idxCol', iState);

    % Group
    tbl.Group = ones(height(tbl), 1) * 1;
    tbl.Group(ismember(tbl.Name, cfg.miceMCU), :) = 2;
    tbl.Group = categorical(tbl.Group, [1, 2], cfg.lbl.grp);

    tblCell{iState} = tbl;
end

tblSs = vertcat(tblCell{:});

% Reorder columns
varOrder = {'Group', 'Name', 'File', 'UnitID', 'State'};
tblSs = movevars(tblSs, varOrder, 'Before', 1);

% Plot
hFig = tblGUI_bar(tblSs, 'yVar', 'StatePrct', 'xVar', 'State', 'GrpVar', 'Group');



% -------------------------------------------------------------------------
% LME

% Formula
frml = 'StatePrct ~ Group * State + (1|Name)';

% Check best model
% statsPark = lme_parkTest(tblSs, frml)
% statsDist = lme_compareDists(tblSs, frml)

% Run LME
% Comparison reveals InverseGaussian
cfgLme.dist = 'Gamma';
cfgLme.contrasts = [1 : 9];
[lmeStats, lmeMdl] = lme_analyse(tblSs, frml, cfgLme);


% Prism
idxRow = tblSs.State == 'REM';
prismMat = tbl2prism(tblSs(idxRow, :), 'yVar', 'StatePrct');




