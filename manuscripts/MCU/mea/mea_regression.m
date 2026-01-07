
%% ========================================================================
%  LOAD AND PREP
%  ========================================================================

[tbl, xVec, basepaths, v] = mcu_tblMea(basepaths, v);

% tblGUI_xy(xVec, tbl);
% tblGUI_scatHist(tbl, 'xVar', 'bFrac', 'yVar', 'rcvTime', 'grpVar', 'Group');
% tblGUI_bar(tbl, 'yVar', 'bFrac', 'xVar', 'Group');
% tblGUI_raster(tbl, 'grpVar', 'Name', 'grpVal', 'ctrl1')


% -------------------------------------------------------------------------
% PREPS

tblLme = tbl;
tblLme.Genotype = tblLme.Group == "Control";

% Recovered units
idxUnits = tbl.uRcv;

% List of possible predictors
listPrdct = {'pertDepth', 'fr', 'bFrac', 'Group', '(1|Name)'};
listRspns = {'uRcv', 'rcvGain', 'rcvErr', 'rcvTime', 'spkDfct', 'rcvSlope', 'Genotype'};


%% ========================================================================
%  ABLATION
%  ========================================================================

nFolds = 4;
nReps = 10;

% Recovery
varsFxd  = {'pertDepth', 'bFrac', 'fr', 'Group'};
frml = sprintf('uRcv ~ %s + (1|Name)', strjoin(varsFxd, ' + '));
res = lme_ablation(tblLme, frml, 'nFolds', nFolds, 'nReps', nReps);
mdl = lme_analyse(tblLme, frml);

% Steady-state firing (w/o fr, which is too predictive)
varsFxd  = {'pertDepth', 'bFrac', 'Group'};
frml = sprintf('frSs ~ %s + (1|Name)', strjoin(varsFxd, ' + '));
res = lme_ablation(tblLme, frml, 'nFolds', nFolds, 'nReps', nReps);

% Predict Genotype
varsFxd  = {'pertDepth', 'bFrac', 'fr'};
frml = sprintf('Genotype ~ %s', strjoin(varsFxd, ' + '));
res = lme_ablation(tblLme, frml, 'nFolds', nFolds, 'nReps', nReps);



frml = sprintf('uRcv ~ %s + (1|Name)', strjoin(varsFxd, ' + '));
res = lme_mediation(tblLme, frml, 'Group', 'bFrac');





%% ========================================================================
%  NOTE: DATA SPECIFICS
%  ========================================================================
%
% Insisted to use rcvGain and spkDfct on logarithm scale so they can be
% analyzed with fitlme. rcvTime (and MF) probably still needs glme.
%
% Recovery slope depends too much on the selected model (e.g. sigmoid vs.
% exponential).
%
% Recovery time includes units that reached their threshold value but
% didn't manage to maintain it.
%
% Predicting Genotype is only meaningful without the random effect of Name.
%  ========================================================================








