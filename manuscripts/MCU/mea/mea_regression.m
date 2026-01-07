
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
% BslFr is positively correlated with recovery time (model and model free).
% This makes sense considering most units drop to zero, and thus despite
% normalizing the target value, it is still largely influenced by BslFr. A
% similar result is obtained for SpkDftc. Because of this, and because it
% does not predict uRcv, it is omitted from subsequent models.
%
% Recovery slope depends too much on the selected model (e.g. sigmoid vs.
% exponential).
%
% A discripency between model-based and model-free parameters is that in
% the latter there is no correlation between frBsl and pertDepth because
% many values are clamped to c. Hence pertDepth should only be from the
% model.
%
% Recovery time includes units that reached their threshold value but
% didn't manage to maintain it.
%
% Predicting Genotype is only meaningful without the random effect of Name.
%  ========================================================================








