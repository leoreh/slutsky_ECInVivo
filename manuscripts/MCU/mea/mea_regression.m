
%% ========================================================================
%  LOAD AND PREP
%  ========================================================================

presets = {'rcv'};
[tbl, xVec, basepaths, v] = mcu_tblMea('presets', presets);
tblLme = tbl;

% tblGUI_xy(xVec, tbl);
% tblGUI_scatHist(tbl, 'xVar', 'pBspk', 'yVar', 'rcvTime', 'grpVar', 'Group');
% tblGUI_bar(tbl, 'yVar', 'pBspk', 'xVar', 'Group');
% tblGUI_raster(tbl, 'grpVar', 'Name', 'grpVal', 'ctrl1')

% List of possible predictors
listPrdct = {'pertDepth', 'fr', 'pBspk', 'Group', '(1|Name)'};
listRspns = {'uRcv', 'rcvBsl', 'rcvTime', 'spkDfct', 'Genotype', 'frSs', 'rcvWork'};

%% ========================================================================
%  ABLATION
%  ========================================================================


% Recovery
frml = 'frSs ~ frTrough + fr * Group + pBspk * Group + (1|Name)';

res = lme_ablation(tblLme, frml, 'dist', 'log-normal');

[uniqueNames, firstIdx] = unique(tblLme.Name, 'stable');
groupLabels = tblLme.Group(firstIdx);

[lmeMdl, lmeStats, lmeInfo, tblMdl] = lme_analyse(tblLme, frml, 'dist', 'log-normal');






[pd, grp, pBspk_pred] = partialDependence(mdl, {'Group', 'pBspk'});

hFig = figure;
plot(pBspk_pred, pd)

plotPartialDependence(mdl, {'Group', 'pBspk'})

% Predict Genotype (w/o Name as random effect)
tblLme.Genotype = tblLme.Group == "Control";
varsFxd  = {'pertDepth', 'pBspk', 'fr'};
frml = sprintf('Genotype ~ %s', strjoin(varsFxd, ' + '));
res = lme_ablation(tblLme, frml, 'nFolds', nFolds, 'nReps', nReps);






% Medation
frml = sprintf('uRcv ~ %s + (1|Name)', strjoin(varsFxd, ' + '));
res = lme_mediation(tblLme, frml, 'Group', 'pBspk');

% Limit tbl to control
tblCtrl = tblLme(tblLme.Group == "Control", :);
tblCtrl(:, 'Group') = [];
frml = 'frSs ~ frTrough + fr + pBspk + (1|Name)';

res = lme_ablation(tblCtrl, frml, 'dist', 'log-normal');


varsFxd  = {'pertDepth', 'fr', 'pBspk'};
frml = sprintf('rcvBsl ~ %s + (1|Name)', strjoin(varsFxd, ' + '));
res = lme_ablation(tblLme, frml, 'nFolds', nFolds, 'nReps', nReps);




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
% 
% If rcvBsl is the response and fr is a predictor, the results are exactly
% the same as when frSs is the response and fr is not a predictor
% (obvsiouly)
% 
% Burstiness shows stronger predictive capacity when frTrough is provided
% instead of pertDepth
% 
% Burstiness shows stronger predictive capacity when NOT log- or
% logit-transformed. It's skeweness is only ~1.6 so it's okay. 
% 
% Gamma distribution fails to converge in the ablation analysis
%  ========================================================================








