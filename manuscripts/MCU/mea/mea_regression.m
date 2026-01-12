
%% ========================================================================
%  LOAD AND PREP
%  ========================================================================

presets = {'rcv', 'frNet'};
[tbl, xVec, basepaths, v] = mcu_tblMea('presets', presets);
tblLme = tbl;

% Clean units that were not perturbed
tblLme(~tblLme.uPert, :) = [];
tblLme = removevars(tblLme, 'uPert');

% tblGUI_xy(xVec, tbl);
% tblGUI_scatHist(tbl, 'xVar', 'pBspk', 'yVar', 'conn', 'grpVar', 'Group');
% tblGUI_bar(tbl, 'yVar', 'pBspk', 'xVar', 'Group');
% tblGUI_raster(tbl, 'grpVar', 'Name', 'grpVal', 'ctrl1')


% List of possible vars
listPrdct = {'pertDepth', 'fr', 'pBspk', 'Group', '(1|Name)'};
listRspns = {'uRcv', 'rcvBsl', 'rcvTime', 'spkDfct', 'Genotype', 'frSs', 'rcvWork'};

%% ========================================================================
%  ABLATION
%  ========================================================================


% Recovery
frml = 'frSs ~ frTrough + fr * Group + pBspk * Group + (1|Name)';

res = lme_ablation(tblLme, frml, 'dist', 'log-normal');


frml = 'frSs ~ frTrough + fr * Group + pBspk * Group + conn * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, tblMdl] = lme_analyse(tblLme, frml, 'dist', 'log-normal');

% Partial Dependence
vars = {'Group', 'pBspk'};
hFig = plot_axSize('flgFullscreen', true, 'flgPos', true);
hAx = nexttile;
[pdRes, hFig] = lme_pd(lmeMdl, vars, 'transParams', lmeInfo.transParams, ...
    'hAx', hAx);

% Partial Residuals
hAx = nexttile;
[pdRes, hFig] = lme_pr(lmeMdl, vars, 'transParams', lmeInfo.transParams, ...
    'hAx', hAx);





tblGUI_scatHist(tblMdl, 'xVar', 'pBspk', 'yVar', 'frSs', 'grpVar', 'Group');




%% ========================================================================
%  OTHER
%  ========================================================================


% Predict Genotype (w/o Name as random effect)
tblLme.Genotype = tblLme.Group == "Control";
varsFxd  = {'pertDepth', 'pBspk', 'fr'};
frml = sprintf('Genotype ~ %s', strjoin(varsFxd, ' + '));
res = lme_ablation(tblLme, frml, 'nFolds', nFolds, 'nReps', nReps);






% Mediation
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








