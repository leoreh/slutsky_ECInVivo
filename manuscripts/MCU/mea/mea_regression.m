
%% ========================================================================
%  LOAD AND PREP
%  ========================================================================

presets = {'frNet'};
[tbl, xVec, basepaths, v] = mcu_tblMea('presets', presets);

% tblGUI_xy(xVec, tbl);
% tblGUI_scatHist(tbl, 'xVar', 'pBspk', 'yVar', 'funcon', 'grpVar', 'Group');
% tblGUI_bar(tbl, 'yVar', 'pBspk', 'xVar', 'Group');
% tblGUI_raster(tbl, 'grpVar', 'Name', 'grpVal', 'ctrl1')

tblLme = tbl;

% Clean units that were not perturbed
tblLme(~tblLme.uPert, :) = [];
tblLme = removevars(tblLme, 'uPert');

% Add logit pBspk
tblTrans = tbl_trans(tblLme, 'varsInc', {'pBspk'}, 'logBase', 'logit');
tblLme.pBspk_trans = tblTrans.pBspk;



tblGUI_scatHist(tblLme, 'xVar', 'pBspk', 'yVar', 'funcon', 'grpVar', 'Group');
tblGUI_bar(tblLme, 'yVar', 'pBspk', 'xVar', 'Group');

grpIdx = tblLme.Group == "MCU-KO";
tblLme(grpIdx, {'pBspk_trans', 'funcon'})

% Exclude non-bursting
tblLme = tblLme(tblLme.pBspk > 0, :);



%% ========================================================================
%  ABLATION
%  ========================================================================

frml = 'frSs ~ (frBspk + frSspk)';

tblWt = tblLme(tblLme.Group == 'Control', :);

abl = lme_ablation(tblWt, frml, 'dist', 'gamma', ...
    'flgBkTrans', false, 'partitionMode', 'split');

tblMcu = tblLme(tblLme.Group == 'MCU-KO', :);

abl = lme_ablation(tblMcu, frml, 'dist', 'log-normal', ...
    'flgBkTrans', false, 'partitionMode', 'split');

%% ========================================================================
%  Full Model
%  ========================================================================

xVar = 'pBspk';

frml = sprintf('frSs ~ (fr + %s) * Group + (1|Name)', xVar);

[lmeMdl, lmeStats, lmeInfo, tblMdl] = lme_analyse(tblLme, frml, 'dist', 'gamma');


% PLOT INTERACTION
vars = {xVar, 'Group'};
hFig = plot_axSize('flgFullscreen', true, 'flgPos', true);

% Partial Residuals
hAx = nexttile;
[prRes, hFig] = lme_pr(lmeMdl, vars, 'transParams', lmeInfo.transParams, ...
    'hAx', hAx);

% Partial Dependence
hAx = nexttile;
[pdRes, hFig] = lme_pd(lmeMdl, vars, 'transParams', lmeInfo.transParams, ...
    'hAx', hAx);






% To Prism
grpIdx = pdRes.Group == "MCU-KO";
[pdRes(grpIdx, {'pBspk_trans'}), ...
    pdRes(grpIdx, {'frSs_pred', 'frSs_upper', 'frSs_lower'})]








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
% Gamma distribution fails to converge in the ablation analysis
% =========================================================================








