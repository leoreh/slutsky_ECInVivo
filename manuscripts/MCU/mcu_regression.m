%% ========================================================================
%  MCU REGRESSION MODELS
%  ========================================================================
% Various LME models for analyzing FRH recovery dynamics, speficially the
% relationship between burstiness (BSL), burst- and single-firing (SS).
% Including:
% State Space Trajectory 
% Residuals (partial regression)
% LSmeans (partial dependence)
% Mediation
% Ablation
% 
% Can work on both in vivo and mea tables (tbl)


%% ========================================================================
%  LOAD DATA
%  ========================================================================

% Load with steady state variables
presets = {'steadyState'};
[tbl, xVec, basepaths, v] = mcu_tblMea('presets', presets, 'flgOtl', true);

% Add logit pBspk
tblTrans = tbl_trans(tbl, 'varsInc', {'pBspk', 'ss_pBspk'}, 'logBase', 'logit');
tbl.pBspk_trans = tblTrans.pBspk;
tbl.ss_pBspk_trans = tblTrans.ss_pBspk;


%% ========================================================================
%  PRE-PROCESS
%  ========================================================================

% Relative
tbl.dBrst_rel = log((tbl.ss_frBspk) ./ (tbl.frBspk));
tbl.dSngl_rel = log((tbl.ss_frSspk) ./ (tbl.frSspk));
tbl.dFr = log((tbl.frSs) ./ (tbl.fr));
tbl.dpBspk = (tbl.ss_pBspk_trans) - (tbl.pBspk_trans);

% Absolute
tbl.dBrst_abs = (tbl.ss_frBspk - tbl.frBspk);
tbl.dSngl_abs = (tbl.ss_frSspk - tbl.frSspk);
tbl.dFr_abs = (tbl.frSs - tbl.fr);

% % Prism
% idxGrp = tbl.Group == 'MCU-KO';
% prismVars = {'Name', 'dSngl_rel', 'dBrst_rel', 'dSngl_abs', 'dBrst_abs', ...
%     'dFr', 'dpBspk', 'pBspk', 'pBspk_trans'};
% tblPrism = tbl(idxGrp, prismVars);
% string(prismVars)




%% ========================================================================
%  PLOTTING
%  ========================================================================

hFig = mcu_rcvSpace(tbl);

frml = 'dSngl_rel ~ dBrst_rel * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'verbose', true);

frml = 'dSngl_abs ~ dBrst_abs * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'verbose', true);


% tblGUI_scatHist(tbl, 'xVar', 'pBspk_trans', 'yVar', 'dBrst_rel', 'grpVar', 'Group');
% tblGUI_bar(tbl, 'yVar', 'pBspk', 'xVar', 'Group');

%% ========================================================================
%  RESIDUAL ANALYSIS
%  ========================================================================
% Partial Regression


hFig = figure;
hTile = tiledlayout(hFig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

hAx = nexttile;
frml = 'dSngl_rel ~ (pBspk + fr + dBrst_rel) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'verbose', false);

[tblRes, hFig] = lme_pr(lmeMdl, 'dBrst_rel', ...
    'flgMode', 'regression', 'hAx', hAx, ...
    'varGrp',  'Group', 'transParams', lmeInfo.transParams);

hAx = nexttile;
frml = 'dSngl_abs ~ (pBspk + fr + dBrst_abs) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'verbose', false);

[tblRes, hFig] = lme_pr(lmeMdl, 'dBrst_abs', ...
    'flgMode', 'regression', 'hAx', hAx, ...
    'varGrp',  'Group', 'transParams', lmeInfo.transParams);

% Prism
idxGrp = tblRes.Group == 'Control';
tblRes.ResidY(idxGrp);

% % Obsolete
% hFig = mcu_rcvRes(tbl);








%% ========================================================================
%  PARTIAL DEPENDENCE (LSMEANS)
%  ========================================================================

% Steady-state raw values
dist = 'log-normal';

frml = 'ss_frBspk ~ (fr + pBspk) * Group + (1|Name)';
[lmeMdl1, lmeStats, lmeInfo1, ~] = lme_analyse(tbl, frml, ...
    'dist', dist);

frml = 'ss_frSspk ~ (fr + pBspk) * Group + (1|Name)';
[lmeMdl2, lmeStats, lmeInfo2, ~] = lme_analyse(tbl, frml, ...
    'dist', dist);

frml = 'frSs ~ (fr + pBspk) * Group + (1|Name)';
[lmeMdl3, lmeStats, lmeInfo3, ~] = lme_analyse(tbl, frml, ...
    'dist', dist);

% ---
% Plot Interaction
hFig = plot_axSize('flgFullscreen', true, 'flgPos', true);

vars = {'pBspk', 'Group'};
hAx = nexttile;
[pdRes, hFig] = lme_pd(lmeMdl1, vars, 'transParams', lmeInfo1.transParams, ...
    'hAx', hAx);
set(gca, "YScale", "log")
set(gca, "XScale", "log")

hAx = nexttile;
[pdRes, hFig] = lme_pd(lmeMdl2, vars, 'transParams', lmeInfo2.transParams, ...
    'hAx', hAx);
set(gca, "YScale", "log")

hAx = nexttile;
[pdRes, hFig] = lme_pd(lmeMdl3, vars, 'transParams', lmeInfo3.transParams, ...
    'hAx', hAx);
set(gca, "YScale", "log")

vars = {'fr', 'Group'};
hAx = nexttile;
[pdRes, hFig] = lme_pd(lmeMdl1, vars, 'transParams', lmeInfo1.transParams, ...
    'hAx', hAx);
set(gca, "YScale", "log")
set(gca, "XScale", "log")

hAx = nexttile;
[pdRes, hFig] = lme_pd(lmeMdl2, vars, 'transParams', lmeInfo2.transParams, ...
    'hAx', hAx);
set(gca, "YScale", "log")
set(gca, "XScale", "log")

hAx = nexttile;
[pdRes, hFig] = lme_pd(lmeMdl3, vars, 'transParams', lmeInfo3.transParams, ...
    'hAx', hAx);
set(gca, "YScale", "log")
set(gca, "XScale", "log")




% To Prism
grpIdx = pdRes.Group == "Control";
prismMat = [pdRes(grpIdx, vars(1)), ...
    pdRes(grpIdx, {'frSs_pred', 'frSs_upper', 'frSs_lower'})];
grpIdx = pdRes.Group == "MCU-KO";
prismMat = [pdRes(grpIdx, vars(1)), ...
    pdRes(grpIdx, {'frSs_pred', 'frSs_upper', 'frSs_lower'})];

%% ========================================================================
%  MEDIATION
%  ========================================================================

frml = 'ss_frSspk ~ pBspk + fr + (1|Name)';
xVar = 'pBspk';
mVar = 'ss_frBspk';
distM = 'log-normal';
distY = distM;

% WT
tblWt = tbl(tbl.Group == 'Control', :);
resWt = lme_mediation(tblWt, frml, xVar, mVar, 'distM', distM, 'distY', distY);

resWt.plot.X = tblWt.pBspk_trans;
lme_mediationPlot(resWt)


% MCU
tblMcu = tbl(tbl.Group == 'MCU-KO', :);
resMcu = lme_mediation(tblMcu, frml, xVar, mVar, 'distM', distM, 'distY', distY);

resMcu.plot.X = tblMcu.pBspk_trans;
lme_mediationPlot(resMcu)





%% ========================================================================
%  ABLATION
%  ========================================================================

frml = 'dFr ~ (dBrst_rel + dSngl_rel)';
dist = 'normal';

frml = 'frSs ~ (frBspk + frSspk) + (1 | Name)';
% frml = 'frSs ~ (frBspk + frSspk)';
% frml = 'frSs ~ (fr + pBspk) + (1 | Name)';
dist = 'log-normal';

nRep = 5;
partMode = 'split';

tblWt = tbl(tbl.Group == 'Control', :);
abl = lme_ablation(tblWt, frml, 'dist', dist, ...
    'flgBkTrans', false, 'partitionMode', partMode, 'nrep', nRep);

tblMcu = tbl(tbl.Group == 'MCU-KO', :);
abl = lme_ablation(tblMcu, frml, 'dist', dist, ...
    'flgBkTrans', false, 'partitionMode', partMode, 'nrep', nRep);


% With Group
tblLme = tbl;

frml = 'frSs ~ (frBspk + frSspk) * Group + (1 | Name)';
frml = 'frSs ~ (fr + pBspk) * Group + (1 | Name)';

nRep = 5;
dist = 'log-normal';

abl = lme_ablation(tblLme, frml, 'dist', dist, ...
    'flgBkTrans', false, 'partitionMode', 'batch', 'nrep', nRep);


% tblLme.Group = double(tbl.Group) - 1;
% dist = 'binomial';
% frml = 'Group ~ (frSs + fr + pBspk) + (1 | Name)';

%% ========================================================================
%  CONTRIBUTION ANALYSIS: BSL vs SS
%  ========================================================================
%  Test if the relative contribution of burst spikes (pBspk) changes between
%  Baseline and Steady State (SS), and if this change differs by Group.

fprintf('\n================================================================\n');
fprintf(' STATISTICS: CONTRIBUTION ANALYSIS (pBspk: BSL vs SS)\n');
fprintf('================================================================\n');

% Construct Tall Table

idxGrp = tbl.Group == 'Control';
idxGrp = true(height(tbl), 1);

% Baseline
varsTbl = {'Group', 'Name', 'frBspk', 'frSspk', 'fr'};
tBsl = tbl(idxGrp, varsTbl);
tBsl.Timepoint = repmat({'BSL'}, height(tBsl), 1);

% Steady State
varsSs = {'Group', 'Name', 'ss_frBspk', 'ss_frSspk', 'frSs'};
tSs = tbl(idxGrp, varsSs);
tSs.Properties.VariableNames = varsTbl;
tSs.Timepoint = repmat({'SS'}, height(tSs), 1);

% Concatenate
tblLong = [tBsl; tSs];
tblLong.Timepoint = categorical(tblLong.Timepoint, {'BSL', 'SS'});

% LME Analysis
frml = 'fr ~ (frSspk + frBspk) * Timepoint + (frSspk + frBspk) * Group + (1|Name)';

[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tblLong, frml, ...
    'dist', 'log-normal');


% Per Group
idxGrp = tblLong.Group == 'MCU-KO';
tblGrp = tblLong(idxGrp, :);

frml = 'fr ~ (frSspk + frBspk) * Timepoint + (1|Name)';

[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tblGrp, frml, ...
    'dist', 'log-normal');
