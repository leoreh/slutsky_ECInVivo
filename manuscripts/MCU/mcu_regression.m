%% ========================================================================
%  MCU REGRESSION MODELS
%  ========================================================================
% Various LME models for analyzing FRH recovery dynamics, speficially the
% relationship between burstiness (BSL), burst- and single-firing (SS).
% Including:
% State Space Trajectory 
% Residuals (partial regression)
% LSmeans (least-squares means)
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
tbl.bGain = log((tbl.ss_frBspk) ./ (tbl.frBspk));
tbl.sGain = log((tbl.ss_frSspk) ./ (tbl.frSspk));
tbl.frGain = log((tbl.ss_fr) ./ (tbl.fr));
tbl.pBspkGain = (tbl.ss_pBspk_trans) - (tbl.pBspk_trans);

% Absolute
tbl.bDelta = (tbl.ss_frBspk - tbl.frBspk);
tbl.sDelta = (tbl.ss_frSspk - tbl.frSspk);
tbl.frDelta = (tbl.ss_fr - tbl.fr);



%% ========================================================================
%  EXCHANGE RATE
%  ========================================================================
% Use raw values (not standardized effects)

frml = 'sDelta ~ (bDelta) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'verbose', true, 'flgStnd', false);

frml = 'sGain ~ bGain * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'verbose', true, 'flgStnd', false);

%% ========================================================================
%  PLOTTING
%  ========================================================================

hFig = mcu_rcvSpace(tbl);

frml = 'sGain ~ bGain * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'verbose', true);

frml = 'sDelta ~ bDelta * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'verbose', true);


% tblGUI_scatHist(tbl, 'xVar', 'pBspk_trans', 'yVar', 'bGain', 'grpVar', 'Group');
% tblGUI_bar(tbl, 'yVar', 'pBspk', 'xVar', 'Group');

%% ========================================================================
%  RESIDUAL ANALYSIS
%  ========================================================================
% Partial Regression


hFig = figure;
hTile = tiledlayout(hFig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

hAx = nexttile;
frml = 'sGain ~ (pBspk + fr + bGain) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'verbose', false);

[tblRes, hFig] = lme_pr(lmeMdl, 'bGain', ...
    'flgMode', 'regression', 'hAx', hAx, ...
    'varGrp',  'Group', 'transParams', lmeInfo.transParams);

% frml = 'ResidY ~ ResidX * Group + (1|Name)';
% [lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tblRes, frml, ...
%     'dist', 'normal', 'verbose', true, 'flgStnd', false);

hAx = nexttile;
frml = 'sDelta ~ (pBspk * fr * bDelta) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'verbose', false);

[tblRes, hFig] = lme_pr(lmeMdl, 'bDelta', ...
    'flgMode', 'regression', 'hAx', hAx, ...
    'varGrp',  'Group', 'transParams', lmeInfo.transParams);

% Prism
idxGrp = tblRes.Group == 'Control';
tblRes.ResidY(idxGrp);

% % Obsolete
% hFig = mcu_rcvRes(tbl);








%% ========================================================================
%  LEAST-SQUARES MEANS (LSMEANS)
%  ========================================================================

% Steady-state raw values
dist = 'log-normal';

frml = 'ss_frBspk ~ (fr + pBspk) * Group + (1|Name)';
[lmeMdl1, lmeStats, lmeInfo1, ~] = lme_analyse(tbl, frml, ...
    'dist', dist);

frml = 'ss_frSspk ~ (fr + pBspk) * Group + (1|Name)';
[lmeMdl2, lmeStats, lmeInfo2, ~] = lme_analyse(tbl, frml, ...
    'dist', dist);

frml = 'ss_fr ~ (fr + pBspk) * Group + (1|Name)';
[lmeMdl3, lmeStats, lmeInfo3, ~] = lme_analyse(tbl, frml, ...
    'dist', dist);

% ---
% Plot Interaction
hFig = plot_axSize('flgFullscreen', true, 'flgPos', true);

vars = {'pBspk', 'Group'};
hAx = nexttile;
[pdRes, hFig] = lme_lsmeans(lmeMdl1, vars, 'transParams', lmeInfo1.transParams, ...
    'hAx', hAx);
set(gca, "YScale", "log")

hAx = nexttile;
[pdRes, hFig] = lme_lsmeans(lmeMdl2, vars, 'transParams', lmeInfo2.transParams, ...
    'hAx', hAx);
set(gca, "YScale", "log")

hAx = nexttile;
[pdRes, hFig] = lme_lsmeans(lmeMdl3, vars, 'transParams', lmeInfo3.transParams, ...
    'hAx', hAx);
set(gca, "YScale", "log")

vars = {'fr', 'Group'};
hAx = nexttile;
[pdRes, hFig] = lme_lsmeans(lmeMdl1, vars, 'transParams', lmeInfo1.transParams, ...
    'hAx', hAx);
set(gca, "YScale", "log")
set(gca, "XScale", "log")

hAx = nexttile;
[pdRes, hFig] = lme_lsmeans(lmeMdl2, vars, 'transParams', lmeInfo2.transParams, ...
    'hAx', hAx);
set(gca, "YScale", "log")
set(gca, "XScale", "log")

hAx = nexttile;
[pdRes, hFig] = lme_lsmeans(lmeMdl3, vars, 'transParams', lmeInfo3.transParams, ...
    'hAx', hAx);
set(gca, "YScale", "log")
set(gca, "XScale", "log")



% To Prism
grpIdx = pdRes.Group == "Control";
prismMat = [pdRes(grpIdx, vars(1)), ...
    pdRes(grpIdx, {'ss_fr_pred', 'ss_fr_upper', 'ss_fr_lower'})];
grpIdx = pdRes.Group == "MCU-KO";
prismMat = [pdRes(grpIdx, vars(1)), ...
    pdRes(grpIdx, {'ss_fr_pred', 'ss_fr_upper', 'ss_fr_lower'})];


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

% COMBINED MODELS (unstandardized)

% X -> M
frml = 'ss_frBspk ~ (fr + pBspk) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', dist, 'verbose', true, 'flgStnd', false);

% X -> Y
frml = 'ss_frSspk ~ (fr + pBspk) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tblWt, frml, ...
    'dist', dist, 'verbose', true, 'flgStnd', false);

% X -> Y | M 
frml = 'ss_frSspk ~ (fr + pBspk + ss_frBspk) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'log-normal', 'verbose', true, 'flgStnd', false);



%% ========================================================================
%  ABLATION
%  ========================================================================


frml = 'ss_fr ~ (frBspk + frSspk) + (1 | Name)';
% frml = 'ss_fr ~ (frBspk + frSspk)';
% frml = 'ss_fr ~ (fr + pBspk) + (1 | Name)';
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

frml = 'ss_fr ~ (frBspk + frSspk) * Group + (1 | Name)';
frml = 'ss_fr ~ (fr + pBspk) * Group + (1 | Name)';

nRep = 5;
dist = 'log-normal';

abl = lme_ablation(tblLme, frml, 'dist', dist, ...
    'flgBkTrans', false, 'partitionMode', 'batch', 'nrep', nRep);


% tblLme.Group = double(tbl.Group) - 1;
% dist = 'binomial';
% frml = 'Group ~ (ss_fr + fr + pBspk) + (1 | Name)';

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
varsSs = {'Group', 'Name', 'ss_frBspk', 'ss_frSspk', 'ss_fr'};
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
