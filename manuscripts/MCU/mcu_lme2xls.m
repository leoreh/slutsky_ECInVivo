
% File Names
pathName = 'D:\OneDrive - Tel-Aviv University\PhD\Slutsky\Manuscripts\MCU\Results';
xlsName = 'mcu_suppTbl.xlsx';

% Load MEA ----------------------------------------------------------------
presets = {'steadyState'};
tblMea = mcu_tblMea('presets', presets, 'flgOtl', true);
tblTrans = tbl_trans(tblMea, 'varsInc', {'pBspk', 'ss_pBspk'}, 'logBase', 'logit');
tblMea.pBspk_trans = tblTrans.pBspk;

% Load In Vivo ------------------------------------------------------------
basepaths = [mcu_basepaths('wt'), mcu_basepaths('mcu')];
presets = {'brst'};
tblVivo = mcu_tblVivo('basepaths', basepaths, 'presets', presets, 'flgClean', true);


%% ========================================================================
% FIGURE 2
% =========================================================================


% A -----------------------------------------------------------------------

% FR
frml = 'fr ~ Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblMea, frml, 'dist', 'log-normal');
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lme_save('2A (FR)', lmeTbls, 'pathName', pathName, 'xlsName', xlsName)

% P_burst
frml = 'pBspk ~ Group * fr + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblMea, frml, 'dist', 'logit-normal');
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lme_save('2A (P_burst)', lmeTbls, 'pathName', pathName, 'xlsName', xlsName)


% B -----------------------------------------------------------------------

% FR
tblLme = tblVivo(tblVivo.Day == 'BSL', :);
frml = 'fr ~ Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLme, frml, 'dist', 'gamma');
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lme_save('2B (FR)', lmeTbls, 'pathName', pathName, 'xlsName', xlsName)

% P_burst
frml = 'pBspk ~ Group * fr + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLme, frml, 'dist', 'logit-normal');
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lme_save('2B (P_burst)', lmeTbls, 'pathName', pathName, 'xlsName', xlsName)


% C -----------------------------------------------------------------------
basepaths = [mcu_basepaths('wt_bsl_ripp'), mcu_basepaths('mcu_bsl')];
tblRipp = mcu_tblVivo('basepaths', basepaths, 'presets', {'ripp'});
frml = 'amp ~ Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblRipp, frml, 'dist', 'log-normal');
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lme_save('2C', lmeTbls, 'pathName', pathName, 'xlsName', xlsName)


% D -----------------------------------------------------------------------
frml = 'freq ~ Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblRipp, frml, 'dist', 'normal');
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lme_save('2D', lmeTbls, 'pathName', pathName, 'xlsName', xlsName)


% E -----------------------------------------------------------------------
frml = 'dur ~ Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblRipp, frml, 'dist', 'log-normal');
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lme_save('2E', lmeTbls, 'pathName', pathName, 'xlsName', xlsName)


% G -----------------------------------------------------------------------
tblRipp = mcu_tblVivo('basepaths', basepaths, 'presets', {'rippSpks', 'brst'}, 'flgClean', true);
frml = 'com ~ (fr + pBspk) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblRipp, frml, 'dist', 'normal', 'flgStnd', false);
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lme_save('2G', lmeTbls, 'pathName', pathName, 'xlsName', xlsName)

hFig = figure;
hAx = nexttile; lme_lsmeans(lmeMdl, {'pBspk', 'Group'}, 'transParams', lmeInfo.transParams, ...
    'hAx', hAx, 'xLims', {[0, 1], []});
hAx = nexttile; lme_lsmeans(lmeMdl, {'fr', 'Group'}, 'transParams', lmeInfo.transParams, ...
    'hAx', hAx); 


%% ========================================================================
% FIGURE 3
% =========================================================================

% G -----------------------------------------------------------------------
frml = 'fr ~ Group * Day + (Day|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblVivo, frml, 'dist', 'gamma');
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lme_save('3G', lmeTbls, 'pathName', pathName, 'xlsName', xlsName)

% J -----------------------------------------------------------------------

% FR
tblLme = tblVivo(tblVivo.Day == 'BAC3', :);
frml = 'fr ~ Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLme, frml, 'dist', 'gamma');
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lme_save('3J (FR)', lmeTbls, 'pathName', pathName, 'xlsName', xlsName)

% P_burst
frml = 'pBspk ~ Group * fr + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLme, frml, 'dist', 'logit-normal');
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lme_save('3J (P_burst)', lmeTbls, 'pathName', pathName, 'xlsName', xlsName)


%% ========================================================================
% FIGURE 4
% =========================================================================

% D -----------------------------------------------------------------------
tblMea.bGain = log((tblMea.ss_frBspk) ./ (tblMea.frBspk));
tblMea.sGain = log((tblMea.ss_frSspk) ./ (tblMea.frSspk));
frml = 'sGain ~ bGain * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblMea, frml, ...
    'dist', 'normal', 'flgStnd', false);
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lme_save('4D', lmeTbls, 'pathName', pathName, 'xlsName', xlsName)

% F -----------------------------------------------------------------------
frml = 'ss_frSspk ~ pBspk + fr + (1|Name)';
xVar = 'pBspk';
mVar = 'ss_frBspk';
distM = 'log-normal';
distY = distM;
transX = [];

% Per Group (standardized)
tblWt = tblMea(tblMea.Group == 'Control', :);
[~, tmpl] = tbl_trans(tblWt, 'varsInc', {'fr', 'pBspk', 'ss_frBspk'}, 'logBase', 10, 'skewThr', 2, 'flgZ', false);
tmpl.varsTrans.pBspk.logBase = transX;
tmpl.varsTrans.ss_frBspk.logBase = 'e'; % Force ln to match Path A response scaling
resWt = lme_mediation(tblWt, frml, xVar, mVar, 'distM', distM, 'distY', distY, 'transTemplate', tmpl);

tblMcu = tblMea(tblMea.Group == 'MCU-KO', :);
[~, tmpl] = tbl_trans(tblMcu, 'varsInc', {'fr', 'pBspk', 'ss_frBspk'}, 'logBase', 10, 'skewThr', 2, 'flgZ', false);
resMcu = lme_mediation(tblMcu, frml, xVar, mVar, 'distM', distM, 'distY', distY, 'transTemplate', tmpl);

% Consolidate Mediation Summaries into one sheet
medTbls = struct('Title', {}, 'Table', {});
medTbls(1).Title = 'Mediation Summary: Control';
medTbls(1).Table = resWt.paths;
medTbls(2).Title = 'Mediation Summary: MCU-KO';
medTbls(2).Table = resMcu.paths;
lme_save('Mediation', medTbls, 'pathName', pathName, 'xlsName', xlsName)

resWt.plot.X = tblWt.pBspk_trans;
resWt.plot = tbl_trans(resWt.plot, 'varsInc', {'M'}, 'logBase', 10, 'skewThr', 2, 'flgZ', false);
lme_mediationPlot(resWt)
resMcu.plot.X = tblMcu.pBspk_trans;
resMcu.plot = tbl_trans(resMcu.plot, 'varsInc', {'M'}, 'logBase', 10, 'skewThr', 2, 'flgZ', false);
lme_mediationPlot(resMcu)


% Combined Models (unstandardized)
[~, tmpl] = tbl_trans(tblMea, 'varsInc', {'fr', 'pBspk', 'ss_frBspk'}, 'logBase', 10, 'skewThr', 2, 'flgZ', false);
tmpl.varsTrans.pBspk.logBase = transX;
tmpl.varsTrans.ss_frBspk.logBase = 'e'; % Force ln for consistency with mediation

% X -> M 
frml = 'ss_frBspk ~ (fr + pBspk) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblMea, frml, 'dist', 'log-normal', 'flgStnd', false, 'transTemplate', tmpl);
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lme_save('4F (X->M)', lmeTbls, 'pathName', pathName, 'xlsName', xlsName)

% X -> Y
frml = 'ss_frSspk ~ (fr + pBspk) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblMea, frml, 'dist', 'log-normal', 'flgStnd', false, 'transTemplate', tmpl);
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lme_save('4F (X->Y)', lmeTbls, 'pathName', pathName, 'xlsName', xlsName)

% X -> Y | M 
frml = 'ss_frSspk ~ (fr + pBspk + ss_frBspk) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblMea, frml, 'dist', 'log-normal', 'flgStnd', false, 'transTemplate', tmpl);
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lme_save('4F (X+M->Y)', lmeTbls, 'pathName', pathName, 'xlsName', xlsName)

%% ========================================================================
%  ABLATION
%  ========================================================================

frml = 'ss_fr ~ (frBspk + frSspk)';
partMode = 'split';
tblWt = tbl(tbl.Group == 'Control', :);
tblMcu = tbl(tbl.Group == 'MCU-KO', :);

abl = lme_ablation(tblWt, frml, 'dist', 'log-normal', 'partitionMode', partMode, 'nrep', 10);
ablTbl.Title = 'ABLATION SUMMARY (Pooled R2_OOS)';
ablTbl.Table = table(abl.vars', round(abl.pR2', 4), 'VariableNames', {'Ablated_Feature', 'Ctrl'});
abl = lme_ablation(tblMcu, frml, 'dist', 'log-normal', 'partitionMode', partMode, 'nrep', 10);
ablTbl.Table{:, 'KO'} = round(abl.pR2', 4);

headWt.Title = 'FULL MODEL: CONTROL';
headWt.Table = [];
headMcu.Title = 'FULL MODEL: MCU-KO';
headMcu.Table = [];

[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblWt, frml, 'dist', 'log-normal', 'flgStnd', false);
tblWt = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblMcu, frml, 'dist', 'log-normal', 'flgStnd', false);
tblMcu = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lmeTbls = [ablTbl, headWt, tblWt, headMcu, tblMcu];
lme_save('Ablation', lmeTbls, 'pathName', pathName, 'xlsName', xlsName)


%% ========================================================================
% FIGURE S3
% =========================================================================

% C -----------------------------------------------------------------------
frml = 'sGain ~ (pBspk + fr + bGain) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tblMea, frml, ...
    'dist', 'normal', 'flgStnd', false);
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lme_save('S3C', lmeTbls, 'pathName', pathName, 'xlsName', xlsName)


%% ========================================================================
% FIGURE S4
% =========================================================================

% A,B ---------------------------------------------------------------------
frml = 'ss_frBspk ~ (fr + pBspk) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblMea, frml, 'dist', 'log-normal', 'flgStnd', false);
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lme_save('S4A,B (fr_burst)', lmeTbls, 'pathName', pathName, 'xlsName', xlsName)

frml = 'ss_frSspk ~ (fr + pBspk) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblMea, frml, 'dist', 'log-normal', 'flgStnd', false);
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lme_save('S4A,B (fr_single)', lmeTbls, 'pathName', pathName, 'xlsName', xlsName)

frml = 'ss_fr ~ (fr + pBspk) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblMea, frml, 'dist', 'log-normal', 'flgStnd', false);
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lme_save('S4A,B (fr)', lmeTbls, 'pathName', pathName, 'xlsName', xlsName)


% EOF