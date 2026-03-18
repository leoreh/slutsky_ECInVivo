
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

flgPlot = false;

%% ========================================================================
% Table S1
% =========================================================================
tblIdx = 1;
sheetNames{tblIdx} = ['S' num2str(tblIdx)];
tblInfo{tblIdx} = 'BSL Firing';
dataSet{tblIdx} = 'MEA';
tblPnls{tblIdx} = '1K-N';

frml = 'fr ~ Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblMea, frml, 'dist', 'log-normal');
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);

frml = 'bRate ~ Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblMea, frml, 'dist', 'log-normal');
lmeTbls = [lmeTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];

frml = 'nBspk ~ Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblMea, frml);
lmeTbls = [lmeTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];

frml = 'pBspk ~ Group * fr + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblMea, frml, 'dist', 'logit-normal');
lmeTbls = [lmeTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];
lme_save(sheetNames{1}, lmeTbls, 'pathName', pathName, 'xlsName', xlsName)

if flgPlot
    tblGUI_bar(tblMea, 'yVar', 'pBspk', 'xVar', 'Group');
    tblGUI_scatHist(tblMea, 'xVar', 'fr', 'yVar', 'bRate', 'grpVar', 'Group');
end


%% ========================================================================
% Table S2
% =========================================================================
tblIdx = 2;
sheetNames{tblIdx} = ['S' num2str(tblIdx)];
tblInfo{tblIdx} = 'BSL Firing';
dataSet{tblIdx} = 'In Vivo';
tblPnls{tblIdx} = '2B';

tblLme = tblVivo(tblVivo.Day == 'BSL', :);
frml = 'fr ~ Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLme, frml, 'dist', 'gamma');
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);

frml = 'pBspk ~ Group * fr + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLme, frml, 'dist', 'logit-normal');
lmeTbls = [lmeTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];
lme_save(sheetNames{2}, lmeTbls, 'pathName', pathName, 'xlsName', xlsName)

frml = 'bFreq ~ Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLme, frml);

frml = 'nBspk ~ Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLme, frml);

if flgPlot
    tblGUI_bar(tblLme, 'yVar', 'pBspk', 'xVar', 'Group');
    tblGUI_scatHist(tblLme, 'xVar', 'fr', 'yVar', 'bRate', 'grpVar', 'Group');
end

%% ========================================================================
% Table S3
% =========================================================================
tblIdx = 3;
sheetNames{tblIdx} = ['S' num2str(tblIdx)];
tblInfo{tblIdx} = 'SWR Properties';
dataSet{tblIdx} = 'In Vivo';
tblPnls{tblIdx} = '2C-E';

basepaths = [mcu_basepaths('wt_bsl_ripp'), mcu_basepaths('mcu_bsl')];
tblRipp = mcu_tblVivo('basepaths', basepaths, 'presets', {'ripp'});

frml = 'amp ~ Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblRipp, frml, 'dist', 'log-normal');
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);

frml = 'freq ~ Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblRipp, frml, 'dist', 'normal');
lmeTbls = [lmeTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];

frml = 'dur ~ Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblRipp, frml, 'dist', 'log-normal');
lmeTbls = [lmeTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];
lme_save(sheetNames{3}, lmeTbls, 'pathName', pathName, 'xlsName', xlsName)

%% ========================================================================
% Table S4
% =========================================================================
tblIdx = 4;
sheetNames{tblIdx} = ['S' num2str(tblIdx)];
tblInfo{tblIdx} = 'Spike COM during SWR';
dataSet{tblIdx} = 'In Vivo';
tblPnls{tblIdx} = '2G,H';

basepaths = [mcu_basepaths('wt_bsl_ripp'), mcu_basepaths('mcu_bsl')];
tblRipp = mcu_tblVivo('basepaths', basepaths, 'presets', {'rippSpks', 'brst'}, 'flgClean', true);

frml = 'com ~ (fr + pBspk) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblRipp, frml, 'dist', 'normal', 'flgStnd', false);
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lme_save(sheetNames{4}, lmeTbls, 'pathName', pathName, 'xlsName', xlsName)

if flgPlot
    hFig = figure;
    hAx = nexttile; lme_lsmeans(lmeMdl, {'pBspk', 'Group'}, 'transParams', lmeInfo.transParams, ...
        'hAx', hAx, 'xLims', {[0, 1], []});
    hAx = nexttile; lme_lsmeans(lmeMdl, {'fr', 'Group'}, 'transParams', lmeInfo.transParams, ...
        'hAx', hAx);
end

%% ========================================================================
% Table S5
% =========================================================================
tblIdx = 5;
sheetNames{tblIdx} = ['S' num2str(tblIdx)];
tblInfo{tblIdx} = 'FR during BAC';
dataSet{tblIdx} = 'In Vivo';
tblPnls{tblIdx} = '3H';

frml = 'fr ~ Group * Day + (Day|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblVivo, frml, 'dist', 'gamma');
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lme_save(sheetNames{tblIdx}, lmeTbls, 'pathName', pathName, 'xlsName', xlsName)

%% ========================================================================
% Table S6
% =========================================================================
tblIdx = 6;
sheetNames{tblIdx} = ['S' num2str(tblIdx)];
tblInfo{tblIdx} = 'SS Firing';
dataSet{tblIdx} = 'In Vivo';
tblPnls{tblIdx} = '3J';

tblLme = tblVivo(tblVivo.Day == 'BAC3', :);
frml = 'fr ~ Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLme, frml, 'dist', 'gamma');
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);

frml = 'pBspk ~ Group * fr + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLme, frml, 'dist', 'logit-normal');
lmeTbls = [lmeTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];
lme_save(sheetNames{tblIdx}, lmeTbls, 'pathName', pathName, 'xlsName', xlsName)


%% ========================================================================
% Table S7
% =========================================================================
tblIdx = 7;
sheetNames{tblIdx} = ['S' num2str(tblIdx)];
tblInfo{tblIdx} = 'Firing Gain during FRH';
dataSet{tblIdx} = 'MEA';
tblPnls{tblIdx} = '4D, S3C';

tblMea.bGain = log((tblMea.ss_frBspk) ./ (tblMea.frBspk));
tblMea.sGain = log((tblMea.ss_frSspk) ./ (tblMea.frSspk));
frml = 'sGain ~ bGain * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblMea, frml, 'dist', 'normal', 'flgStnd', false);
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);

frml = 'sGain ~ (pBspk + fr + bGain) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tblMea, frml, 'dist', 'normal', 'flgStnd', false);
lmeTbls = [lmeTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];
lme_save(sheetNames{tblIdx}, lmeTbls, 'pathName', pathName, 'xlsName', xlsName)

%% ========================================================================
% Table S8
% =========================================================================
tblIdx = 8;
sheetNames{tblIdx} = ['S' num2str(tblIdx)];
tblInfo{tblIdx} = 'Mediation Analysis';
dataSet{tblIdx} = 'MEA';
tblPnls{tblIdx} = '4F, S4';

frml = 'ss_frSspk ~ pBspk + fr + (1|Name)';
xVar = 'pBspk';
mVar = 'ss_frBspk';
distM = 'log-normal';
distY = distM;
transX = [];

% Per Group 
tblWt = tblMea(tblMea.Group == 'Control', :);
[~, tmpl] = tbl_trans(tblWt, 'varsInc', {'fr', 'pBspk', 'ss_frBspk'}, 'logBase', 10, 'skewThr', 2, 'flgZ', false);
tmpl.varsTrans.pBspk.logBase = transX;
tmpl.varsTrans.ss_frBspk.logBase = 'e'; % Force ln to match Path A response scaling
resWt = lme_mediation(tblWt, frml, xVar, mVar, 'distM', distM, 'distY', distY, 'transTemplate', tmpl);

tblMcu = tblMea(tblMea.Group == 'MCU-KO', :);
[~, tmpl] = tbl_trans(tblMcu, 'varsInc', {'fr', 'pBspk', 'ss_frBspk'}, 'logBase', 10, 'skewThr', 2, 'flgZ', false);
resMcu = lme_mediation(tblMcu, frml, xVar, mVar, 'distM', distM, 'distY', distY, 'transTemplate', tmpl);

% Consolidate Summaries
medTbls = struct('Title', {}, 'Table', {});
medTbls(1).Title = 'MEDIATION SUMMARY: CONTROL';
medTbls(1).Table = resWt.paths;
medTbls(2).Title = 'MEDIATION SUMMARY: MCU-KO';
medTbls(2).Table = resMcu.paths;

% Plot
if flgPlot
    resWt.plot.X = tblWt.pBspk_trans;
    resWt.plot = tbl_trans(resWt.plot, 'varsInc', {'M'}, 'logBase', 10, 'skewThr', 2, 'flgZ', false);
    lme_mediationPlot(resWt)
    resMcu.plot.X = tblMcu.pBspk_trans;
    resMcu.plot = tbl_trans(resMcu.plot, 'varsInc', {'M'}, 'logBase', 10, 'skewThr', 2, 'flgZ', false);
    lme_mediationPlot(resMcu)
end

% Combined Models
[~, tmpl] = tbl_trans(tblMea, 'varsInc', {'fr', 'pBspk', 'ss_frBspk'}, 'logBase', 10, 'skewThr', 2, 'flgZ', false);
tmpl.varsTrans.pBspk.logBase = transX;
tmpl.varsTrans.ss_frBspk.logBase = 'e'; % Force ln for consistency with mediation

frml = 'ss_frBspk ~ (fr + pBspk) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblMea, frml, 'dist', 'log-normal', 'flgStnd', false, 'transTemplate', tmpl);
lmeTbls = [medTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];

frml = 'ss_frSspk ~ (fr + pBspk) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblMea, frml, 'dist', 'log-normal', 'flgStnd', false, 'transTemplate', tmpl);
lmeTbls = [lmeTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];

frml = 'ss_frSspk ~ (fr + pBspk + ss_frBspk) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblMea, frml, 'dist', 'log-normal', 'flgStnd', false, 'transTemplate', tmpl);
lmeTbls = [lmeTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];

frml = 'ss_fr ~ (fr + pBspk) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblMea, frml, 'dist', 'log-normal', 'flgStnd', false);
lmeTbls = [lmeTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];

lme_save(sheetNames{tblIdx}, lmeTbls, 'pathName', pathName, 'xlsName', xlsName)


%% ========================================================================
% Table S9
% =========================================================================
tblIdx = 9;
sheetNames{tblIdx} = ['S' num2str(tblIdx)];
tblInfo{tblIdx} = 'Feature Abalation';
dataSet{tblIdx} = 'MEA';
tblPnls{tblIdx} = '4G';

frml = 'ss_fr ~ (frBspk + frSspk)';
partMode = 'split';
tblWt = tblMea(tblMea.Group == 'Control', :);
tblMcu = tblMea(tblMea.Group == 'MCU-KO', :);

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

lme_save(sheetNames{tblIdx}, lmeTbls, 'pathName', pathName, 'xlsName', xlsName)


%% ========================================================================
% TOC
% =========================================================================

nTbls = length(sheetNames);
indexData = cell(nTbls + 1, 3);
indexData(1, :) = {'Table', 'Description', 'Figure Panels'};

for iTbl = 1:nTbls
    linkFrml = sprintf('=HYPERLINK("#''%s''!A1", "%s")', sheetNames{iTbl}, sheetNames{iTbl});
    indexData{iTbl + 1, 1} = linkFrml;
    indexData{iTbl + 1, 2} = tblInfo{iTbl};
    indexData{iTbl + 1, 3} = dataSet{iTbl};
    indexData{iTbl + 1, 4} = tblPnls{iTbl};
end

writecell(indexData, fullfile(pathName, xlsName), 'Sheet', 'TOC');
mcu_xlsFormat(xlsName, 'pathName', pathName)
system('taskkill /f /im excel.exe')
