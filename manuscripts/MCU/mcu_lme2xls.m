
% File Names
pathName = 'D:\OneDrive - Tel-Aviv University\PhD\Slutsky\Manuscripts\MCU\Results';
xlsName = 'mcu_suppTbl.xlsx';

% Load MEA ----------------------------------------------------------------
presets = {'steadyState'};
tblMea = mcu_tblMea('presets', presets, 'flgOtl', true);
tblTrans = tbl_trans(tblMea, 'varsInc', {'pBurst', 'ss_pBurst'}, 'logBase', 'logit');
tblMea.pBurst_trans = tblTrans.pBurst;

% Load In Vivo ------------------------------------------------------------
presets = {'burst'};
tblVivo = mcu_tblVivo('presets', presets, 'flgClean', true);
tblTrans = tbl_trans(tblVivo, 'varsInc', {'pBurst'}, 'logBase', 'logit');
tblVivo.pBurst_trans = tblTrans.pBurst;

% Load Imaging -----------------------------------------------------------
spDir = fileparts(which('spontCa_detect'));
load(fullfile(spDir, 'cache', 'spontCa_tbl.mat'), 'tblCell', 'tblEvent', 'fs');
[tblCell, tblEvent] = spontCa2_metrics(tblCell, tblEvent, fs, 'aggFcn', 'mean');

flgPlot = false;

%% ========================================================================
% Table S1
% =========================================================================
tblIdx = 1;
sheetNames{tblIdx} = ['S' num2str(tblIdx)];
tblInfo{tblIdx} = 'Spontaneous Ca Transients';
dataSet{tblIdx} = 'Imaging';
tblPnls{tblIdx} = '1E-F; S1B-C';

% Fig. 1E - per-event transfer function (paired cyto events)
tblLme = tblEvent(tblEvent.compartment == 'Cyto' & tblEvent.paired, :);
frml = 'pairAmp ~ amp * genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLme, frml, 'dist', 'log-normal', 'flgStnd', false);
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);

% Fig. 1F - cell-level amplitude by compartment
frml = 'amp ~ compartment * genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblCell, frml, 'dist', 'log-normal', 'flgStnd', false);
lmeTbls = [lmeTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];

% Fig. S1B - cell-level event rate by compartment
frml = 'rate ~ compartment * genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblCell, frml, 'dist', 'log-normal', 'flgStnd', false);
lmeTbls = [lmeTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];

% Fig. S1C - per-event mito coupling probability (binomial GLMM)
tblLme = tblEvent(tblEvent.compartment == 'Mito', :);
tblLme.paired = double(tblLme.paired);
frml = 'paired ~ genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLme, frml, 'dist', 'binomial', 'flgStnd', false);
lmeTbls = [lmeTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];

lme_save(sheetNames{tblIdx}, lmeTbls, 'pathName', pathName, 'xlsName', xlsName, ...
    'tblInfo', tblInfo{tblIdx}, 'dataSet', dataSet{tblIdx}, 'tblPnls', tblPnls{tblIdx})

if flgPlot
    tblGUI_bar(tblCell, 'yVar', 'amp', 'xVar', 'compartment', 'grpVar', 'genotype');
    tblGUI_bar(tblCell, 'yVar', 'rate', 'xVar', 'compartment', 'grpVar', 'genotype');
    tblGUI_scatHist(tblEvent(tblEvent.compartment == 'Cyto' & tblEvent.paired, :), ...
        'xVar', 'amp', 'yVar', 'pairAmp', 'grpVar', 'genotype');
end


%% ========================================================================
% Table S2
% =========================================================================
tblIdx = 2;
sheetNames{tblIdx} = ['S' num2str(tblIdx)];
tblInfo{tblIdx} = 'BSL Firing';
dataSet{tblIdx} = 'MEA';
tblPnls{tblIdx} = '1J-M';

frml = 'fr ~ genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblMea, frml, 'dist', 'log-normal');
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);

frml = 'br ~ genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblMea, frml, 'dist', 'log-normal');
lmeTbls = [lmeTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];

frml = 'bSize ~ genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblMea, frml, 'dist', 'log-normal');
lmeTbls = [lmeTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];

frml = 'pBurst ~ genotype * fr + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblMea, frml, 'dist', 'logit-normal');
lmeTbls = [lmeTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];
lme_save(sheetNames{tblIdx}, lmeTbls, 'pathName', pathName, 'xlsName', xlsName, ...
    'tblInfo', tblInfo{tblIdx}, 'dataSet', dataSet{tblIdx}, 'tblPnls', tblPnls{tblIdx})

if flgPlot
    tblGUI_bar(tblMea, 'yVar', 'pBurst', 'xVar', 'genotype');
    tblGUI_scatHist(tblMea, 'xVar', 'fr', 'yVar', 'br', 'grpVar', 'genotype');
end


%% ========================================================================
% Table S3
% =========================================================================
tblIdx = 3;
sheetNames{tblIdx} = ['S' num2str(tblIdx)];
tblInfo{tblIdx} = 'BSL Firing';
dataSet{tblIdx} = 'In Vivo';
tblPnls{tblIdx} = '2D-E; S2F-G';

tblLme = tblVivo(tblVivo.day == 'BSL', :);
frml = 'fr ~ genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLme, frml, 'dist', 'gamma');
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);

frml = 'pBurst ~ genotype * fr + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLme, frml, 'dist', 'logit-normal');
lmeTbls = [lmeTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];

frml = 'br ~ genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLme, frml, 'dist', 'gamma');
lmeTbls = [lmeTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];

frml = 'bSize ~ genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLme, frml, 'dist', 'log-normal');
lmeTbls = [lmeTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];
lme_save(sheetNames{tblIdx}, lmeTbls, 'pathName', pathName, 'xlsName', xlsName, ...
    'tblInfo', tblInfo{tblIdx}, 'dataSet', dataSet{tblIdx}, 'tblPnls', tblPnls{tblIdx})

if flgPlot
    tblGUI_bar(tblLme, 'yVar', 'pBurst', 'xVar', 'genotype');
    tblGUI_scatHist(tblLme, 'xVar', 'fr', 'yVar', 'br', 'grpVar', 'genotype');
end


%% ========================================================================
% Table S4
% =========================================================================
tblIdx = 4;
sheetNames{tblIdx} = ['S' num2str(tblIdx)];
tblInfo{tblIdx} = 'Spike COM during SWR';
dataSet{tblIdx} = 'In Vivo';
tblPnls{tblIdx} = '2F; S2G';

basepaths = [mcu_basepaths('wt_bsl_ripp'), mcu_basepaths('mcu_bsl')];
[tblRipp, ~, ~, xVec] = mcu_tblVivo('basepaths', basepaths, 'presets', {'rippSpks', 'burst'}, 'flgClean', true);
tblTrans = tbl_trans(tblRipp, 'varsInc', {'pBurst'}, 'logBase', 'logit');
tblRipp.pBurst_trans = tblTrans.pBurst;

frml = 'com ~ genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblRipp, frml, 'dist', 'normal', 'flgStnd', false);
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);

frml = 'com ~ (fr + pBurst) + genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblRipp, frml, 'dist', 'normal', 'flgStnd', false);
lmeTbls = [lmeTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];
lme_save(sheetNames{tblIdx}, lmeTbls, 'pathName', pathName, 'xlsName', xlsName, ...
    'tblInfo', tblInfo{tblIdx}, 'dataSet', dataSet{tblIdx}, 'tblPnls', tblPnls{tblIdx})

if flgPlot
    hFig = figure;
    hAx = nexttile; pdRes = lme_lsmeans(lmeMdl, {'pBurst', 'genotype'}, 'transParams', lmeInfo.transParams, ...
        'hAx', hAx, 'xLims', {[0, 1], []});
    hAx = nexttile; pdRes = lme_lsmeans(lmeMdl, {'fr', 'genotype'}, 'transParams', lmeInfo.transParams, ...
        'hAx', hAx); 
    tblGUI_xy(xVec, tblRipp, 'grpVar', 'genotype');
    tblGUI_scatHist(tblRipp, 'grpVar', 'genotype');
end


%% ========================================================================
% Table S5
% =========================================================================
tblIdx = 5;
sheetNames{tblIdx} = ['S' num2str(tblIdx)];
tblInfo{tblIdx} = 'SWR Properties';
dataSet{tblIdx} = 'In Vivo';
tblPnls{tblIdx} = '2H; S2H';

basepaths = [mcu_basepaths('wt_bsl_ripp'), mcu_basepaths('mcu_bsl')];
tblRipp = mcu_tblVivo('basepaths', basepaths, 'presets', {'ripp'});

frml = 'amp ~ genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblRipp, frml, 'dist', 'log-normal');
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);

frml = 'freq ~ genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblRipp, frml, 'dist', 'normal');
lmeTbls = [lmeTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];

frml = 'dur ~ genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblRipp, frml, 'dist', 'log-normal');
lmeTbls = [lmeTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];
lme_save(sheetNames{tblIdx}, lmeTbls, 'pathName', pathName, 'xlsName', xlsName, ...
    'tblInfo', tblInfo{tblIdx}, 'dataSet', dataSet{tblIdx}, 'tblPnls', tblPnls{tblIdx})

%% ========================================================================
% Table S6
% =========================================================================
tblIdx = 6;
sheetNames{tblIdx} = ['S' num2str(tblIdx)];
tblInfo{tblIdx} = 'FRH during BAC';
dataSet{tblIdx} = 'In Vivo';
tblPnls{tblIdx} = '4C';

tblLme = tblVivo; tblLme(tblLme.sbjID == 'lh137', :) = [];
frml = 'fr ~ genotype * day + (day|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblVivo, frml, 'dist', 'gamma');
lmeStats = lme_postHoc(lmeMdl, 'contrasts', [1 : 9, 12, 15, 17 : 19]);
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);

lme_save(sheetNames{tblIdx}, lmeTbls, 'pathName', pathName, 'xlsName', xlsName, ...
    'tblInfo', tblInfo{tblIdx}, 'dataSet', dataSet{tblIdx}, 'tblPnls', tblPnls{tblIdx})


%% ========================================================================
% Table S7
% =========================================================================
tblIdx = 7;
sheetNames{tblIdx} = ['S' num2str(tblIdx)];
tblInfo{tblIdx} = 'SS Burstiness';
dataSet{tblIdx} = 'In Vivo';
tblPnls{tblIdx} = 'S4D';

tblLme = tblVivo(tblVivo.day == 'BAC3', :);
frml = 'pBurst ~ genotype * fr + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLme, frml, 'dist', 'logit-normal');
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lme_save(sheetNames{tblIdx}, lmeTbls, 'pathName', pathName, 'xlsName', xlsName, ...
    'tblInfo', tblInfo{tblIdx}, 'dataSet', dataSet{tblIdx}, 'tblPnls', tblPnls{tblIdx})

if flgPlot
    tblGUI_bar(tblLme, 'yVar', 'pBurst', 'xVar', 'genotype');
end

%% ========================================================================
% Table S8
% =========================================================================
tblIdx = 8;
sheetNames{tblIdx} = ['S' num2str(tblIdx)];
tblInfo{tblIdx}    = 'FR Component x Time during FRH';
dataSet{tblIdx}    = 'MEA';
tblPnls{tblIdx}    = 'S5C';

% Reshape tblMea to long format: one row per (unit × component × time).
% Each unit contributes 4 rows crossing:
%   component : {'bSpk', 'sSpk'} — burst vs. single spike firing rate
%   time     : {'BSL',  'SS'}   — baseline vs. steady state
tblLong = stack(tblMea, {'frBurst', 'frSingle', 'ss_frBurst', 'ss_frSingle'}, ...
    'NewDataVariableName', 'fr', ...
    'IndexVariableName', 'SourceVar', ...
    'ConstantVariables', {'genotype', 'sbjID', 'unitID'});
tblLong.time = categorical(tblLong.SourceVar, ...
    {'frBurst', 'frSingle', 'ss_frBurst', 'ss_frSingle'}, ...
    {'BSL',    'BSL',    'SS',        'SS'});
tblLong.component = categorical(tblLong.SourceVar, ...
    {'frBurst', 'frSingle', 'ss_frBurst', 'ss_frSingle'}, ...
    {'bSpk',   'sSpk',   'bSpk',      'sSpk'});
tblLong.SourceVar = [];

frml = 'fr ~ component * time * genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblLong, frml, 'dist', 'log-normal');
lmeStats = lme_postHoc(lmeMdl, 'contrasts', [1 : 9, 32 : 39]);
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);
lme_save(sheetNames{tblIdx}, lmeTbls, 'pathName', pathName, 'xlsName', xlsName, ...
    'tblInfo', tblInfo{tblIdx}, 'dataSet', dataSet{tblIdx}, 'tblPnls', tblPnls{tblIdx})

%% ========================================================================
% Table S9
% =========================================================================
tblIdx = 9;
sheetNames{tblIdx} = ['S' num2str(tblIdx)];
tblInfo{tblIdx} = 'Firing Gain during FRH';
dataSet{tblIdx} = 'MEA';
tblPnls{tblIdx} = '4H,K; S5D';

tblMea.bGain = log((tblMea.ss_frBurst) ./ (tblMea.frBurst));
tblMea.sGain = log((tblMea.ss_frSingle) ./ (tblMea.frSingle));
tblMea.frGain = log((tblMea.ss_fr) ./ (tblMea.fr));

% Model 1: Proportional allocation (Figure 4H)
frml = 'sGain ~ bGain * genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblMea, frml, 'dist', 'normal', 'flgStnd', false);
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);

% Model 2: Total FR gain vs burstiness (Figure 4K)
frml = 'frGain ~ (pBurst + fr) * genotype  + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tblMea, frml, 'dist', 'normal', 'flgStnd', false);
lmeTbls = [lmeTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];
lme_save(sheetNames{tblIdx}, lmeTbls, 'pathName', pathName, 'xlsName', xlsName, ...
    'tblInfo', tblInfo{tblIdx}, 'dataSet', dataSet{tblIdx}, 'tblPnls', tblPnls{tblIdx})

% Partial regression for proportional allocation, controlling for baseline
% burstiness. This model tests whether the allocation coefficient (beta)
% differs between genotypes after accounting for baseline firing pattern.
% Serves as the statistical basis for the added-variable plot (Figure S5D).
frml = 'sGain ~ (pBurst + bGain) * genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblMea, frml, ...
    'dist', 'normal', 'flgStnd', false);
lmeTbls = [lmeTbls, lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo)];
lme_save(sheetNames{tblIdx}, lmeTbls, 'pathName', pathName, 'xlsName', xlsName, ...
    'tblInfo', tblInfo{tblIdx}, 'dataSet', dataSet{tblIdx}, 'tblPnls', tblPnls{tblIdx})

if flgPlot
    tblGUI_scatHist(tblMea, 'grpVar', 'genotype');
end




%% ========================================================================
% Table S10
% =========================================================================
tblIdx = 10;
sheetNames{tblIdx} = ['S' num2str(tblIdx)];
tblInfo{tblIdx} = 'Feature Ablation';
dataSet{tblIdx} = 'MEA';
tblPnls{tblIdx} = '4I';

frml = 'ss_fr ~ (frBurst + frSingle)';
partMode = 'split';
tblWt = tblMea(tblMea.genotype == 'Control', :);
tblMcu = tblMea(tblMea.genotype == 'MCU-KO', :);

abl = lme_ablation(tblWt, frml, 'dist', 'log-normal', 'partitionMode', partMode, 'nrep', 10, 'flgPlot', flgPlot);
ablTbl.Title = 'ABLATION SUMMARY (Pooled R2_OOS)';
ablTbl.Table = table(abl.vars', round(abl.pR2', 4), 'VariableNames', {'Ablated_Feature', 'Ctrl'});
abl = lme_ablation(tblMcu, frml, 'dist', 'log-normal', 'partitionMode', partMode, 'nrep', 10, 'flgPlot', flgPlot);
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

lme_save(sheetNames{tblIdx}, lmeTbls, 'pathName', pathName, 'xlsName', xlsName, ...
    'tblInfo', tblInfo{tblIdx}, 'dataSet', dataSet{tblIdx}, 'tblPnls', tblPnls{tblIdx})


%% ========================================================================
% TOC
% =========================================================================

nTbls = length(sheetNames);
indexData = cell(nTbls + 1, 4);
indexData(1, :) = {'Table', 'Description', 'Data Set', 'Figure Panels'};

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
