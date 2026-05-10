
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

flgPlot = false;

%% ========================================================================
% Table S1
% =========================================================================
tblIdx = 1;
sheetNames{tblIdx} = ['S' num2str(tblIdx)];
tblInfo{tblIdx} = 'BSL Firing';
dataSet{tblIdx} = 'MEA';
tblPnls{tblIdx} = '1K-N';

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
% Table S2
% =========================================================================
tblIdx = 2;
sheetNames{tblIdx} = ['S' num2str(tblIdx)];
tblInfo{tblIdx} = 'BSL Firing';
dataSet{tblIdx} = 'In Vivo';
tblPnls{tblIdx} = '2D-E; S2E-F';

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
% Table S3
% =========================================================================
tblIdx = 3;
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
tblInfo{tblIdx} = 'FRH during BAC';
dataSet{tblIdx} = 'In Vivo';
tblPnls{tblIdx} = '3G';

tblLme = tblVivo; tblLme(tblLme.sbjID == 'lh137', :) = [];
frml = 'fr ~ genotype * day + (day|sbjID)';
[lmeMdl, lmeStats, lmeInfo] = lme_analyse(tblVivo, frml, 'dist', 'gamma');
lmeStats = lme_postHoc(lmeMdl, 'contrasts', [1 : 9, 12, 15, 17 : 19]);
lmeTbls = lme_mdl2tbls(lmeMdl, lmeStats, lmeInfo);

lme_save(sheetNames{tblIdx}, lmeTbls, 'pathName', pathName, 'xlsName', xlsName, ...
    'tblInfo', tblInfo{tblIdx}, 'dataSet', dataSet{tblIdx}, 'tblPnls', tblPnls{tblIdx})

% --- Distribution shape across days (per genotype) -----------------------
% Cross-sectional quantification of the RS firing rate distribution per
% day. Mean and median are reported in Hz for interpretability; variance
% and IQR are computed on log10(FR), where the distribution is
% approximately symmetric. Shape is compared against BSL with a two-
% sample KS test on log-FR. Units are sorted independently per 24-h
% session, so this describes the per-session population state, not
% trajectories of the same neurons across days. Kept in the workspace
% only (distRows) -- not written to the supp Excel file.

% Reload table to include washout, then remove bad and FS units, remove bac on, bac off
tblVivo = mcu_tblVivo('presets', presets, 'flgClean', false);
tblLme = tblVivo; tblLme(tblLme.sbjID == 'lh137', :) = [];
tblLme(tblLme.unitType == 'Other', :) = [];
tblLme(tblLme.unitType == 'FS', :) = [];
tblLme.unitType = removecats(tblLme.unitType, {'Other', 'FS'});
tblLme.unitType = [];
tblLme(tblLme.day == 'BAC_ON', :) = [];
tblLme(tblLme.day == 'BAC_OFF', :) = [];
tblLme.day = removecats(tblLme.day, {'BAC_ON', 'BAC_OFF'});

logFR = log10(tblVivo.fr); logFR(isinf(logFR)) = NaN;
grps = categories(removecats(tblVivo.genotype));
days = categories(removecats(tblVivo.day));
distRows = table;
for iGrp = 1:numel(grps)
    idxGrp = tblVivo.genotype == grps{iGrp};
    vBslLog = logFR(idxGrp & tblVivo.day == 'BSL' & ~isnan(logFR));
    for iDay = 1:numel(days)
        idxDay = idxGrp & tblVivo.day == days{iDay};
        vLog   = logFR(idxDay & ~isnan(logFR));
        vHz    = tblVivo.fr(idxDay);
        vHz    = vHz(~isnan(vHz) & vHz > 0);
        if isempty(vLog), continue; end
        if strcmp(days{iDay}, 'BSL')
            dKS = NaN; pKS = NaN;
        else
            [~, pKS, dKS] = kstest2(vBslLog, vLog);
        end
        distRows = [distRows; table(string(grps{iGrp}), string(days{iDay}), ...
            numel(vLog), mean(vHz), median(vHz), var(vLog), iqr(vLog), ...
            dKS, pKS, ...
            'VariableNames', {'Genotype','Day','n','Mean_Hz','Median_Hz', ...
            'Variance_log','IQR_log','KS_D','KS_p'})]; %#ok<AGROW>
    end
end
disp(distRows);

% Number of units per day
tblN = groupsummary(tblVivo, {'sbjID', 'genotype', 'day'});

% Publication figure (thesis response R2-08): per-genotype RS FR distributions
%   kdMat = mcu_frDist_export(tblLme, pathName)


%% ========================================================================
% Table S6
% =========================================================================
tblIdx = 6;
sheetNames{tblIdx} = ['S' num2str(tblIdx)];
tblInfo{tblIdx} = 'SS Burstiness';
dataSet{tblIdx} = 'In Vivo';
tblPnls{tblIdx} = '3H';

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
% Table S7
% =========================================================================
tblIdx = 7;
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
% Table S8
% =========================================================================
tblIdx = 8;
sheetNames{tblIdx} = ['S' num2str(tblIdx)];
tblInfo{tblIdx} = 'Firing Gain during FRH';
dataSet{tblIdx} = 'MEA';
tblPnls{tblIdx} = '4D, S5D';

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

tblSum = tblMea(tblMea.genotype == 'Control', {'sGain', 'frGain', 'bGain'});
sumStats = groupsummary(tblSum, [], {"mean", "std"}, {'sGain', 'frGain', 'bGain'});



%% ========================================================================
% Table S9
% =========================================================================
tblIdx = 9;
sheetNames{tblIdx} = ['S' num2str(tblIdx)];
tblInfo{tblIdx} = 'Feature Ablation';
dataSet{tblIdx} = 'MEA';
tblPnls{tblIdx} = '4G';

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
