%% ========================================================================
%  MEA ANALYSE BURSTS
%  ========================================================================

% Load table
basepaths = [mcu_basepaths('mea_bac'), mcu_basepaths('mea_mcuko')];
presets = {'spktimes', 'rcv', 'frNet'};
[tblFull, ~, ~, v] = mcu_tblMea('basepaths', basepaths, 'presets', presets([1, 3]));

% Prepare subtable of what's needed
tblBrst = tblFull(:, {'Name', 'Group', 'UnitID', 'fr', 'ss_fr', 'frAcute', ...
    'rcvBsl', 'spktimes', 'funcon'});

% Spike times from all sessions
spktimes = tblBrst.spktimes;

% Control units
idxWt = tblBrst.Group == 'Control';
idxMcu = tblBrst.Group == 'MCU-KO';

% Baseline Window
rcv = catfields([v(:).rcv], 1);
winBsl = rcv.info.winBsl;
winBsl = [0, min(winBsl(:, 2))];
winBsl = [0, 4000];

% Sweeping Params
isiSweep = [0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.1];
isiSweep = [0.01, 0.02, 0.03, 0.05, 0.1];
spkSweep = 3:5;


%% ========================================================================
%  BURST DETECTION SWEEP
%  ========================================================================
fprintf('[BRST_SWEEP] Starting Parameter Sweep (Detection)...\n');

for iIsi = 1 : length(isiSweep)
    for iSpk = 1 : length(spkSweep)

        spkThr = spkSweep(iSpk);
        isiThr = isiSweep(iIsi);

        isiEnd = isiThr * 2;
        minIbi = isiEnd * 1;
        minDur = 0;

        % Dynamic Field Name
        fNameSib = sprintf('sib_s%d_i%03d', spkThr, round(isiThr*1000));
        fNameFrB = sprintf('frB_s%d_i%03d', spkThr, round(isiThr*1000));
        fNameFrS = sprintf('frS_s%d_i%03d', spkThr, round(isiThr*1000));

        % Burst detection
        brst = brst_detect(spktimes, ...
            'minSpks', spkThr, ...
            'isiStart', isiThr, ...
            'isiEnd', isiEnd, ...
            'minDur', minDur, ...
            'minIBI', minIbi, ...
            'flgForce', true, 'flgSave', false, 'flgPlot', false);

        % Burst statistics
        stats = brst_stats(brst, spktimes, 'winCalc', winBsl, ...
            'flgSave', false);

        % Store in Table
        tblBrst.(fNameSib) = stats.pBspk;
        tblBrst.(fNameFrB) = stats.frBspk;
        tblBrst.(fNameFrS) = stats.frSspk;
    end

    fprintf('[BRST_SWEEP] Detection: Finished ISI %.3f (%d/%d)...\n', ...
        isiThr, iIsi, length(isiSweep));
end


%% ========================================================================
%  VALIDATION SWEEP
%  ========================================================================

% Logit transfrom burstiness
% tblVars = tblBrst.Properties.VariableNames;
% bVarsIdx = contains(tblVars, 'sib');
% tblLme = tbl_trans(tblBrst, 'varsInc', tblVars(bVarsIdx), 'logBase', 'logit');
tblLme = tblBrst;

% Grab WT data
tblWt = tblLme(idxWt, :);
tblMcu = tblLme(idxMcu, :);


%% ========================================================================
%  PREDICTIVE POWER SWEEP
%  ========================================================================
tblRes = table();
fprintf('[BRST_SWEEP] Starting Parameter Sweep (Analysis)...\n');

for iIsi = 1 : length(isiSweep)
    for iSpk = 1 : length(spkSweep)

        spkThr = spkSweep(iSpk);
        isiThr = isiSweep(iIsi);

        row = struct();
        row.spkThr = spkThr;
        row.isiThr = isiThr;

        % -----------------------------------------------------------------
        % Predictive Power (LME) - Ablation (WT Only)
        % -----------------------------------------------------------------
        % Formula: ss_fr ~ frB + frS

        fNameSib = sprintf('sib_s%d_i%03d', spkThr, round(isiThr*1000));
        fNameFrB = sprintf('frB_s%d_i%03d', spkThr, round(isiThr*1000));
        fNameFrS = sprintf('frS_s%d_i%03d', spkThr, round(isiThr*1000));

        frml = sprintf('ss_fr ~ %s + %s', fNameFrB, fNameFrS);

        % Run ablation (CV) - Control
        ablWt = lme_ablation(tblWt, frml, 'dist', 'log-normal', ...
            'nReps', 10, 'nFolds', 5, ...
            'flgBkTrans', false, 'partitionMode', 'split', ...
            'flgPlot', false);

        idxFrB = find(strcmp(ablWt.vars(2:end), fNameFrB));
        idxFrS = find(strcmp(ablWt.vars(2:end), fNameFrS));
        vB = NaN; vS = NaN;
        if ~isempty(idxFrB), vB = ablWt.dR2(idxFrB); end
        if ~isempty(idxFrS), vS = ablWt.dR2(idxFrS); end
        row.dR2_wt = [vB, vS, ablWt.dR2(end)];

        % Run ablation (CV) - MCU-KO
        ablMcu = lme_ablation(tblMcu, frml, 'dist', 'log-normal', ...
            'nReps', 10, 'nFolds', 5, ...
            'flgBkTrans', false, 'partitionMode', 'split', ...
            'flgPlot', false);

        idxFrB = find(strcmp(ablMcu.vars(2:end), fNameFrB));
        idxFrS = find(strcmp(ablMcu.vars(2:end), fNameFrS));
        vB = NaN; vS = NaN;
        if ~isempty(idxFrB), vB = ablMcu.dR2(idxFrB); end
        if ~isempty(idxFrS), vS = ablMcu.dR2(idxFrS); end
        row.dR2_mcu = [vB, vS, ablMcu.dR2(end)];


        % Zeros
        % -----------------------------------------------------------------
        currSib = tblBrst.(fNameSib);
        currSib = currSib(idxWt);
        row.pZero = sum(currSib == 0) / height(currSib) * 100;


        % Predictive Power (LME) - Interaction
        % -----------------------------------------------------------------
        
        % Formula: ss_fr ~ (fr + sib) * Group + (1|Name)
        frml = sprintf('ss_fr ~ (fr + %s) * Group + (1|Name)', fNameSib);
        [lmeMdl, ~, ~, ~] = lme_analyse(tblLme, frml, ...
            'dist', 'log-normal', 'fitMethod', 'ML', ...
            'flgPlot', false, 'verbose', false);

        fxdEffect = lmeMdl.Coefficients;

        % Find indices by name
        idxInt  = find(contains(fxdEffect.Name, ':') & contains(fxdEffect.Name, fNameSib));

        if ~isempty(idxInt),  row.tStatInt  = fxdEffect.tStat(idxInt);  else, row.tStatInt = NaN; end

        % AIC
        row.AIC = lmeMdl.ModelCriterion.AIC;


        % Group Effect (LME)
        % -----------------------------------------------------------------
        % Formula: sib ~ Group + (1|Name)
        frmlGrp = sprintf('%s ~ Group + (1|Name)', fNameSib);
        [lmeGrp, lmeStats, ~, ~] = lme_analyse(tblLme, frmlGrp, ...
            'dist', 'logit-normal', 'fitMethod', 'ML', ...
            'flgPlot', false, 'verbose', false);

        fxdGrp = lmeGrp.Coefficients;
        idxGrp = find(strncmpi(fxdGrp.Name, 'Group', 5)); 
        if ~isempty(idxGrp)
            row.tStatGroup = fxdGrp.tStat(idxGrp(1));
        else
            row.tStatGroup = NaN;
        end

        % Correlation (Baseline)
        % -----------------------------------------------------------------
        row.corr = corr(tblWt.fr, tblWt.(fNameSib), ...
            'Type', 'Spearman', 'Rows', 'complete');

        % Store
        tblRes = [tblRes; struct2table(row)];

    end

    fprintf('[BRST_SWEEP] Analysis: Finished ISI %.3f (%d/%d)...\n', ...
        isiThr, iIsi, length(isiSweep));
end


%% ========================================================================
%  PLOT RESULTS
%  ========================================================================

% tblGUI_scatHist(tblBrst, 'xVar', 'pBspk', 'yVar', 'fr', 'grpVar', 'Group');
% tblGUI_bar(tblBrst, 'xVar', 'Group', 'yVar', 'fr');


% Convert Table to Matrices for Heatmaps
matTStat     = unstack(tblRes(:, {'spkThr', 'isiThr', 'dR2_wt'}), 'dR2_wt', 'spkThr');
matTStat     = cell2mat(table2array(matTStat(:, 2:end)));
matTStat     = matTStat(:, 1:3:end); % Extract Unique Burst (Component 1) for each spkThr group
matTStatGrp  = unstack(tblRes(:, {'spkThr', 'isiThr', 'tStatGroup'}), 'tStatGroup', 'spkThr');
matTStatInt  = unstack(tblRes(:, {'spkThr', 'isiThr', 'tStatInt'}), 'tStatInt', 'spkThr');
matAIC       = unstack(tblRes(:, {'spkThr', 'isiThr', 'AIC'}), 'AIC', 'spkThr');
matCorr      = unstack(tblRes(:, {'spkThr', 'isiThr', 'corr'}), 'corr', 'spkThr');
mat0         = unstack(tblRes(:, {'spkThr', 'isiThr', 'pZero'}), 'pZero', 'spkThr');

% Extract matrix data (remove first col which is isiThr label)
matTStatGrp  = table2array(matTStatGrp(:, 2:end));
matTStatInt  = table2array(matTStatInt(:, 2:end));
matAIC       = table2array(matAIC(:, 2:end));
matCorr      = table2array(matCorr(:, 2:end));
mat0         = table2array(mat0(:, 2:end));



figure('Name', 'Burst Detection Optimization', 'Color', 'w', 'Position', [100 100 1200 800]);
tiledlayout(2, 3, 'TileSpacing', 'compact');

% 1. Predictive Power (frB vs ss_fr)
nexttile;
heatmap(spkSweep, isiSweep, matTStat, 'ColorMap', parula);
xlabel('Min Spikes'); ylabel('ISI Threshold (s)');
title('Predictive Power: \Delta R^2 (frB)');

% % 2. T-Statistic (Group Effect)
nexttile;
heatmap(spkSweep, isiSweep, matTStatGrp, 'ColorMap', parula);
xlabel('Min Spikes'); ylabel('ISI Threshold (s)');
title('Difference: t-stat (Group)');


% 3. T-Statistic (Interaction)
nexttile;
heatmap(spkSweep, isiSweep, matTStatInt, 'ColorMap', parula);
xlabel('Min Spikes'); ylabel('ISI Threshold (s)');
title('Inference: t-stat (Interaction)');

% 4. Correlation
nexttile;
heatmap(spkSweep, isiSweep, matCorr, 'ColorMap', parula);
xlabel('Min Spikes'); ylabel('ISI Threshold (s)');
title('Correlation (sib vs fr)');

% 5. Zeros
nexttile;
heatmap(spkSweep, isiSweep, mat0, 'ColorMap', parula);
xlabel('Min Spikes'); ylabel('ISI Threshold (s)');
title('Percent Zeros');

% 6. Model Fit: AIC
nexttile;
heatmap(spkSweep, isiSweep, matAIC, 'ColorMap', flipud(parula));
xlabel('Min Spikes'); ylabel('ISI Threshold (s)');
title('Model Fit: AIC (Lower is Better)');


% -------------------------------------------------------------------------
% STACKED BAR PLOTS (Variance Partitioning)
% -------------------------------------------------------------------------


% Clip negative shared variance (can happen if predictors are correlated in complex ways)
tblRes.dR2_wt(tblRes.dR2_wt(:, 3) < 0, 3) = 0;
tblRes.dR2_mcu(tblRes.dR2_mcu(:, 3) < 0, 3) = 0;

uSpk = unique(tblRes.spkThr);
clrs = [0.8 0.3 0.3; 0.3 0.3 0.8; 0.7 0.7 0.7]; % Red (Burst), Blue (Single), Gray (Shared)

figure('Name', 'Burst Sweeep - Variance Partition', 'Color', 'w', 'Position', [100 100 1200 800]);
tiledlayout(2, length(uSpk), 'TileSpacing', 'compact');

for iGrp = 1:2
    for iThr = 1 : length(uSpk)
        nexttile;

        idx = tblRes.spkThr == uSpk(iThr);
        subTbl = tblRes(idx, :);

        % Ensure unique X-values & Sort by ISI for consistent plotting
        [~, idxUnq] = unique(subTbl.isiThr);
        subTbl = subTbl(idxUnq, :);
        subTbl = sortrows(subTbl, 'isiThr');

        % Data for bar
        if iGrp == 1
            yData = subTbl.dR2_wt;
            grpName = 'Control';
        else
            yData = subTbl.dR2_mcu;
            grpName = 'MCU-KO';
        end

        % Evenly spaced bars (Categorical axis)
        xData = 1:height(subTbl);
        b = bar(xData, yData, 'stacked');

        % Colors
        for k = 1:3, b(k).FaceColor = clrs(k, :); end

        if iGrp == 1
            title(sprintf('Min Spikes: %d', uSpk(iThr)));
        else
            title(sprintf('%s (Spikes: %d)', grpName, uSpk(iThr)));
        end
        
        if iThr == 1
            ylabel({grpName, 'R^2'});
        end
        if iGrp == 2
            xlabel('ISI Threshold (s)');
        end

        xticks(xData);
        xticklabels(string(subTbl.isiThr));

        if iGrp == 1 && iThr == 1
            legend({'Unique Burst', 'Unique Single', 'Shared'}, 'Location', 'northwest');
        end
        ylim([0 0.5]);
    end
end




