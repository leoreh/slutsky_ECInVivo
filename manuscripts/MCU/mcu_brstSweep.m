%% ========================================================================
%  MEA ANALYSE BURSTS
%  ========================================================================

% Load table
basepaths = [mcu_basepaths('mea_bac'), mcu_basepaths('mea_mcuko')];
presets = {'spktimes', 'rcv', 'frNet'};
[tblFull, ~, ~, v] = mcu_tblMea('basepaths', basepaths, 'presets', presets([1, 3]));

% Prepare subtable of what's needed
tblBrst = tblFull(:, {'Name', 'Group', 'UnitID', 'fr', 'frSs', 'frTrough', ...
    'rcvBsl', 'uPert', 'spktimes', 'funcon'});

% Clean units that were not perturbed
tblBrst(~tblBrst.uPert, :) = [];
tblBrst = removevars(tblBrst, 'uPert');

% Spike times from all sessions
spktimes = tblBrst.spktimes;

% Control units
idxWt = tblBrst.Group == 'Control';

% Baseline Window
rcv = catfields([v(:).rcv], 1);
winBsl = rcv.info.winBsl;
winBsl = [0, min(winBsl(:, 2))];
winBsl = [0, 4300];

% Sweeping Params
isiSweep = [0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2];
spkSweep = 2:6;



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
        % Formula: frSs ~ frB + frS

        fNameSib = sprintf('sib_s%d_i%03d', spkThr, round(isiThr*1000));
        fNameFrB = sprintf('frB_s%d_i%03d', spkThr, round(isiThr*1000));
        fNameFrS = sprintf('frS_s%d_i%03d', spkThr, round(isiThr*1000));

        frml = sprintf('frSs ~ %s + %s', fNameFrB, fNameFrS);

        % Run lightweight ablation (CV)
        abl = lme_ablation(tblWt, frml, 'dist', 'gamma', ...
            'nReps', 5, 'nFolds', 5, ...
            'flgBkTrans', false, 'partitionMode', 'split', ...
            'flgPlot', false);

        % Full R2
        row.fullR2 = abl.pR2(1);

        % Extract dR2 (Unique contribution)
        idxFrB = find(strcmp(abl.vars(2:end), fNameFrB));
        idxFrS = find(strcmp(abl.vars(2:end), fNameFrS));

        if ~isempty(idxFrB), row.dR2_frB = abl.dR2(idxFrB); else, row.dR2_frB = NaN; end
        if ~isempty(idxFrS), row.dR2_frS = abl.dR2(idxFrS); else, row.dR2_frS = NaN; end

        % Clean Outliers / Numerical Instability
        % Detect and remove extreme values (e.g. -1e153) caused by LME
        % convergence failures or singularities.
        varsCheck = {'fullR2', 'dR2_frB', 'dR2_frS'};
        for iVar = 1:length(varsCheck)
            fName = varsCheck{iVar};
            if abs(row.(fName)) > 100
                row.(fName) = NaN;
            end
        end

        % Zeros
        % -----------------------------------------------------------------
        currSib = tblBrst.(fNameSib);
        currSib = currSib(idxWt);
        row.pZero = sum(currSib == 0) / height(tblBrst) * 100;


        % Predictive Power (LME) - Interaction
        % -----------------------------------------------------------------
        % Formula: frSs ~ (fr + sib) * Group + (1|Name)
        frml = sprintf('frSs ~ (fr + %s) * Group + (1|Name)', fNameSib);
        [lmeMdl, ~, ~, ~] = lme_analyse(tblLme, frml, ...
            'dist', 'gamma', ...
            'flgPlot', false, 'verbose', false);

        fxdEffect = lmeMdl.Coefficients;

        % Find indices by name
        idxInt  = find(contains(fxdEffect.Name, ':') & contains(fxdEffect.Name, fNameSib));

        if ~isempty(idxInt),  row.tStatInt  = fxdEffect.tStat(idxInt);  else, row.tStatInt = NaN; end

        % AIC
        row.AIC = lmeMdl.ModelCriterion.AIC;


        % Group Effect (LME)
        % -----------------------------------------------------------------
        % Formula: sib ~ Group * fr + (1|Name)
        frmlGrp = sprintf('%s ~ Group * fr + (1|Name)', fNameSib);
        [lmeGrp, lmeStats, ~, ~] = lme_analyse(tblLme, frmlGrp, ...
            'dist', 'logit-normal', ...
            'flgPlot', false, 'verbose', false);

        fxdGrp = lmeGrp.Coefficients;
        idxGrp = find(strncmpi(fxdGrp.Name, 'Group', 5)); % Find 'Group_...'

        if ~isempty(idxGrp), row.tStatGroup = fxdGrp.tStat(idxGrp(1)); else, row.tStatGroup = NaN; end


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
matTStat     = unstack(tblRes(:, {'spkThr', 'isiThr', 'dR2_frB'}), 'dR2_frB', 'spkThr');
matTStatGrp  = unstack(tblRes(:, {'spkThr', 'isiThr', 'tStatGroup'}), 'tStatGroup', 'spkThr');
matTStatInt  = unstack(tblRes(:, {'spkThr', 'isiThr', 'tStatInt'}), 'tStatInt', 'spkThr');
matAIC       = unstack(tblRes(:, {'spkThr', 'isiThr', 'AIC'}), 'AIC', 'spkThr');
matCorr      = unstack(tblRes(:, {'spkThr', 'isiThr', 'corr'}), 'corr', 'spkThr');
mat0         = unstack(tblRes(:, {'spkThr', 'isiThr', 'pZero'}), 'pZero', 'spkThr');

% Extract matrix data (remove first col which is isiThr label)
matTStat     = table2array(matTStat(:, 2:end));
matTStatGrp  = table2array(matTStatGrp(:, 2:end));
matTStatInt  = table2array(matTStatInt(:, 2:end));
matAIC       = table2array(matAIC(:, 2:end));
matCorr      = table2array(matCorr(:, 2:end));
mat0         = table2array(mat0(:, 2:end));



figure('Name', 'Burst Detection Optimization', 'Color', 'w', 'Position', [100 100 1200 800]);
tiledlayout(2, 3, 'TileSpacing', 'compact');

% 1. Predictive Power (frB vs frSs)
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
% Calculate Shared Variance
% Full = UniqueA + UniqueB + Shared
% Shared = Full - UniqueA - UniqueB
tblRes.dR2_Shared = tblRes.fullR2 - tblRes.dR2_frB - tblRes.dR2_frS;

% Clip negative shared variance (can happen if predictors are correlated in complex ways)
tblRes.dR2_Shared(tblRes.dR2_Shared < 0) = 0;

uSpk = unique(tblRes.spkThr);
clrs = [0.8 0.3 0.3; 0.3 0.3 0.8; 0.7 0.7 0.7]; % Red (Burst), Blue (Single), Gray (Shared)

figure('Name', 'Burst Sweeep - Variance Partition', 'Color', 'w', 'Position', [100 100 1200 400]);
tiledlayout(2, 3, 'TileSpacing', 'compact');

for iThr = 1 : length(uSpk)
    nexttile;

    idx = tblRes.spkThr == uSpk(iThr);
    subTbl = tblRes(idx, :);

    % Ensure unique X-values & Sort by ISI for consistent plotting
    [~, idxUnq] = unique(subTbl.isiThr);
    subTbl = subTbl(idxUnq, :);
    subTbl = sortrows(subTbl, 'isiThr');

    % Data for bar
    yData = [subTbl.dR2_frB, subTbl.dR2_frS, subTbl.dR2_Shared];

    % Evenly spaced bars (Categorical axis)
    xData = 1:height(subTbl);
    b = bar(xData, yData, 'stacked');

    % Colors
    for k = 1:3, b(k).FaceColor = clrs(k, :); end

    title(sprintf('Min Spikes: %d', uSpk(iThr)));
    xlabel('ISI Threshold (s)');
    ylabel('R^2');

    xticks(xData);
    xticklabels(string(subTbl.isiThr));

    if iThr == 1
        legend({'Unique Burst', 'Unique Single', 'Shared'}, 'Location', 'northwest');
    end
    ylim([0 1]);
end




%% ========================================================================
%  FULL MODEL
%  ========================================================================


% Exclude non-bursting
% Note: fNameSib was just a local var in the loop, logic needs to re-define or assume example
fNameFrB = 'frB_s3_i050';
fNameFrS = 'frS_s3_i050';


xVar = 'sib_s3_i050';

frml = sprintf('frSs ~ (fr + %s) * Group + (1|Name)', xVar);

[lmeMdl, lmeStats, lmeInfo, tblMdl] = lme_analyse(tblLme, frml, 'dist', 'gamma');



% Recovery

frml = sprintf('rcvBsl ~ (%s) * Group', fNameSib);

[lmeMdl, lmeStats, lmeInfo, tblMdl] = lme_analyse(tblLme, frml, 'dist', 'log-normal');
lmeMdl

abl = lme_ablation(tblLme, frml, 'dist', 'log-normal', ...
    'flgBkTrans', false, 'partitionMode', 'split');

tblMcu = tblLme(tblLme.Group == 'MCU-KO', :);

frml = sprintf('rcvBsl ~ (fr + %s)', fNameSib);

abl = lme_ablation(tblWt, frml, 'dist', 'log-normal', ...
    'flgBkTrans', false, 'partitionMode', 'split');



%% ========================================================================
%  MEA ISI VALLEY
%  ========================================================================
% Completely useless. No bi-modality. Also, optimization via predictive
% power is preferred.

% presets = {'spktimes'};
% [tblMea, ~, ~, ~] = mcu_tblMea('presets', presets([1]));
%
% idxWt = tblMea.Group == 'Control';
% spktimesMea = tblMea.spktimes(idxWt);
%
% % ISI VALLEY AS THRESHOLD
% isiValMea = brst_isiValley(spktimesMea, 'nSpks', 3);


%% ========================================================================
%  IN VIVO
%  ========================================================================
%
% basepaths = [mcu_basepaths('wt_bsl')];
% [tblVivo, ~, ~, ~] = mcu_tblVivo('basepaths', basepaths, 'presets', {'spktimes'});
%
% idxUnit = tblVivo.UnitType == 'RS';
% spktimesVivo = tblVivo.spktimes(idxUnit);
%
% % ISI VALLEY AS THRESHOLD
% isiValVivo = brst_isiValley(spktimesVivo, 'nSpks', 2);