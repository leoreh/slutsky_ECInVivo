%% ========================================================================
%  MEA ANALYSE BURSTS
%  ========================================================================

% Load table
basepaths = [mcu_basepaths('mea_bac')];
presets = {'spktimes', 'rcv', 'frNet'};
[tblFull, ~, ~, v] = mcu_tblMea('basepaths', basepaths, 'presets', presets([1, 3]));

% Prepare subtable of what's needed
tbl = tblFull(:, {'Name', 'Group', 'UnitID', 'fr', 'frSs', 'frTrough', ...
    'uPert', 'spktimes', 'funcon'});

% Clean units that were not perturbed
tbl(~tbl.uPert, :) = [];
tbl = removevars(tbl, 'uPert');

% Spike times from all sessions
spktimes = tbl.spktimes;

% Baseline Window
rcv = catfields([v(:).rcv], 1);
winBsl = rcv.info.winBsl;
winBsl = [min(winBsl(:)), min(winBsl(:, 2))];

% Sweeping Params
isiSweep = [0.005, 0.01, 0.02, 0.05, 0.1];
spkSweep = 2:4;

% Output Matrices
matTStat = nan(length(isiSweep), length(spkSweep));
matCorr  = nan(length(isiSweep), length(spkSweep));
mat0     = nan(length(isiSweep), length(spkSweep));


%% ========================================================================
%  BURST DETECTION SWEEP
%  ========================================================================
fprintf('[BRST_SWEEP] Starting Parameter Sweep (Detection)...\n');

for iIsi = 1 : length(isiSweep)
    for iSpk = 1 : length(spkSweep)

        spkThr = spkSweep(iSpk);
        isiThr = isiSweep(iIsi);
        
        isiEnd = isiThr * 4;
        minIbi = isiEnd * 2;

        % Dynamic Field Name
        fName = sprintf('sib_s%d_i%03d', spkThr, round(isiThr*1000));

        % Burst detection
        brst = brst_detect(spktimes, ...
            'minSpks', spkThr, ...
            'isiStart', isiThr, ...
            'isiEnd', isiEnd, ...
            'minDur', 0, ...
            'minIBI', minIbi, ...
            'flgForce', true, 'flgSave', false, 'flgPlot', false);

        % Burst statistics
        stats = brst_stats(brst, spktimes, 'winCalc', winBsl, ...
            'flgSave', false);

        % Store in Table
        tbl.(fName) = stats.pBspk;
    end

    fprintf('[BRST_SWEEP] Detection: Finished ISI %.3f (%d/%d)...\n', ...
        isiThr, iIsi, length(isiSweep));
end


%% ========================================================================
%  VALIDATION SWEEP
%  ========================================================================

for iIsi = 1 : length(isiSweep)
    for iSpk = 1 : length(spkSweep)

        spkThr = spkSweep(iSpk);
        isiThr = isiSweep(iIsi);

        fName = sprintf('sib_s%d_i%03d', spkThr, round(isiThr*1000));

        % Calculate percent of units with zero burstiness
        currSib = tbl.(fName);
        pZero = sum(currSib == 0) / height(tbl) * 100;

        mat0(iIsi, iSpk) = pZero;

    end
end

% Logit transfrom burstiness
tblVars = tbl.Properties.VariableNames;
bVarsIdx = contains(tblVars, 'sib');
tblLme = tbl_trans(tbl, 'varsInc', tblVars(bVarsIdx), 'logBase', 'logit');


%% ========================================================================
%  PREDICTIVE POWER SWEEP
%  ========================================================================
fprintf('[BRST_SWEEP] Starting Parameter Sweep (Analysis)...\n');

for iIsi = 1 : length(isiSweep)
    for iSpk = 1 : length(spkSweep)

        spkThr = spkSweep(iSpk);
        isiThr = isiSweep(iIsi);

        % Validation Check
        % -----------------------------------------------------------------
        if mat0(iIsi, iSpk) > 50
            continue;
        end

        % -----------------------------------------------------------------
        % Predictive Power (LME)
        % -----------------------------------------------------------------
        fName = sprintf('sib_s%d_i%03d', spkThr, round(isiThr*1000));
        frml = sprintf('frSs ~ %s + (1|Name)', fName);

        [lmeMdl, ~, ~, ~] = lme_analyse(tblLme, frml, ...
            'dist', 'log-normal', ...
            'flgPlot', false, 'verbose', false);

        % Get t-statistic for the specific term
        fxdEffect = lmeMdl.Coefficients;
        matTStat(iIsi, iSpk) = fxdEffect.tStat(2);

        % -----------------------------------------------------------------
        % Correlation (Baseline)
        % -----------------------------------------------------------------
        matCorr(iIsi, iSpk) = corr(tblLme.fr, tblLme.(fName), ...
            'Type', 'Spearman', 'Rows', 'complete');

    end

    fprintf('[BRST_SWEEP] Analysis: Finished ISI %.3f (%d/%d)...\n', ...
        isiThr, iIsi, length(isiSweep));
end


%% ========================================================================
%  PLOT RESULTS
%  ========================================================================

figure('Name', 'Burst Detection Optimization', 'Color', 'w');

% T-Statistic
subplot(1, 3, 1);
heatmap(spkSweep, isiSweep, matTStat);
xlabel('Min Spikes');
ylabel('ISI Threshold (s)');
title('LME Predictive Power (t-stat)');
colormap(gca, parula);

% Correlation
subplot(1, 3, 2);
heatmap(spkSweep, isiSweep, matCorr);
xlabel('Min Spikes');
ylabel('ISI Threshold (s)');
title('Correlation (sib vs fr)');
colormap(gca, parula);

% Zeros
subplot(1, 3, 3);
heatmap(spkSweep, isiSweep, mat0);
xlabel('Min Spikes');
ylabel('ISI Threshold (s)');
title('Percent Zeros');
colormap(gca, parula);

% tblGUI_scatHist(tblLme, 'xVar', 'pBspk', 'yVar', 'fr', 'grpVar', 'Group');


%% ========================================================================
%  FULL MODEL
%  ========================================================================

% % Recovery
% fName = 'sib_s2_i010';
% frml = sprintf('frSs ~ fr + %s + frTrough + funcon + (1|Name)', fName);
% lmeMdl = lme_analyse(tblLme, frml, 'dist', 'log-normal');
% 
% abl = lme_ablation(tblLme, frml, 'dist', 'log-normal', 'flgBkTrans', false);


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