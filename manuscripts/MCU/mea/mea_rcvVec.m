%% ========================================================================
%  MEA RECOVERY VECTOR ANALYSIS (State Space Trajectory)
%  ========================================================================
%  Analyzes the recovery "strategy" of neurons by decomposing firing rate
%  recovery into Single-Spike and Burst-Spike components.

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
%  PSEUDO-TRACKING
%  ========================================================================
% Validate "Pseudo-Tracking" (Quantile Matching) on MEA data

flgPseudo = true; 
if flgPseudo
    
    nBins = 50;

    % --- ASSERT EQUAL DISTRIBUTIONS ---
    cutoff_z = -1.5; % User defined cut-off (Visual Inspection)
    % cutoff_z = -6;

    hFig = mcu_rcvQq(tbl, 'var', 'frSspk', 'dataSet', 'mea', 'cutoff_z', cutoff_z);

    % Define threshold based on Control Baseline
    frRef = log(tbl.fr(tbl.Group == 'Control'));
    mu = mean(frRef);
    sigma = std(frRef);
    cutoff_Hz = exp(mu + cutoff_z * sigma);

    cutoff_p = normcdf(cutoff_z); % For display only

    % Apply Filter
    nBefore = height(tbl);

    % Identify Units to keep (Must have BSL FR > cutoff)
    goodIdx = tbl.fr >= cutoff_Hz;

    % Filter original table
    tbl = tbl(goodIdx, :);

    fprintf('\n================================================================\n');
    fprintf(' FILTERING LOW FR UNITS\n');
    fprintf('================================================================\n');
    fprintf('Cut-off Z-score: %.2f (%.1f%%)\n', cutoff_z, cutoff_p*100);
    fprintf('Cut-off FR     : %.4f Hz\n', cutoff_Hz);
    fprintf('Removed %d units (%.1f%%). Remaining: %d units\n', ...
        sum(~goodIdx), sum(~goodIdx)/length(goodIdx)*100, sum(goodIdx));
    fprintf('================================================================\n\n');


    % --- CONVERT TO LONG FORMAT ---
    
    % Baseline Table
    varsTbl = {'fr', 'frBspk', 'frSspk', 'pBspk', 'pBspk_trans'};
    tBsl = tbl(:, [{'Group', 'Name'}, varsTbl]);
    tBsl.Day = repmat({'BSL'}, height(tBsl), 1);
    tBsl.Day = categorical(tBsl.Day);

    % Steady State Table (Rename ss_ vars)
    varsSs = {'frSs', 'ss_frBspk', 'ss_frSspk', 'ss_pBspk', 'ss_pBspk_trans'}; 
    tSs = tbl(:, [{'Group', 'Name'}, varsSs]);
    tSs = renamevars(tSs, {'frSs', 'ss_frBspk', 'ss_frSspk', 'ss_pBspk', 'ss_pBspk_trans'}, ...
                           {'fr',    'frBspk',    'frSspk',    'pBspk', 'pBspk_trans'});
    tSs.Day = repmat({'BAC3'}, height(tSs), 1); % match_qntl hardcodes 'BAC3' as the 2nd day
    tSs.Day = categorical(tSs.Day);

    % Combine
    tblLong = [tBsl; tSs];
    
    % --- RUN MATCHING ---
    fprintf('VALIDATION: Running match_qntl on MEA data (nBins=%d)...\n', nBins);
    tbl = match_qntl(tblLong, nBins, 'flgPool', true, ...
        'var', 'fr', 'avgType', 'mean');
    
    tbl.frSs = tbl.ss_fr;
    fprintf('VALIDATION: Replaced real units with %d Synthetic Units.\n', height(tbl));
end

%% ========================================================================
%  PRE-PROCESS
%  ========================================================================

tbl.dBrst_rel = log((tbl.ss_frBspk) ./ (tbl.frBspk));
tbl.dSngl_rel = log((tbl.ss_frSspk) ./ (tbl.frSspk));
tbl.dFr = log((tbl.frSs) ./ (tbl.fr));
tbl.dpBspk = (tbl.ss_pBspk_trans) - (tbl.pBspk_trans);

tbl.dBrst_abs = (tbl.ss_frBspk - tbl.frBspk);
tbl.dSngl_abs = (tbl.ss_frSspk - tbl.frSspk);

% PRISM

idxGrp = tbl.Group == 'MCU-KO';
prismVars = {'Name', 'dSngl_rel', 'dBrst_rel', 'dSngl_abs', 'dBrst_abs', 'dFr', 'dpBspk', 'pBspk', 'pBspk_trans'};
tblPrism = tbl(idxGrp, prismVars);
string(prismVars)

%% ========================================================================
%  PLOTTING
%  ========================================================================

hFig = mcu_rcvSpace(tbl);

% Resdidual analysis
hFig = mcu_rcvRes(tbl);

frml = 'dSngl_rel ~ (dBrst_rel + pBspk + fr) * Group + (1|Name)';
lmeMdl = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'verbose', true);

[tblRes, hFig] = lme_pr(lmeMdl, 'dBrst_rel', ...
    'flgMode', 'regression', ...
    'varGrp',  'Group');


% tblGUI_scatHist(tbl, 'xVar', 'pBspk_trans', 'yVar', 'dBrst_rel', 'grpVar', 'Group');
% tblGUI_bar(tbl, 'yVar', 'pBspk', 'xVar', 'Group');



%% ========================================================================
%  Single Group Models
%  ========================================================================

idxWt = tbl.Group == 'Control';
tblWt = tbl(idxWt, :);

frml = 'dBrst_rel ~ (fr + pBspk * dSngl_rel) + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tblWt, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');

frml = 'dSngl_rel ~ (fr + pBspk * dBrst_rel) + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tblWt, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');

frml = 'dFr ~ (frBspk + frSspk) + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');

abl = lme_ablation(tblWt, frml, 'dist', 'normal', ...
    'flgBkTrans', false, 'partitionMode', 'split', 'nrep', 5);


%% ========================================================================
%  MEDIATION
%  ========================================================================


% Relative

frml = 'dBrst_rel ~ (pBspk + fr) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');


% PLOT INTERACTION
vars = {'pBspk', 'Group'};
hFig = plot_axSize('flgFullscreen', true, 'flgPos', true);

% Partial Dependence
hAx = nexttile;
[pdRes, hFig] = lme_pd(lmeMdl, vars, 'transParams', lmeInfo.transParams, ...
    'hAx', hAx);


frml = 'dSngl_rel ~ (pBspk + fr) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');


frml = 'dSngl_rel ~ (pBspk + fr) + (1|Name)';
tblWt = tbl(tbl.Group == 'Control', :);
tblMcu = tbl(tbl.Group == 'MCU-KO', :);
xVar = 'pBspk';
mVar = 'dBrst_rel';
distM = 'normal';
distY = 'normal';
res = lme_mediation(tblWt, frml, xVar, mVar, 'distM', distM, 'distY', distY);



tblMcu.frSspk





% Absolute

frml = 'dSngl_abs ~ (dBrst_abs) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'fitMethod', 'REML', 'flgPlot', false);


frml = 'dBrst_abs ~ (fr + pBspk) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');

frml = 'dSngl_abs ~ (fr + pBspk + dBrst_abs) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');




frml = 'dSngl_abs ~ frSspk + pBspk + (1|Name)';
tblWt = tbl(tbl.Group == 'Control', :);
tblMcu = tbl(tbl.Group == 'MCU-KO', :);
xVar = 'pBspk';
mVar = 'dBrst_abs';
distM = 'normal';
distY = 'normal';
res = lme_mediation(tblWt, frml, xVar, mVar, 'distM', distM, 'distY', distY);


frml = 'dSngl_abs ~ (dBrst_abs * pBspk) + fr + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tblMcu, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');

frml = 'dSngl_rel ~ (dBrst_rel * pBspk) + fr + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tblWt, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');













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
