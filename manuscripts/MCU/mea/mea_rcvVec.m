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

flgPseudo = false; 
if flgPseudo
    % --- CONVERT TO LONG FORMAT ---
    % match_qntl expects: [Group, Name, Day, vars...]
    
    % Baseline Table
    varsTbl = {'fr', 'frBspk', 'frSspk', 'pBspk'};
    tBsl = tbl(:, [{'Group', 'Name'}, varsTbl]);
    tBsl.Day = repmat({'BSL'}, height(tBsl), 1);
    tBsl.Day = categorical(tBsl.Day);

    % Steady State Table (Rename ss_ vars)
    varsSs = {'ss_fr', 'ss_frBspk', 'ss_frSspk', 'pBspk'}; % pBspk roughly same, or use ss_pBspk if available
    % Note: mcu_tblMea provides ss_pBspk? Let's assume we map standard vars.
    % Actually, let's just grab the ss columns and rename.
    tSs = tbl(:, {'Group', 'Name', 'frSs', 'ss_frBspk', 'ss_frSspk', 'ss_pBspk'});
    tSs = renamevars(tSs, {'frSs', 'ss_frBspk', 'ss_frSspk', 'ss_pBspk'}, ...
                           {'fr',    'frBspk',    'frSspk',    'pBspk'});
    tSs.Day = repmat({'BAC3'}, height(tSs), 1); % match_qntl hardcodes 'BAC3' as the 2nd day
    tSs.Day = categorical(tSs.Day);

    % Combine
    tblLong = [tBsl; tSs];
    
    % --- RUN MATCHING ---
    nBins = 80;
    fprintf('VALIDATION: Running match_qntl on MEA data (nBins=%d)...\n', nBins);
    tblSynth = match_qntl(tblLong, nBins, 'flgPool', true, ...
        'vars', {'fr', 'frBspk', 'frSspk', 'pBspk'});
    
    % --- MAP BACK TO MEA FORMAT ---
    % MEA Script Expects: frBspk (BSL), ss_frBspk (SteadyState)
    
    tbl = tblSynth;
    
    % Baseline Mappings
    tbl.fr        = tbl.fr_BSL;      % For Size Scaling
    tbl.pBspk     = tbl.pBspk_BSL;   % For Color
    tbl.frBspk    = tbl.frBspk_BSL;
    tbl.frSspk    = tbl.frSspk_BSL;
    
    % Steady State Mappings
    tbl.ss_frBspk = tbl.frBspk_BAC3;
    tbl.ss_frSspk = tbl.frSspk_BAC3;
    
    fprintf('VALIDATION: Replaced real units with %d Synthetic Units.\n', height(tbl));
end

%% ========================================================================
%  PRE-PROCESS
%  ========================================================================

% Pseudocount
c = 1 / 3600;          

tbl.dBrst_rel = log((tbl.ss_frBspk + c) ./ (tbl.frBspk + c));
tbl.dSngl_rel = log((tbl.ss_frSspk + c) ./ (tbl.frSspk + c));
tbl.dFr = log((tbl.frSs) ./ (tbl.fr));
tbl.dpBspk = (tbl.ss_pBspk_trans) - (tbl.pBspk_trans);

tbl.dBrst_abs = (tbl.ss_frBspk - tbl.frBspk);
tbl.dSngl_abs = (tbl.ss_frSspk - tbl.frSspk);


%% ========================================================================
%  PLOTTING
%  ========================================================================

hFig = mcu_rcvSpace(tbl);

% Resdidual analysis
hFig = mcu_rcvRes(tbl);

% QQ plot
flgQq = false;
if flgQq
    hFig = mcu_frQq(tbl, 'dataSet', 'mea');
end


% tblGUI_scatHist(tbl, 'xVar', 'pBspk_trans', 'yVar', 'dBrst_rel', 'grpVar', 'Group');
% tblGUI_bar(tbl, 'yVar', 'pBspk', 'xVar', 'Group');


%% ========================================================================
%  MEDIATION
%  ========================================================================


% Relative

frml = 'dBrst_rel ~ (pBspk + fr) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'fitMethod', 'ML');


frml = 'dSngl_rel ~ (pBspk + dBrst_rel + fr) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');


frml = 'dSngl_rel ~ pBspk + fr + (1|Name)';
tblWt = tbl(tbl.Group == 'Control', :);
tblMcu = tbl(tbl.Group == 'MCU-KO', :);
xVar = 'pBspk';
mVar = 'dBrst_rel';
distM = 'normal';
distY = 'normal';
res = lme_mediation(tblMcu, frml, xVar, mVar, 'distM', distM, 'distY', distY);





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
frml = 'fr ~ (frSspk + frBspk) * Timepoint * Group + (1|Name)';

[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tblLong, frml, ...
    'dist', 'log-normal');
