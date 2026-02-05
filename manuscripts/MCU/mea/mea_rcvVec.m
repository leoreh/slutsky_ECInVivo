%% ========================================================================
%  MEA RECOVERY VECTOR ANALYSIS (State Space Trajectory)
%  ========================================================================
%  Analyzes the recovery "strategy" of neurons by decomposing firing rate
%  recovery into Single-Spike and Burst-Spike components.
%
%  Calculates a vector from Baseline (0,0) to Steady State (dSingle, dBurst)
%  and analyzes the Magnitude (Strength of change) and Angle (Strategy).

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


flgPseudo = false; % Validate "Pseudo-Tracking" (Quantile Matching) on MEA data
flgQq = false;


%% ========================================================================
%  PRE-PROCESS: CALCULATE VECTORS
%  ========================================================================

% 1. Calculate Deltas (Steady State - Baseline)
% ---------------------------------------------

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

% frBspk = Baseline Burst Rate
% ss_frBspk = Steady State Burst Rate
tbl.dBrst = tbl.ss_frBspk - tbl.frBspk;
tbl.dSngl = tbl.ss_frSspk - tbl.frSspk;

% Pseudocount
c = 1 / 3600;          

tbl.dBrst = log((tbl.ss_frBspk + c) ./ (tbl.frBspk + c));
tbl.dSngl = log((tbl.ss_frSspk + c) ./ (tbl.frSspk + c));
tbl.dFr = log((tbl.frSs) ./ (tbl.fr));
tbl.dpBspk = (tbl.ss_pBspk_trans) - (tbl.pBspk_trans);

tbl.absBrst = symlog(tbl.ss_frBspk - tbl.frBspk);
tbl.absSngl = symlog(tbl.ss_frSspk - tbl.frSspk);

%% ========================================================================
%  PLOTTING
%  ========================================================================

hFig = mcu_rcvSpace(tbl);

% Resdidual analysis
hFig = mcu_rcvRes(tbl);

% QQ plot
if flgQq
    hFig = mcu_frQq(tbl, 'dataSet', 'mea');
end


% tblGUI_scatHist(tbl, 'xVar', 'pBspk_trans', 'yVar', 'dBrst', 'grpVar', 'Group');
% tblGUI_bar(tbl, 'yVar', 'pBspk', 'xVar', 'Group');


%% ========================================================================
%  ANALYSIS
%  ========================================================================

frml = 'dSngl ~ (dBrst) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');

% Step 1: Regress dSngl against dBrst alone (The mechanical link)
tbl.snglRes = residuals(lmeMdl);

% Step 2: Ask if Baseline Burstiness predicts these residuals
lm_resid = fitlme(tblMcu, 'Resid_Sngl ~ pBspk + (1|Name)');
disp(lm_resid);


tbl.absBrst = (tbl.ss_frBspk - tbl.frBspk);
tbl.absSngl = (tbl.ss_frSspk - tbl.frSspk);
figure; histogram(tbl.absSngl)

frml = 'absSngl ~ (absBrst) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'fitMethod', 'REML', 'flgPlot', false);


frml = 'absBrst ~ (fr + pBspk) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');

frml = 'absSngl ~ (fr + pBspk + absBrst) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');


frml = 'absSngl ~ (absBrst * pBspk_trans) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');


frml = 'absSngl ~ pBspk + fr + (1|Name)';
tblWt = tbl(tbl.Group == 'Control', :);
tblMcu = tbl(tbl.Group == 'MCU-KO', :);
xVar = 'pBspk';
mVar = 'absBrst';
distM = 'normal';
distY = 'normal';
res = lme_mediation(tblMcu, frml, xVar, mVar, 'distM', distM, 'distY', distY);


frml = 'absSngl ~ (absBrst * pBspk) + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tblMcu, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');




frml = 'dFr ~ (pBspk) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');
% MCU-KO neurons have a significant recovery deficit that cannot be
% explained by their change in burstiness.



frml = 'dBrst ~ (pBspk + fr) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'fitMethod', 'ML');
% MCU-KO neurons have a significant recovery deficit that cannot be
% explained by their change in burstiness.



frml = 'dSngl ~ (pBspk + fr + dBrst) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');



frml = 'dSngl ~ (pBspk) * dBrst * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');


frml = 'absSngl ~ (absBrst) * Group * pBspk + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');

frml = 'absSngl ~ (absBrst) * pBspk + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tblWt, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');

[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tblMcu, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');




frml = 'dBrst ~ (pBspk + fr) + Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');

frml = 'dpBspk ~ (pBspk + fr + ss_pBspk) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');




frml = 'frSs ~ (fr + pBspk) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'gamma');

frml = 'frSs ~ (frBspk + frSspk) + Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'gamma');


frml = 'dSngl ~ (pBspk + dBrst + fr) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'fitMethod', 'REML');


frml = 'dSngl ~ pBspk + fr + (1|Name)';
tblWt = tbl(tbl.Group == 'Control', :);
tblMcu = tbl(tbl.Group == 'MCU-KO', :);
xVar = 'pBspk';
mVar = 'dBrst';
distM = 'normal';
distY = 'normal';
res = lme_mediation(tblMcu, frml, xVar, mVar, 'distM', distM, 'distY', distY);


frml = 'dFr ~ (Group) + (1|Name)';
xVar = 'Group';
mVar = 'pBspk';
distM = 'logit-normal';
distY = 'normal';
res = lme_mediation(tbl, frml, xVar, mVar, 'distM', distM, 'distY', distY);


sum(tbl.pBspk(tbl.Group == 'Control') == 0)
sum(tbl.ss_pBspk(tbl.Group == 'Control') == 0)



frml = 'dSngl ~ (pBspk + fr + dBrst) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal');

frml = 'dBrst ~ (pBspk + fr) * Group + (1|Name)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal');



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
