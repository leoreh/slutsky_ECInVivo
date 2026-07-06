%% ========================================================================
%  MCU REGRESSION MODELS
%  ========================================================================
% Various LME models for analyzing FRH recovery dynamics, speficially the
% relationship between burstiness (BSL), burst- and single-firing (SS).
% Including:
% State Space Trajectory 
% Residuals (partial regression)
% LSmeans (least-squares means)
% Mediation
% Ablation
% 
% Can work on both in vivo and mea tables (tbl)


%% ========================================================================
%  LOAD DATA
%  ========================================================================

% Load with steady state variables
presets = {'steadyState'};
[tbl, xVec, basepaths, v] = mcu_tblMea('presets', presets, 'flgOtl', true);

% Add logit pBurst
tblTrans = tbl_trans(tbl, 'varsInc', {'pBurst', 'ss_pBurst'}, 'logBase', 'logit');
tbl.pBurst_trans = tblTrans.pBurst;
tbl.ss_pBurst_trans = tblTrans.ss_pBurst;


%% ========================================================================
%  PRE-PROCESS
%  ========================================================================

% Relative
tbl.bGain = log((tbl.ss_frBurst) ./ (tbl.frBurst));
tbl.sGain = log((tbl.ss_frSingle) ./ (tbl.frSingle));
tbl.frGain = log((tbl.ss_fr) ./ (tbl.fr));
tbl.pBurstGain = (tbl.ss_pBurst_trans) - (tbl.pBurst_trans);

% Absolute
tbl.bDelta = (tbl.ss_frBurst - tbl.frBurst);
tbl.sDelta = (tbl.ss_frSingle - tbl.frSingle);
tbl.frDelta = (tbl.ss_fr - tbl.fr);



%% ========================================================================
%  EXCHANGE RATE
%  ========================================================================
% Use raw values (not standardized effects)

frml = 'sDelta ~ (bDelta) * genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'verbose', true, 'flgStnd', false);

frml = 'sGain ~ bGain * genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'verbose', true, 'flgStnd', false);

%% ========================================================================
%  PLOTTING
%  ========================================================================

hFig = mcu_rcvSpace(tblMea);

frml = 'sGain ~ bGain * genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'verbose', true);

frml = 'sDelta ~ bDelta * genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'verbose', true);


% guiTbl_scatHist(tbl, 'xVar', 'pBurst_trans', 'yVar', 'bGain', 'grpVar', 'genotype');
% guiTbl_bar(tbl, 'yVar', 'pBurst', 'xVar', 'genotype');

%% ========================================================================
%  RESIDUAL ANALYSIS
%  ========================================================================
% Partial Regression


hFig = figure;
hTile = tiledlayout(hFig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
hAx = nexttile;

frml = 'sGain ~ (pBurst + bGain) * genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'verbose', false, 'flgStnd', false);

[tblRes, hFig] = lme_pr(lmeMdl, 'bGain', ...
    'flgMode', 'regression', 'hAx', hAx, ...
    'varGrp',  'genotype', 'transParams', lmeInfo.transParams);

% frml = 'ResidY ~ ResidX * genotype + (1|sbjID)';
% [lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tblRes, frml, ...
%     'dist', 'normal', 'verbose', true, 'flgStnd', false);

hAx = nexttile;
frml = 'sDelta ~ (pBurst * fr * bDelta) * genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'verbose', false);

[tblRes, hFig] = lme_pr(lmeMdl, 'bDelta', ...
    'flgMode', 'regression', 'hAx', hAx, ...
    'varGrp',  'genotype', 'transParams', lmeInfo.transParams);

% Prism
idxGrp = tblRes.genotype == 'Control';
tblRes.ResidY(idxGrp);

% % Obsolete
% hFig = mcu_rcvRes(tbl);








%% ========================================================================
%  LEAST-SQUARES MEANS (LSMEANS)
%  ========================================================================

% Steady-state raw values
dist = 'log-normal';

frml = 'ss_frBurst ~ (fr + pBurst) * genotype + (1|sbjID)';
[lmeMdl1, lmeStats, lmeInfo1, ~] = lme_analyse(tbl, frml, ...
    'dist', dist, 'flgStnd', false);

frml = 'ss_frSingle ~ (fr + pBurst) * genotype + (1|sbjID)';
[lmeMdl2, lmeStats, lmeInfo2, ~] = lme_analyse(tbl, frml, ...
    'dist', dist, 'flgStnd', false);

frml = 'ss_fr ~ (fr + pBurst) * genotype + (1|sbjID)';
[lmeMdl3, lmeStats, lmeInfo3, ~] = lme_analyse(tbl, frml, ...
    'dist', dist, 'flgStnd', false);

% ---
% Plot Interaction
hFig = plot_axSize('flgFullscreen', true, 'flgPos', true);

vars = {'pBurst', 'genotype'};
hAx = nexttile;
[pdRes, hFig] = lme_lsmeans(lmeMdl1, vars, 'transParams', lmeInfo1.transParams, ...
    'hAx', hAx);
set(gca, "YScale", "log")

hAx = nexttile;
[pdRes, hFig] = lme_lsmeans(lmeMdl2, vars, 'transParams', lmeInfo2.transParams, ...
    'hAx', hAx);
set(gca, "YScale", "log")

hAx = nexttile;
[pdRes, hFig] = lme_lsmeans(lmeMdl3, vars, 'transParams', lmeInfo3.transParams, ...
    'hAx', hAx);
set(gca, "YScale", "log")

vars = {'fr', 'genotype'};
hAx = nexttile;
[pdRes, hFig] = lme_lsmeans(lmeMdl1, vars, 'transParams', lmeInfo1.transParams, ...
    'hAx', hAx);
set(gca, "YScale", "log")
set(gca, "XScale", "log")

hAx = nexttile;
[pdRes, hFig] = lme_lsmeans(lmeMdl2, vars, 'transParams', lmeInfo2.transParams, ...
    'hAx', hAx);
set(gca, "YScale", "log")
set(gca, "XScale", "log")

hAx = nexttile;
[pdRes, hFig] = lme_lsmeans(lmeMdl3, vars, 'transParams', lmeInfo3.transParams, ...
    'hAx', hAx);
set(gca, "YScale", "log")
set(gca, "XScale", "log")



% To Prism
grpIdx = pdRes.genotype == "Control";
prismMat = [pdRes(grpIdx, vars(1)), ...
    pdRes(grpIdx, {'ss_fr_pred', 'ss_fr_upper', 'ss_fr_lower'})];
grpIdx = pdRes.genotype == "MCU-KO";
prismMat = [pdRes(grpIdx, vars(1)), ...
    pdRes(grpIdx, {'ss_fr_pred', 'ss_fr_upper', 'ss_fr_lower'})];


%% ========================================================================
%  MEDIATION
%  ========================================================================
% NOTE: Formal mediation analysis was removed. The decomposition of
% pBurst effects into component-specific paths (burst-spike vs single-spike)
% is compositionally confounded: pBurst defines the baseline split
% (frBurst = fr * pBurst), so any model predicting a component from pBurst
% conflates biological effects with compositional arithmetic. The
% localization of deficits to firing components relies instead on the
% factorial analysis (S7), feature ablation (S10), and the allocation
% formalization (Supp. Note 1, Equations 1-4).


frml = 'ss_frSingle ~ pBurst + fr + (1|sbjID)';
xVar = 'pBurst';
mVar = 'ss_frBurst';
distM = 'log-normal';
distY = distM;

% WT
tblWt = tbl(tbl.genotype == 'Control', :);
resWt = lme_mediation(tblWt, frml, xVar, mVar, 'distM', distM, 'distY', distY);

resWt.plot.X = tblWt.pBurst_trans;
lme_mediationPlot(resWt)

% MCU
tblMcu = tbl(tbl.genotype == 'MCU-KO', :);
resMcu = lme_mediation(tblMcu, frml, xVar, mVar, 'distM', distM, 'distY', distY);

resMcu.plot.X = tblMcu.pBurst_trans;
lme_mediationPlot(resMcu)

% COMBINED MODELS (unstandardized)

% X -> M
frml = 'ss_frBurst ~ (fr + pBurst) * genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', dist, 'verbose', true, 'flgStnd', false);

% X -> Y
frml = 'ss_frSingle ~ (fr + pBurst) * genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tblWt, frml, ...
    'dist', dist, 'verbose', true, 'flgStnd', false);

% X -> Y | M 
frml = 'ss_frSingle ~ (fr + pBurst + ss_frBurst) * genotype + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'log-normal', 'verbose', true, 'flgStnd', false);


%% ========================================================================
%  MEDIATION GAIN
%  ========================================================================

hFig = plot_axSize('flgFullscreen', true, 'flgPos', true);

[mdlBC, statsBC, infoBC] = lme_analyse(tblMea, ...
    'ss_frSingle ~ (pBurst + ss_frBurst) * genotype  + (1|sbjID)', ...
    'dist', 'log-normal', 'flgStnd', false, 'transTemplate', tmpl);

hAx = nexttile;
[pdRes, hFig] = lme_lsmeans(mdlBC, {'pBurst', 'genotype'}, ...
    'hAx', hAx);

hFig = plot_axSize('flgFullscreen', true, 'flgPos', true);

% Model
[mdlBC, statsBC, infoBC] = lme_analyse(tblMea, ...
    'sGain ~ (pBurst + bGain) * genotype  + (1|sbjID)', ...
    'dist', 'log-normal', 'flgStnd', false, 'transTemplate', tmpl);
hAx = nexttile;
[pdRes, hFig] = lme_lsmeans(mdlBC, {'pBurst', 'genotype'}, ...
    'hAx', hAx);

hFig = plot_axSize('flgFullscreen', true, 'flgPos', true);

% Model
[mdlBC, statsBC, infoBC] = lme_analyse(tblMea, ...
    'sGain ~ (pBurst + bGain) * genotype  + (1|sbjID)', ...
    'dist', 'log-normal', 'flgStnd', false, 'transTemplate', tmpl);
hAx = nexttile;
[pdRes, hFig] = lme_lsmeans(mdlBC, {'pBurst', 'genotype'}, ...
    'hAx', hAx);

% Model
[mdlBC, statsBC, infoBC] = lme_analyse(tblMea, ...
    'bGain ~ (pBurst + fr) * genotype  + (1|sbjID)', ...
    'dist', 'log-normal', 'flgStnd', false, 'transTemplate', tmpl);
hAx = nexttile;
[pdRes, hFig] = lme_lsmeans(mdlBC, {'pBurst', 'genotype'}, ...
    'hAx', hAx);

% Model
hFig = plot_axSize('flgFullscreen', true, 'flgPos', true);
[mdlBC, statsBC, infoBC] = lme_analyse(tblMea, ...
    'frGain ~ (pBurst + fr) * genotype  + (1|sbjID)', ...
    'dist', 'normal', 'flgStnd', false, 'transTemplate', tmpl);
hAx = nexttile;
[pdRes, hFig] = lme_lsmeans(mdlBC, {'fr', 'genotype'}, ...
    'hAx', hAx);



%%% ---

hFig = plot_axSize('flgFullscreen', true, 'flgPos', true);

% Model
[mdlBC, statsBC, infoBC] = lme_analyse(tblMea, ...
    'ss_frBurst ~ (pBurst + fr + ss_frSingle) * genotype  + (1|sbjID)', ...
    'dist', 'log-normal', 'flgStnd', false, 'transTemplate', tmpl);
hAx = nexttile;
[pdRes, hFig] = lme_lsmeans(mdlBC, {'pBurst', 'genotype'}, ...
    'hAx', hAx);

% Model
[mdlBC, statsBC, infoBC] = lme_analyse(tblMea, ...
    'ss_frSingle ~ (pBurst + fr + ss_frBurst) * genotype  + (1|sbjID)', ...
    'dist', 'log-normal', 'flgStnd', false, 'transTemplate', tmpl);
hAx = nexttile;
[pdRes, hFig] = lme_lsmeans(mdlBC, {'pBurst', 'genotype'}, ...
    'hAx', hAx);

% Model
[mdlBC, statsBC, infoBC] = lme_analyse(tblMea, ...
    'ss_frSingle ~ (pBurst + fr) * genotype  + (1|sbjID)', ...
    'dist', 'log-normal', 'flgStnd', false, 'transTemplate', tmpl);
hAx = nexttile;
[pdRes, hFig] = lme_lsmeans(mdlBC, {'fr', 'genotype'}, ...
    'hAx', hAx);


%% ========================================================================
%  ABLATION
%  ========================================================================


frml = 'ss_fr ~ (frBurst + frSingle) + (1 | sbjID)';
% frml = 'ss_fr ~ (frBurst + frSingle)';
% frml = 'ss_fr ~ (fr + pBurst) + (1 | sbjID)';
dist = 'log-normal';

nRep = 5;
partMode = 'split';

tblWt = tbl(tbl.genotype == 'Control', :);
abl = lme_ablation(tblWt, frml, 'dist', dist, ...
    'flgBkTrans', false, 'partitionMode', partMode, 'nrep', nRep);

tblMcu = tbl(tbl.genotype == 'MCU-KO', :);
abl = lme_ablation(tblMcu, frml, 'dist', dist, ...
    'flgBkTrans', false, 'partitionMode', partMode, 'nrep', nRep);


% With Group
tblLme = tbl;

frml = 'ss_fr ~ (frBurst + frSingle) * genotype + (1 | sbjID)';
frml = 'ss_fr ~ (fr + pBurst) * genotype + (1 | sbjID)';

nRep = 5;
dist = 'log-normal';

abl = lme_ablation(tblLme, frml, 'dist', dist, ...
    'flgBkTrans', false, 'partitionMode', 'split', 'nrep', nRep);


% tblLme.genotype = double(tbl.genotype) - 1;
% dist = 'binomial';
% frml = 'genotype ~ (ss_fr + fr + pBurst) + (1 | sbjID)';

%% ========================================================================
%  CONTRIBUTION ANALYSIS: BSL vs SS
%  ========================================================================
%  Test if the relative contribution of burst spikes (pBurst) changes between
%  Baseline and Steady State (SS), and if this change differs by Group.

fprintf('\n================================================================\n');
fprintf(' STATISTICS: CONTRIBUTION ANALYSIS (pBurst: BSL vs SS)\n');
fprintf('================================================================\n');

% Construct Tall Table

idxGrp = tbl.genotype == 'Control';
idxGrp = true(height(tbl), 1);

% Baseline
varsTbl = {'genotype', 'sbjID', 'frBurst', 'frSingle', 'fr'};
tBsl = tbl(idxGrp, varsTbl);
tBsl.Timepoint = repmat({'BSL'}, height(tBsl), 1);

% Steady State
varsSs = {'genotype', 'sbjID', 'ss_frBurst', 'ss_frSingle', 'ss_fr'};
tSs = tbl(idxGrp, varsSs);
tSs.Properties.VariableNames = varsTbl;
tSs.Timepoint = repmat({'SS'}, height(tSs), 1);

% Concatenate
tblLong = [tBsl; tSs];
tblLong.Timepoint = categorical(tblLong.Timepoint, {'BSL', 'SS'});

% LME Analysis
frml = 'fr ~ (frSingle + frBurst) * Timepoint + (frSingle + frBurst) * genotype + (1|sbjID)';

[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tblLong, frml, ...
    'dist', 'log-normal');


% Per Group
idxGrp = tblLong.genotype == 'MCU-KO';
tblGrp = tblLong(idxGrp, :);

frml = 'fr ~ (frSingle + frBurst) * Timepoint + (1|sbjID)';

[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tblGrp, frml, ...
    'dist', 'log-normal');


%% ========================================================================
%  ENDPOINT MODEL: ACUTE BURST ACTIVITY
%  ========================================================================
% Tests whether burst-spike firing during the acute suppression phase
% predicts total FR recovery, controlling for suppression depth.
%
%   frGain ~ acBGain * genotype + acFrGain + (1|sbjID)
%
% All quantities are log-ratios relative to baseline:
%   frGain   = ln(ss_fr / fr)              Recovery endpoint (BSL -> SS)
%   acBGain  = ln(ac_frBurst / frBurst)    Burst FR change (BSL -> Acute)
%   acFrGain = ln(ac_fr / frTot)           Total FR change (BSL -> Acute)
%
% acFrGain controls for suppression depth. By identity,
% Delta FR_burst = Delta FR + Delta P_burst, so the coefficient on
% acBGain isolates the pattern-specific (burstiness) component:
%   Positive -> induction command (retained bursts drive plasticity)
%   Negative -> error signal (burst loss encodes larger error)


% Load Data
presets = {'steadyState', 'acute'};
[tbl, ~, basepaths, v] = mcu_tblMea('presets', presets, 'flgOtl', true);

% Log-Ratios: BSL -> SS (recovery endpoint)
tbl.frGain = log((tbl.ss_fr) ./ (tbl.fr));
tbl.frGain = log((tbl.ss_fr) ./ (tbl.ac_fr));

% Log-Ratios: BSL -> Acute (suppression state)
tbl.acFrGain = log((tbl.ac_fr) ./ (tbl.frTot));
tbl.acBGain  = log((tbl.ac_frBurst) ./ (tbl.frBurst));
% tbl.acSGain  = log((tbl.ac_frSingle) ./ (tbl.frSingle));

% Endpoint Model
frml = 'frGain ~ acBGain * genotype + acFrGain + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'verbose', true, 'flgStnd', false);

frml = 'frGain ~ (acBGain + pBurst) * genotype + acFrGain + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'verbose', true, 'flgStnd', false);

% Partial Regression (AVP) and LS-Means
hFig = figure;
% hTile = tiledlayout(hFig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% hAx = nexttile;
% [tblRes, hFig] = lme_pr(lmeMdl, 'acBGain', ...
%     'flgMode', 'regression', 'hAx', hAx, ...
%     'varGrp', 'genotype', 'transParams', lmeInfo.transParams);

hAx = nexttile;
[pdRes, hFig] = lme_lsmeans(lmeMdl, {'acBGain', 'genotype'}, ...
    'transParams', lmeInfo.transParams, 'hAx', hAx);


%% ========================================================================
%  BURST EXCESS: PROPORTIONAL ALLOCATION RESIDUAL
%  ========================================================================
% The proportional allocation equation (sGain = alpha + beta * bGain,
% S8 Model 1) defines the expected burst-single coupling for any FR
% change. acBurstExcess is the orthogonal distance from this line in
% (acBGain, acSGain) space — positive means the neuron was burstier
% during acute than proportional allocation predicts. By construction,
% this variable is orthogonal to the proportional allocation confound:
% recovery-driven changes lie ON the line and produce zero excess.

% Exchange Rate Coefficients (pooled across genotypes)
tbl.bGain = log((tbl.ss_frBurst) ./ (tbl.frBurst));
tbl.sGain = log((tbl.ss_frSingle) ./ (tbl.frSingle));

[mdlER, ~, ~, ~] = lme_analyse(tbl, 'sGain ~ bGain + (1|sbjID)', ...
    'dist', 'normal', 'verbose', false, 'flgStnd', false);
alpha = mdlER.Coefficients.Estimate(1);
beta  = mdlER.Coefficients.Estimate(2);

% Burst Excess (orthogonal residual from PA line)
tbl.acSGain = log((tbl.ac_frSingle) ./ (tbl.frSingle));
tbl.acBurstExcess = -(tbl.acSGain - alpha - beta .* tbl.acBGain) / sqrt(1 + beta^2);

% Endpoint Model with Burst Excess
frml = 'frGain ~ acBurstExcess * genotype + acFrGain + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'verbose', true, 'flgStnd', false);

% With Baseline Burstiness
frml = 'frGain ~ (acBurstExcess + pBurst) * genotype + acFrGain + (1|sbjID)';
[lmeMdl, lmeStats, lmeInfo, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'verbose', true, 'flgStnd', false);

% Partial Regression (AVP)
[tblRes, hFig] = lme_pr(lmeMdl, 'acBurstExcess', ...
    'flgMode', 'regression', ...
    'varGrp', 'genotype', 'transParams', lmeInfo.transParams);

hAx = nexttile;
[pdRes, hFig] = lme_lsmeans(lmeMdl, {'acBurstExcess', 'genotype'}, ...
    'transParams', lmeInfo.transParams, 'hAx', hAx);


%% ========================================================================
%  BASELINE GAP UNDER PROPORTIONAL ALLOCATION
%  ========================================================================
% Does the FRH allocation rule (sGain = alpha + beta * bGain, Supp Note 1)
% also govern the between-genotype baseline gap in firing components? If
% MCU-KO units simply populate a higher region of the same ln(FR_single)
% vs ln(FR_burst) curve that FRH recovery traces out, the genotype main
% effect on FR_single dissolves once ln(FR_burst) is conditioned on. That
% would unify baseline compensation and FRH under a single biophysical
% rule - a sharper claim than either section currently makes alone.
%
% IS THIS WORTH DOING?
% Yes, cautiously. The decomposition is free (data already in hand), and
% either outcome is informative: a clean result merges two descriptive
% findings into one mechanism; a failure rules out the simplest unifying
% hypothesis. Frame the reported result as evidence, not proof.
%
% Three caveats worth stating in the manuscript if the test goes through:
%
% (1) beta_FRH is a within-unit slope across time (gain pairs per unit).
%     Mapping it to a between-unit, between-genotype baseline contrast
%     assumes the burst-generation nonlinearity is the same whether
%     driven by FRH plasticity or by constitutive baseline compensation.
%     Theoretically defensible (Supp Note 1 derives the nonlinearity
%     from burst biophysics, not from FRH itself), not guaranteed.
%
% (2) A null for genotype (main or interaction) is failure-to-reject.
%     "Consistent with beta" is weaker than "predicted by beta exactly"
%     and cannot exclude small true offsets below power.
%
% (3) Compositional structure. FR_burst and FR_single partition total
%     firing, so a cross-sectional regression of one on the other at
%     baseline carries arithmetic coupling on top of any biophysical
%     coupling. The test of interest (does genotype add beyond
%     ln(FR_burst)?) remains meaningful under compositionality, but the
%     numerical slope can differ from beta_FRH for reasons unrelated to
%     mechanism. Hence the slope equivalence is a secondary finding; the
%     primary claim is the genotype-term non-significance.
%
% MODEL (MEA baseline, per-unit, animal as random intercept):
%   ln(FR_single) ~ ln(FR_burst) * genotype + (1 | sbjID)
%
%   Target pattern for "baseline gap captured by allocation rule":
%     beta1 (ln(FR_burst) main effect)       ~ beta_FRH
%     beta2 (genotype main effect)           ~ 0 (n.s.)
%     beta3 (genotype x ln(FR_burst))        ~ 0 (n.s.)
%
% ESTIMATOR CHOICE.
% Published FRH beta (Table S8 Model 1, Fig 4H) is OLS-in-LME. Matching
% that estimator here is essential for apples-to-apples slope comparison
% - OLS attenuates slopes toward zero when the predictor has measurement
% noise, so mixing OLS (FRH) and orthogonal (baseline) would bias the
% comparison in a predictable direction. An orthogonal (PCA) slope
% parallel to mea_allocation.m is reported below as a robustness check,
% but it is NOT the primary test.
%
% Animal-level sanity check follows because within-animal unit-level
% clustering inflates effective n. A per-mouse aggregate regression with
% no random effects tells us whether the same qualitative pattern holds
% when each mouse contributes one point.

% Baseline log-transformed FR components (add if missing - tbl may have
% been overwritten by an earlier section with acute preset)
if ~ismember('lnFrBurst', tbl.Properties.VariableNames)
    tbl.lnFrBurst  = log(tbl.frBurst);
    tbl.lnFrSingle = log(tbl.frSingle);
end

% Primary model: OLS-in-LME, matches FRH estimator
frml = 'lnFrSingle ~ lnFrBurst * genotype + (1|sbjID)';
[mdlBsl, statsBsl, infoBsl, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'verbose', true, 'flgStnd', false);

% Guard: auto-log on an already-logged predictor would mangle units.
% lme_analyse re-transforms numeric predictors when skew > skewThr (=2).
if isfield(infoBsl.transParams.varsTrans, 'lnFrBurst') && ...
        ~isempty(infoBsl.transParams.varsTrans.lnFrBurst.logBase)
    warning(['[BSL allocation] lme_analyse re-logged lnFrBurst ', ...
        '(skew exceeded threshold). Slope now in log10-per-ln units; ', ...
        'pass a transTemplate with logBase=[] to disable.']);
end

% Reference: FRH beta (Table S8 Model 1) on the same tbl and estimator
tbl.bGain = log((tbl.ss_frBurst) ./ (tbl.frBurst));
tbl.sGain = log((tbl.ss_frSingle) ./ (tbl.frSingle));
tbl.frGain = log((tbl.ss_fr) ./ (tbl.fr));

frml = 'sGain ~ bGain * genotype + (1|sbjID)';
[mdlFrh, ~, ~, ~] = lme_analyse(tbl, frml, ...
    'dist', 'normal', 'verbose', true, 'flgStnd', false);

% Slope comparison (named lookup, robust to coefficient ordering)
betaBsl = mdlBsl.Coefficients.Estimate(strcmp(mdlBsl.Coefficients.Name, 'lnFrBurst'));
seBsl   = mdlBsl.Coefficients.SE(strcmp(mdlBsl.Coefficients.Name, 'lnFrBurst'));
betaFrh = mdlFrh.Coefficients.Estimate(strcmp(mdlFrh.Coefficients.Name, 'bGain'));
seFrh   = mdlFrh.Coefficients.SE(strcmp(mdlFrh.Coefficients.Name, 'bGain'));

fprintf('\n================================================================\n');
fprintf(' BASELINE ALLOCATION vs FRH BETA\n');
fprintf('================================================================\n');
fprintf('  beta_BSL (baseline)  = %.3f +/- %.3f\n', betaBsl, seBsl);
fprintf('  beta_FRH (recovery)  = %.3f +/- %.3f\n', betaFrh, seFrh);
fprintf('  |Delta|              = %.3f\n', abs(betaBsl - betaFrh));

% Animal-level sanity check: one point per mouse, no random effect
animTbl = groupsummary(tbl, {'sbjID', 'genotype'}, 'mean', ...
    {'lnFrBurst', 'lnFrSingle'});
animTbl.Properties.VariableNames = regexprep( ...
    animTbl.Properties.VariableNames, '^mean_', '');
mdlAnim = fitlm(animTbl, 'lnFrSingle ~ lnFrBurst * genotype');
fprintf('\n--- Animal-level sanity check (one point per mouse) ---\n');
disp(mdlAnim.Coefficients);

% Orthogonal regression robustness check (parallels mea_allocation.m).
% Per-genotype PCA slope in (ln(FR_burst), ln(FR_single)) space.
fprintf('\n--- Orthogonal (PCA) slopes ---\n');
grps = categories(tbl.genotype);
for iGrp = 1:numel(grps)
    idx = (tbl.genotype == grps{iGrp}) & ...
        ~isnan(tbl.lnFrBurst) & ~isnan(tbl.lnFrSingle);
    xy  = [tbl.lnFrBurst(idx), tbl.lnFrSingle(idx)];
    v   = pca(xy);
    slpOrtho = v(2, 1) / v(1, 1);
    fprintf('  %-8s: beta_ortho = %.3f (n = %d units)\n', ...
        char(grps{iGrp}), slpOrtho, sum(idx));
end

% Visual diagnostic
hFig = figure('Position', [100 100 650 550], 'Color', 'w', ...
    'Name', 'Baseline allocation: ln(FRburst) vs ln(FRsingle)');
plot_scat(tbl, 'lnFrBurst', 'lnFrSingle', 'g', 'genotype', ...
    'fitType', 'linear', 'flgStats', true, 'alpha', 0.6);
xlabel('ln(FR_{burst})  [ln Hz]');
ylabel('ln(FR_{single})  [ln Hz]');
title('Baseline component allocation');