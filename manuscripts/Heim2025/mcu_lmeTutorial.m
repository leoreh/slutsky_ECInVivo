% mcu_lmeTutorial

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FR analysis during baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grps = {'wt_bsl'; 'mcu_bsl'};

clear grppaths
for igrp = 1 : length(grps)
    grppaths{igrp} = string(mcu_sessions(grps{igrp})');
end


% FR per unit, WT vs MCU across states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frml = 'FR ~ Group * State + (1|Mouse)';

% organize for lme
[lme_tbl, lme_cfg] = mcu_lmeOrg(grppaths, frml, false);

% run lme
iunit = categorical({'pPYR'});
plot_tbl = lme_tbl(lme_tbl.UnitType == iunit, :);
lme = fitlme(plot_tbl, lme_cfg.frml, 'FitMethod', 'REML');

exlTbl = lme2exl(lme);

% parameter definitions:
% FIXED
% Intercept                     – Baseline FR of WT during AW.
% Group_MCU-KO                  – Difference in FR between MCU-KO and WT during AW.
% State_NREM                    – Difference in FR between NREM and AW in WT.
% State_REM                     – Difference in FR between REM and AW in WT.
% Group_MCU-KO:State_NREM       – Interaction: Difference in FR change from AW to NREM between MCU-KO and WT.
% Group_MCU-KO:State_REM        – Interaction: Difference in FR change from AW to REM between MCU-KO and WT.
% RANDOM
% (Intercept) | Mouse           – Variability in baseline FR across mice.
% State_NREM | Mouse            – Variability in FR change from AW to NREM across mice.
% State_REM | Mouse             – Variability in FR change from AW to REM across mice.
% Corr(State_NREM, Intercept)   – Correlation between baseline FR and FR change from AW to NREM across mice.
% Corr(State_REM, Intercept)    – Correlation between baseline FR and FR change from AW to REM across mice.
% Corr(State_REM, State_NREM)   – Correlation between FR changes from AW to NREM and AW to REM across mice.
% Residual (Error term)         – Unexplained variability in FR.

% FitMethod
% Maximum Likelihood (ML) and Restricted Maximum Likelihood (REML) are two
% methods used to estimate parameters in Linear Mixed-Effects (LME) models,
% particularly the variance components of random effects. While both
% approaches rely on likelihood-based estimation, they differ in how they
% handle variance estimation and are used for different purposes.
% 
% Maximum Likelihood (ML) estimates all parameters—both fixed effects and
% random effects—simultaneously by maximizing the likelihood of the
% observed data. It considers the total variability in the data and
% attributes it to both fixed effects (predictors) and random effects
% (individual variability). Since ML estimates the fixed effect
% coefficients and variance components in one step, it allows for direct
% model comparisons using likelihood ratio tests (LRT), AIC, or BIC, making
% it the preferred method when testing different fixed-effects models.
% However, ML tends to underestimate the variance of random effects,
% particularly in small sample sizes, because it does not account for the
% degrees of freedom lost when estimating fixed effects.
% 
% Restricted Maximum Likelihood (REML), also called Restricted ML or RML,
% adjusts for this bias by separately estimating the variance components of
% random effects before estimating fixed effects. It does this by first
% removing the contribution of the fixed effects from the likelihood
% calculation and then estimating the remaining variance, ensuring that the
% estimation of random effects is unbiased. This makes REML more reliable
% for estimating variance components, particularly when the number of
% groups (e.g., subjects, animals) is small. Unlike ML, REML cannot be used
% for likelihood ratio tests of models with different fixed effects
% structures, since it adjusts for fixed effects before computing
% likelihood.
% 
% The choice between ML and REML depends on the goal of the analysis. If
% the focus is on comparing models with different fixed effects (e.g.,
% testing whether adding an interaction term improves model fit), ML should
% be used. Once the fixed effects structure is finalized, REML is preferred
% for estimating the final model because it provides more accurate variance
% estimates for random effects. In cases where random effect variance is
% small or near zero, as seen in your model, REML can sometimes provide a
% more stable estimate and reduce computational issues such as convergence
% failures.





mean(plot_tbl.FR(plot_tbl.Group == 'MCU-KO' & plot_tbl.State == 'AW'))


fh = mcu_lmePlot(plot_tbl, lme, 'ptype', 'line');


% plot
fh = mcu_lmePlot(plot_tbl, lme, 'ptype', 'bar');

% save
frml = [char(lme.Formula), '_', char(iunit)];
frml = frml2char(frml, 'rm_rnd', false);
th = get(gcf, 'Children');
title(th, frml, 'interpreter', 'none')
axh = get(th, 'Children');
ylim(axh(3), [0 3])
grph_save('fh', fh, 'fname', frml, 'frmt', {'ai', 'jpg'})

prism_data = fh2prism(fh);
