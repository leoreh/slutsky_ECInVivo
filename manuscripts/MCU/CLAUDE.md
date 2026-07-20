# MCU Manuscript Pipeline

Manuscript context: MCU-KO vs Control, firing-rate homeostasis (FRH). Primary datasets are multi-electrode arrays (MEA) from dissociated hippocampal cultures undergoing baclofen-induced suppression and in vivo single-unit recordings from freely-behaving mice across chronic baclofen exposure. Both supply the same per-unit firing metrics: mean firing rate, burstiness (pBurst = spike fraction in bursts), burst rate, burst duration, burst inter-interval, burst size, and the component decomposition into burst-spike and single-spike firing rates. The core biological claim is that MCU-mediated mitochondrial Ca2+ signaling is required for upward FRH, dissociable from the basic capacity to modulate excitability (baseline burstiness is elevated in MCU-KO; acute suppression is intact; recovery fails).

## Data layout on disk

Each recording lives in its own basepath (folder). `basepaths2vars` loads `.mat` files by wildcard — one `.mat` per analysis result. Convention: `<basename>.<varname>.mat` (e.g. `<basename>.fr.mat`, `<basename>.stats.mat`, `<basename>.rcv.mat`, `<basename>.ripp.mat`, `<basename>.units.mat`, `<basename>.spikes.cellinfo.mat`, `<basename>.st_metrics.mat`). Returned as a struct array `v(i).<var>` indexed by basepath.

MEA `stats` fields (`br`, `dur`, `freq`, `ibi`, `pBurst`, `bSize`, `frBurst`, `frSingle`, `fr`) are `[nUnits x 3]` matrices — columns are BSL, Acute, SS time windows. `mcu_tblMea` reads them with `idxCol`: 1 for baseline (the un-prefixed columns like `pBurst`, `frBurst`), 2 for acute (prefix `ac_`), 3 for steady-state (prefix `ss_`). `mcu_tblMea` joins multiple subsets (baseline + steady state + acute + temporal dynamics + network) on `sbjID/unitID` when presets request them.

In vivo stats are `[nUnits x nDays]` across BSL, BAC_ON, BAC1, BAC2, BAC3, BAC_OFF, WASH. `mcu_tblVivo` passes `idxCol=[]`, so each table cell holds a vector and `v2tbl` stacks rows across units-by-days. The resulting table has a `day` column (categorical) in addition to `sbjID/unitID/genotype`.

## Loading tables

`mcu_tblMea(presets)` and `mcu_tblVivo(presets)` are the canonical entry points. Presets control which variables get attached.

MEA presets: `steadyState` (ss_* columns from idxCol=3), `acute` (ac_* from idxCol=2), `time` (t_* dynamics aligned to perturbation onset via `mea_tAlign`), `spktimes` (raw spike + burst times), `frNet` (dim, mcc, funcon, funcon_fish from `frNet.corr`), `rcv` (full recovery metrics: rcvBsl, rcvGain, rcvWork, rcvDiff, pertDepth, spkDfct, uPert, uRcv).

In vivo presets: `swv` (waveform metrics), `burst` (burst stats matching MEA fields), `spktimes`, `prc` (phase-response curve), `frNet`, `rippSpks` (per-unit SWR spike metrics — frRipp, frMod, pFire, com, asym, rankMean), `ripp` (per-event SWR properties; note this REPLACES varMap wholesale), `rippMaps` (event-aligned LFP and PETHs), `rippStates` (rate by brain state), `acg` (narrow + wide autocorrelograms), `spkStates` (per-unit metrics split by vigilance state; REPLACES varMap wholesale).

## State-conditioned spiking

Vigilance state is a *row* factor, never a trailing matrix dimension. A metric function takes one interval set `bouts [n x 2]` and returns `[nUnits x 1]` per field; `spk_byCond(fcn, bouts, lbls)` maps it over labelled interval sets and stacks a long table with `uid` and `state` columns. The same call takes time chunks or drug epochs — only `bouts`/`lbls` change.

`spk_states(basepath)` is the per-session producer. It writes `<basename>.stStates.mat` (ACG/ISI metrics from `spktimes_metrics`) and `<basename>.brstStates.mat` (rate/burst metrics from `burst_stats` with `flgPool=true`) as separate files, so each family is computed and loaded independently. Both share the row order (unit × state, state-major), which is why the `spkStates` preset maps them into one varMap without a join; `mcu_tblVivo` asserts the keys agree and rebuilds `unitID` from `uid` so `(1|unitID)` still groups states of the same unit. `ss.info.names` runs one longer than `ss.bouts.times` (trailing `BIN`), and `spk_states` truncates it. Every row carries `nSpks` and `durState` — REM bouts are short, so filter on exposure before reading a state effect.

Segment awareness lives in the metric, not the caller: ISIs never cross a bout, and ACG lags are divided by `nEff(tau)`, the count of spikes with room for a partner at that lag. Restricting spikes and then calling `diff` — what the old `bins` argument did — inflated `cv` by up to 3.8x on a state made of short bouts.

Both loaders accept `basepaths` and `v` (pre-loaded) to avoid re-reading disk. `mcu_tblMea` runs outlier removal by default (`flgOtl=true`): drops non-perturbed units (`~uPert`) and then drops units with `|Pearson residual| > 3` from `ss_fr ~ frBurst + frSingle + (1|sbjID)`, fit per-genotype. `mcu_tblVivo` has `flgClean=true` that drops FS + Other unit types and drops BAC_ON/BAC_OFF/WASH days.

Default basepaths come from `mcu_basepaths(queryStr)` — a string-keyed path registry (e.g. `mcu`, `wt`, `mea_bac`, `mea_mcuko`, `mcu_bsl`, `wt_bsl_ripp`, individual mouse IDs `lh132` etc.). Auto-substitutes `E:\` with `D:\` and applies natsort.

## Config (mcu_cfg)

Genotype labels `{Control, MCU-KO, CAG-MCU-KO}` (`cfg.lbl.grp`), one color row each in `cfg.clr.grp`. Day labels `{BSL, BAC_ON, BAC1, BAC2, BAC3, BAC_OFF, WASH}` (`cfg.lbl.day`). Unit-type labels `{RS, FS, Other}`. Control mice (`cfg.miceWT`): `lh96, lh100, lh107, lh119, lh122, lh123, lh126, lh142`. MCU-KO mice (`cfg.miceMCU`): `lh132, lh133, lh134, lh136, lh137, lh140`. CAG-MCU-KO mice (`cfg.miceCAG`): `raMCU1–raMCU5` — acute viral KO under the CAG promotor, baseline recordings only, one session per mouse. Core vars always loaded: `{fr, units, st_metrics}`.

`mcu_geno(sbjID)` is the single subject → genotype mapper (used by `mcu_tblVivo` and `ripp_screen`). It keeps only genotypes actually present and warns on a subject registered in no cohort — that case used to fall through to Control silently, which mislabelled all 165 viral-KO units as controls.

Genotype labels must not contain `:` or `_`. `lme_postHoc` splits coefficient names on `:` (interaction) and `_` (factor_level), so a label carrying either is rendered as a bogus term — the earlier `CAG:MCU-KO` came out as `(Control vs CAG) * MCU-KO` in the supp tables.

The CAG cohort is opt-in, not a default basepath. It has baseline only, so any `~ genotype * day` model including it is rank deficient and `fitlme` hard-errors. Use `mcu_basepaths('bsl3')` (= `wt_bsl + mcu_bsl + ra`) for three-genotype baseline comparisons; `mcu_tblVivo` defaults stay two-genotype so Table S6 keeps working.

## Analysis stack

`tbl_trans(tbl, ...)` is the transformation primitive. Three modes: FIT (derive params + apply), APPLY (reuse a template on new data without leakage), INVERSE (back-transform to original scale — used by `lme_lsmeans` and `lme_pr` for display). Supports log (natural or numeric base), logit (with Smithson-Verkuilen squeeze for endpoints), additive offset for zeros, z-scoring, and reference-group normalization. Auto-log is skewness-gated: applied only when `skewness > skewThr` (default 2); logit is unconditional when requested. Returned `transParams` is the template.

`lme_analyse(tbl, frml, ...)` is the orchestrator. Sequence: parse formula via `lme_frml2vars` → transform numeric predictors via `tbl_trans` (log10 if skewed by default, z-score unless `flgStnd=false`) → check collinearity with `lme_vif` → if `dist` not supplied, auto-select across `{Normal, Log-Normal, Logit-Normal, Gamma}` by AIC (`lme_compareDists`) and run Park's test (`lme_parkTest`) for variance-mean scaling → apply response transform matching `dist` (log for Log-Normal, logit for Logit-Normal, zero-offset for Gamma) → fit via `lme_fit` → post-hoc via `lme_postHoc` (ANOVA, coefficients with CI, simple + marginal effects, Holm correction by default, Satterthwaite DF for LME, Residual DF for GLME) → optional residual diagnostics via `lme_plotRes`. Returns `lmeMdl`, `lmeStats` (tidy results), `lmeInfo` (transforms + diagnostics), `lmeTbl` (data actually used).

`lme_fit` wraps `fitlme` / `fitglme` with auto link-function selection — standard MLE / REML / REMPL. No errors-in-variables or orthogonal estimator anywhere in the LME pipeline.

Post-fit utilities. `lme_lsmeans(mdl, vars, ...)` evaluates predictions on a grid of covariates and back-transforms to the response scale when `transParams` is passed. `lme_pr(mdl, varX, ...)` produces partial-regression / added-variable plots (`flgMode=regression`) or component-plus-residual plots (`flgMode=residual`), with grouping support. `lme_ablation(tbl, frml, ...)` runs feature ablation with k-fold out-of-sample R2 (`partitionMode={full, split, kfold}`, `nrep` repeats). `lme_mediation(tbl, frml, xVar, mVar, ...)` tests indirect effects via two nested LMEs; `lme_mediationPlot` visualizes the path decomposition.

Tidy export: `lme_mdl2tbls(mdl, stats, info)` packs coefficients/ANOVA/effects into a struct array of Title+Table pairs. `lme_save(sheet, tbls, ...)` writes them with metadata rows to one Excel sheet.

GUI-style visualizations: `guiTbl_bar`, `guiTbl_scatHist`, `guiTbl_xy` launch interactive widgets over any table with categorical/numeric variables (useful for triaging without writing new plot code).

## Orthogonal regression caveat

The proportional-allocation β ≈ 0.4 reported in Figure 4H and Table S8 Model 1 is fit by OLS via `lme_analyse(..., 'sGain ~ bGain * genotype + (1|sbjID)', 'dist', 'normal', ...)`. A separate orthogonal (PCA) β is computed in `mea_allocation.m:70` for the Supp Note 1 C/β mathematical derivation. `mcu_rcvSpace.m` overlays orthogonal fit lines on the gain-space scatter via `plot_scat(..., 'fitType', 'ortho', ...)`. These are three different estimators of the same slope — numerically close because R2 is high, but not the same number. When quoting a β, name the estimator.

## Key scripts

`mcu_lme2xls.m` builds the canonical supp tables S1–S9. One section per table: build tbl-subset → fit LME via `lme_analyse` → wrap with `lme_mdl2tbls` → write via `lme_save`. Final section assembles a TOC with HYPERLINK formulas and calls `mcu_xlsFormat` (COM automation) to style the workbook. Output: `D:\OneDrive - Tel-Aviv University\PhD\Slutsky\Manuscripts\MCU\Results\mcu_suppTbl.xlsx`. `flgPlot=false` gates the exploratory GUIs in each section.

`mcu_regression.m` is the exploratory bench for FRH dynamics: state-space trajectory, partial regression, LS-means, mediation-style decompositions, feature ablation, endpoint and burst-excess models, baseline-compensation analyses. Sections are meant to be run individually. Variable naming is informal — `tbl` is overwritten across sections with different presets, and both `tbl` and `tblMea` appear interchangeably.

`mcu_tblMea.m` / `mcu_tblVivo.m` are the loaders. `mcu_rcvSpace.m` plots recovery state-space scatters. `mea_allocation.m` computes the empirical orthogonal β feeding Supp Note 1. `mcu_utypes.m` handles unit-type classification.

## Conventions and gotchas

Log-ratio convention: MATLAB `log` is natural log. Gains are ln ratios: `bGain = log(ss_frBurst / frBurst)`, `sGain = log(ss_frSingle / frSingle)`, `frGain = log(ss_fr / fr)`. `pBurst` is on [0,1] and is logit-transformed before LME entry via `tbl_trans(..., 'logBase', 'logit')`. Stored as `pBurst_trans` alongside the raw value.

Empty categorical levels. `fitlme` / `fitglme` reject a fixed-effect predictor holding a declared-but-absent level with `design matrix X must be of full column rank` — the usual trigger is subsetting a table (`tbl(tbl.day == 'BAC3', :)`). `lme_analyse` now calls `removecats` on every categorical model variable after truncating the table. Levels that are present are untouched, and grouping variables fit identically with or without empty levels (verified: same logLik, same estimates), so this only removes a failure mode. Tables built outside `lme_analyse` still need their own `removecats`.

Auto-transform risk. `lme_analyse` applies log10 to any numeric predictor with `skewness > 2` (default `skewThr`) and z-scores every continuous predictor unless `flgStnd=false`. If you pre-compute a log-transformed predictor and then pass it, inspect `lmeInfo.transParams.varsTrans.<var>.logBase` to confirm no double-log. Pass a `transTemplate` to override.

Outlier removal in `mcu_tblMea` is model-based (Pearson residuals from a per-genotype FR-decomposition LME), so changing the outlier formula cascades into every downstream result. If a number moves unexpectedly, check `flgOtl`.

Canonical paths. Code: `D:\Code\slutsky_ECInVivo`. Results: `D:\OneDrive - Tel-Aviv University\PhD\Slutsky\Manuscripts\MCU\Results`. Research notes / manuscript: `D:\OneDrive - Tel-Aviv University\Obsi_Vaults\Obsi_Research\`.
