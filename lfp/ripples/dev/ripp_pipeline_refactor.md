# Ripple / event pipeline refactor — plan

Working plan for converging `ripp_wrapper` and `ripp_screen`, moving QA onto an
`accepted` mask, and cleaning the shared `evt_*` layer. Sequenced so each step is
verified before the next, because the QA change cascades into the manuscript.

## Premise (corrected by the audit)

The `evt_*` family (23 files) is not slop. Almost every file exists to
de-duplicate `ripp_wrapper` and `ed_wrapper`, which share one event-agnostic
spine; deleting them re-introduces the duplication. The honest consolidation is
modest. The real problems the audit found are narrow: one dead function, three
functions with vestigial save parameters, two overlapping channel resolvers, one
control-times bug, and a hardcoded NREM magic index.

## Science decision (validated on data)

The false-positive gate is **MUA gain**, not spectral prominence. Whitened-peak
prominence is largely redundant with detection (every detected event already has
ripple-band power), and on the FP-heavy raMCU session prominence and MUA gain
decouple (Spearman ~0), with high-prominence low-spiking events that only MUA
gain removes. MUA gain is orthogonal to detection, so it is the informative gate.
`.peakProm` is kept as a free per-event property, not the gate.

## Done (safe, verified)

- `ripp_params`: added `.peakProm` (whitened ripple-band peak height above 1/f).
- `evt_ctrlTimes:115` bug fixed (end-crop wrote the start column, making start>end).
- `evt_pethPop` deleted (zero callers; only named in comments).

## Target design

### 1. met-driven config
`ripp_methods()` (rename `current` -> `default`) returns a `met` struct carrying
detection (chMode, passband, detectMet, zMet, thr, limDur) and QA (thrEmg,
gainThr, nremOnly, noiseFloor). `ripp_wrapper` takes `met` (default loaded)
instead of the scattered detection params; keeps win/rippCh/flg* as params so the
existing callers (`mcu_ripples`, `mcu_viralKO`) still work.

### 2. QA as an accepted layer (marks, never removes)
`detect -> QA metrics (emg, spkGain) -> params -> maps -> evt_states (state per
event) -> accepted = evt_qa(inTimes = nremTimes, [emg spkGain], [emg<thrEmg,
gain>gainThr])`. Nothing is dropped.
- Screenable: emg + spkGain are saved per event, so re-deriving `accepted` for a
  new gainThr is one `evt_qa` call, with no re-detection.
- Wake inspectable: all events are on disk; guiPath shows them all; analysis
  filters `accepted` (and state == NREM).
- Manuscript cascade (must verify): `ripp.mat` now holds all events, so
  `mcu_tblVivo` (`ripp`/`rippSpks`/`rippStates`) must filter `accepted`, and the
  rate/density path must reflect `accepted`. S4/S5 will change (NREM-only + gain
  gate); that is intended, but re-run `mcu_lme2xls` and confirm the change is the
  expected one, not a bug.

### 3. Pipeline unification
Factor `ripp_wrapper`'s detect->params->maps->QA core into `ripp_detect(basepath,
met, ...)` (optional injected `rippSig`/`v`), which writes nothing. `ripp_wrapper`
wraps it (spikes/phase/save/plot/NeuroScope/curate). `ripp_screen` calls
`ripp_detect` per met with a signal cache keyed by config; `ripp_screenDetect`
dissolves; `muaGain` folds into the existing `evt_spkGain`.

### 4. Channel resolver (retire channelTags.Ripple)
Merge `evt_pickCh` + `evt_rippCh` + the screen's `bestRippCh` into one
ripple-owned `ripp_pickCh`. Consumer path reads `ripp.info.rippCh`; producer path
(detection) takes an explicit channel, else auto-picks the best NREM ripple-band
channel. Verify channel indexing (0- vs 1-indexed, since `bestRippCh` returns
1-indexed and `evt_loadCh` expects 0-indexed) against a real session before
switching the producer default. Update callers (`ripp_wrapper`, `ripp_sigLoad`,
`ed_wrapper`, `guiPath_presets`); delete `evt_pickCh`, `evt_rippCh`.

### 5. Remaining leaf cleanups
- `ripp_sigPrep`: drop `nremBg` (refuted, unused).
- `ripp_times`: drop the dead peakPower/contPow stats.
- `evt_maps`/`evt_spkPeth`/`evt_spksParams`: strip the vestigial flgSave/basepath
  params (no caller triggers them) and update the call sites; `evt_spkPeth` takes
  `fs` instead of hardcoding 1250.
- `evt_boutTimes`: rewrite the docs; resolve NREM by state name, not `boutTimes{4}`.
- `evt2ns`: simplify the accepted / schema-drift handling now that `accepted` is
  always present.

### 6. ED alignment
`ed_wrapper` switches to `ripp_pickCh` and adapts to any `evt_*` signature change.
ED keeps its own (removal-based) flow; only verify it still runs.

## Status (all complete, each verified)

1. Done: science (`.peakProm` + MUA-gain gate), `evt_ctrlTimes` bug, `evt_pethPop`.
2. Done: `ripp_pickCh` merged resolver; `channelTags.Ripple` retired. NOTE: lh
   sessions whose `.ripp.mat` predate `ripp.info.rippCh` now auto-pick the best
   channel (single, e.g. ch 5) instead of the averaged tag `[5 6 7]` - pass an
   explicit `rippCh` to preserve the shank average.
3. Done: `ripp_methods` (met) + `ripp_detect` core + `ripp_wrapper` QA-`accepted`
   (mark, not remove) + `evt_states` accepted-rate + `mcu_tblVivo` accepted
   filter. Verified on lh100: 7601 events, 1577 accepted (all NREM), spks
   aligned, `mcu_tblVivo('ripp')` = 1577 gated rows. S4/S5 WILL change (NREM-only
   + MUA gate) - re-run `mcu_lme2xls` to compare old-vs-new.
4. Done: `ripp_screen` rebuilt on `ripp_detect` (signal cached across methods);
   `ripp_screenDetect` + `ripp_screenMethods` deleted; `muaGain` -> `evt_spkGain`.
5. Done: `nremBg` dropped; `evt_maps` save scaffolding stripped; `evt_boutTimes`
   state indices named + documented; `evt_ctrlTimes` empty-event guard + plot
   lint. ED verified (empty-guard, normal path with events, channel). Remaining
   minor: `evt_spkPeth`/`evt_spksParams` save scaffolding; `ripp_times` dead
   `contPowAvg`.
6. Done: `evt_doc` + memory updated.

## Verified (final run, all pass)

- checkcode clean on all changed files.
- `ripp_wrapper` fresh on lh100: 7601 events, 1577 accepted, spikes aligned.
- `mcu_tblVivo('ripp')`: 1577 gated rows, amp/freq/dur + frac/asym/com populated.
- `ed_wrapper`: empty-event guard (lh100 0 EDs) + normal path (raMCU1 6 EDs).
- `ripp_screen`: default vs fooof, genotype x method table.

## Addendum (260719b) — spkGain margin + interactive gate GUI

- `evt_spkGain` gained a `'margin'` param (default 5 ms): each event window is
  widened by +/-margin before the in-event rate and the internal matched
  `evt_ctrlTimes`, because ripple-associated MUA leads the detected start and
  lags the detected end (the band-pass envelope clips these shoulders). Empirical
  basis (lh100, 1123 NREM events): pooled-MUA baseline 86 Hz, but the edge-aligned
  rate is ~2.7x baseline within 5 ms of each edge; the broader elevation out to
  40 ms is the sharp-wave / clustering floor and is deliberately excluded. Effect
  on the gate is modest — 0->5 ms shifts accepted 938->941 (10 gained, 7 lost);
  the mean gain drops slightly (5.43->5.27) from shoulder dilution; at 20 ms the
  accepted count turns over, so the margin is kept small. `ripp_detect` picks up
  the default automatically, so `.spkGain`/`accepted` now include the shoulder.
- `ripp_gateFig` (static, buggy rho panel, re-ran the whole detect chain,
  unclear counts) replaced by `ripp_gateGui`: an interactive gate built on
  `guiTbl_xy` (kept vs removed mean +/- spread), live MUA-gain / EMG / prominence
  threshold fields, and a running kept/removed count. Takes the in-memory
  `[ripp, aux]` from `ripp_detect` (the saved `.rippMaps.mat` holds accepted
  events only, so a kept-vs-removed view needs the all-events `aux.rippMaps`),
  else detects once. Writes nothing. Verified: checkcode clean, headless build +
  threshold change on lh100. `.gainThr` in `ripp_methods` now points here.
- `.peakProm` (whitened ripple-band peak height above the 1/f floor, in
  `ripp_params`) is a reported property and a GUI knob, not part of `accepted`.
