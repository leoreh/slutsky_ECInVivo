# ED pipeline rebuild — record and findings

Rebuilt `lfp/ed` to the staged shape the ripple pipeline reached in
`../../ripples/dev/ripp_pipeline_refactor.md`, and calibrated the detector on the
only labelled ED data in the repo. This file is the record of what was measured,
what was decided, and what remains open.

## What the old pipeline was

One pass: `ed_wrapper` loaded a signal, thresholded a moving mean/SD z-score of
the RAW trace (`ed_detect`), measured widths (`ed_params`), ran `evt_qa` +
`evt_subset` to DELETE the failures, then ran matched controls, state labels,
maps and a full spike analysis (`evt_spks`, PETH, rasters).

An audit (5 readers, each finding adversarially verified; 48 findings confirmed)
found these defects in the detection code, all now gone with the rewrite:

- **Self-inflating baseline.** `movmean`/`movstd` over a 5 s window that contains
  the event, so each discharge raised its own threshold and a burst suppressed
  its later members.
- **Dimensionally wrong amplitude gate.** `thrAmp = thr*sigma + mu` is a signal
  LEVEL, compared against `peak2peak`, a DIFFERENCE. The gate's stringency
  tracked the channel's DC offset; with mu ≈ 0 it was near-vacuous.
- **Inert twin-peak merge.** `min(sig(a:b)) > 0.2` on a trace whose sigma is
  ~200, with no distance cap, and unreachable for negative-going discharges —
  which is the polarity these discharges actually have.
- **Merge chain.** Pairwise crossing deletion + `uniquetol` picked survivors by
  sample order, not by peak magnitude, and could delete an event because of an
  unrelated neighbour.
- **NaN widths pass QA.** An event wider than the ±marg clip returned a NaN
  width, which `evt_qa` treated as PASS — so a `durLim` gate silently exempted
  the longest events.
- **`ed_sigLoad` fs.** The 'lfp' branch overwrote `fs` with the LFP rate and
  then windowed the EMG with it, scoring events against a misaligned trace.
- **Stale GUI token.** `guiPath(..., 'preset', 'EDs')` matched no preset; it
  landed on `preset_ed` only through the auto-detect fallback.

## What it is now

`ed_detect` → `ed_curate`, plus `guiPath` for the per-event pass. Same contract
as ripples: detect writes nothing and seeds `.accepted` all-true, the gate is a
separate stage over the saved struct, curation MARKS and never removes, and a
save backs the file up first.

| file | lines | role |
|---|---|---|
| `ed_methods` | 67 | the met config: every detection knob + the default `.qa` spec |
| `ed_sigLoad` | 148 | detection signal + EMG, windowed (`chMode` 'eeg' \| 'ripp') |
| `ed_detect` | 188 | stage 1: signal prep, candidates, every per-event feature |
| `ed_params` | 143 | the features, measured against an event-excluded ring |
| `ed_curate` | 307 | stage 2: the gate, headless or in a bulk-threshold GUI |
| `ed_wrapper` | 179 | the thin chain |
| `ed_tbl` | 108 | counts and rates per session × state — the cross-session product |

`ripp_gate` moved to `lfp/events/evt_gate.m` (a pure rename; the body was already
event-agnostic) and is now shared. `evt_qa` was deleted — `ed_wrapper` was its
only caller. `evt_subset` stays; `mcu_tblVivo` uses it on ripples. `evt_states`
gained a `basename` parameter, without which every session in a flat folder
(the EA cohort) overwrote one `<folder>.edStates.mat`.

No analyze stage, and no spike products. The ED question is counts by state.

## The labelled data

`D:\Data\EA\*.ied.mat` hold legacy `IED.data` objects whose `last_mark` records
how far a curator actually got. Only the prefix is real:

| session | detected | reviewed | accepted | rejected |
|---|---|---|---|---|
| 220611_0750 | 3592 | 14 | 4 | 10 |
| 220615_0801 | 291 | 26 | 24 | 2 |
| 220824_0906 | 7718 | 151 | 1 | 150 |

**The stored `accepted` flag is NOT ground truth** — beyond `last_mark` it is the
constructor default (true). An early calibration run used all 253 "accepted"
events of 220615 as a recall target; that was measuring agreement with the old
detector, not recall of discharges. The 191 reviewed events are the real labels.

`sSig.eeg` is byte-identical to the `LFP` source inside the `.ied` objects, so
the legacy detections map onto the current signal exactly.

## What separates a discharge from a false positive

Confirmed discharges are isolated sharp NEGATIVE deflections (2000–4000 units on
a ~215-unit background) with a following slow wave. The rejected events are
broadband noise bursts and single-sample glitches.

Seven candidate features were measured on the 191 labelled events. Pooled AUCs
looked good (background roughness 0.93, hfRatio 0.92 inverted) but **within a
session none of them separate** — 220611's AUC is 0.50, and the other two
sessions have degenerate class balance (24/26 and 1/151). So the labelled set
demonstrates session-level discrimination only: it distinguishes a session whose
detections are artifacts from one whose detections are discharges. Any claim of
a validated per-event classifier here would be unsupported.

Two features carry that separation and are kept as the gate:

- **`.ampZ`** — peak of the band-passed trace over the robust scale of an
  event-excluded ring (±25 ms cut out, out to ±1 s). Median 7.7 on accepted vs
  4.3 on rejected.
- **`.hfRatio`** — peak in 150–500 Hz over peak in the discharge band. Median
  1.41 accepted vs 5.17 rejected. Orthogonal to amplitude, which is why it adds
  anything: it measures energy faster than any discharge can be.

`ampZ ≥ 6` and `hfRatio ≤ 3` together keep 28/29 accepted and 19/162 rejected.

## The simplification pass

A second pass measured whether each piece of machinery earns its place. Three
did not, and are gone:

- **The block-wise baseline.** `ed_sigPrep` normalised the detection signal by a
  robust scale estimated per 10 s block and interpolated — 45 lines. Swapping it
  for ONE robust scale over the whole recording changes the final counts by a
  few percent and the control-to-epileptic contrast not at all (15.3× vs 15.8×
  on a matched gate). It is gone, and `ed_sigPrep` folded into `ed_detect` as a
  20-line local, since nothing else called it.
- **`.pol`.** Polarity is now the sign of `.amp`, which was unsigned for no
  reason. One field instead of two.
- **The state filter in curation.** ED asks how discharges distribute over
  states, so restricting them is a question for `ed_tbl`, not a curation
  decision. Removing it took the checkbox panel, `qa.states`, `qa.unscored` and
  three helper functions with it. The waveform view still tiles by state.

Also dropped as unrequested: `flgNS`, `flgDetectOnly`, `mapDur`, `flgAll`, and
the `baseWin` / `coreWin` / `ringWin` knobs (now named constants where used).

One simplification was tried and REVERTED on evidence. Making `.ampG` the
band-passed peak over the recording scale would have made it identical to the
detection statistic — elegant, since the detection threshold would then be a
provable lower bound on the gate. It costs too much: contrast falls 27.9× →
15.6×. Most of a discharge's amplitude sits below the 10 Hz corner, so the raw
trace measures something the filtered peak throws away. `.ampG` stays raw.

Net: 8 files → 7, ~1450 lines → 1140, and the contrast improved slightly
(26.3× → 27.9×) because the detection net widened.

## Choosing the band

Measured, not assumed — ring-z of the confirmed discharges on 220615 against
4000 random times:

| band (Hz) | median at events | p99 background | ratio |
|---|---|---|---|
| 5–40 | 7.18 | 4.27 | 1.68 |
| 20–80 | 6.82 | 5.11 | 1.34 |
| **10–100** | **10.13** | **5.22** | **1.94** |
| 25–150 | 10.89 | 5.72 | 1.91 |
| 30–200 | 13.39 | 6.44 | 2.08 |

10–100 Hz was taken: near-best separation, and unlike 30–200 it stays clear of
the band `hfRatio` uses to measure contamination. Smoothing the envelope before
thresholding was dropped — a window wide enough to steady a crossing also halves
a sharp peak, costing about a third of the confirmed discharges.

## What detection cannot do

**No single amplitude or sharpness threshold gets to ~100 events per session
without losing most of the discharges.** Across four bands × seven thresholds:

- 220615 at 10–100 Hz, thr 6: 5053 candidates, 24/24 confirmed discharges kept.
- The same session at thr 12: 62 candidates, 4/24 kept.

The confirmed discharges sit at 6–12× the robust background; over 22 h at
1250 Hz the background crosses that level thousands of times on its own. A
global (recording-wide) normalisation instead of a local one was tried and does
not help — at the threshold that keeps the discharges, the artifact session
still yields 10 763 candidates.

This is why the operating point belongs in `ed_curate` with the waveforms in
view, and not in a constant.

## The shipped defaults, and their evidence

Detection: `chMode 'eeg'`, `passband [10 100]`, `thr 6`, `limDur [4 200 40] ms`.

QA: `ampG ≥ 8`, `ampZ ≥ 7`, `hfRatio ≤ 2`, `emg ≤ 2`.

Every gate criterion was ablated against the control-to-epileptic rate contrast.
Dropping `ampG` costs the most (15.3× → 9.1×), then `ampZ` (→ 12.0×), then `emg`
(→ 7.8×); `hfRatio` moves no rate but halves the survivors on the
artifact-dominated session and removes 6 of the 162 curator-rejected events. All
four stay.

`.ampG` (raw amplitude in recording-wide robust units) is not redundant with
`.ampZ`: the first asks whether a deflection is large for the recording, the
second whether it stands out from its own neighbourhood, and they disagree on a
modest deflection inside a quiet stretch.

`.emg ≤ 2` came from a second, independent comparison — four non-epileptic mice
against the three epileptic sessions, on the reasoning that a control should
have almost no discharges:

| spec | median control /h | median epileptic /h | ratio |
|---|---|---|---|
| no EMG bound | 2.84 | 51.72 | 18.2× |
| **emg ≤ 2** | **1.97** | **51.68** | **26.3×** |
| emg ≤ 1 | 0.29 | 1.73 | 6.1× |

`emg ≤ 1` is where the WAKE/NREM boundary sits (`evt_emgScore`), and it destroys
the wake arm — 220615 falls from 1765 events to 39. At 2 the bound costs 0.8%
and 0.1% on the two discharge-rich sessions and removes 62% of the survivors on
the artifact-dominated one.

## Verified

`checkcode` clean on every touched file. `evt_gate` reproduces the stored ripple
curation bit for bit (read-only check; no ripple product was rewritten). The
`ed_curate` and `guiPath` GUIs build headless, on both a scored session and an
unscored one. ~12 s per 24 h session.

Final counts, whole chain, default met:

| session | | candidates | accepted | rate |
|---|---|---|---|---|
| lh100 | control | 1271 | 39 | 1.7/h |
| lh107 | control | 787 | 50 | 2.1/h |
| lh119 | control | 1021 | 18 | 0.8/h |
| lh132 | control | 2452 | 252 | 10.8/h |
| 220611_0750 | epileptic | 6605 | 1248 | 53.0/h |
| 220615_0801 | epileptic | 3362 | 1531 | 68.1/h |
| 220824_0906 | epileptic | 10721 | 83 | 3.5/h |

Contrast 27.9×. The three genuine controls land under 100 events per session,
which is the scale the pipeline was asked for.

## Open

- **lh132 is an outlier control**: 2703 candidates → 257 accepted (11/h), an
  epileptic-range rate, and the EMG bound barely touches it (253). Either real
  epileptiform activity in an MCU-KO mouse or a bad channel. Worth looking at
  before it enters any group comparison.
- **Unscored events.** On lh100, 32 of the 74 accepted events fall in no scored
  bout, so the state rows of `ed_tbl` sum to 42 while the `ALL` row reads 74.
  The gap is a direct read on scoring completeness. With `qa.states = []` the
  `unscored` flag is a no-op (an empty state list already keeps everything); it
  only bites once states are restricted in the GUI.
- **`chMode` default.** 'eeg' was chosen because it is the only source present
  on every scored session, including EA, which has no `session.mat` and no
  binary `.lfp`. 'ripp' puts EDs and ripples on one electrode and may be the
  better choice for a within-mouse comparison.
- **No per-event validation.** The gate is an artifact-bulk filter. Confirming
  it would need a curator to label a few hundred events *within* one session,
  spanning the accept/reject boundary — which the current labelled prefix does
  not do.
