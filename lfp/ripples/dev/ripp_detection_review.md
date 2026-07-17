# Ripple detection — review and decisions

**Goal.** Improve ripple detection in `lfp/ripples/ripp_wrapper.m`: better start/stop/peak
times, fewer false positives and negatives, on 24 h mouse CA1 s.pyr LFP (16 ch, no s.radiatum),
without per-recording hand-tuning that could confound the Control-vs-MCU-KO comparison.

The detector, per method, is screened by the harness in `lfp/ripples/screen/` (`ripp_screen`),
which also reproduces the manuscript group difference under each method and launches a
per-method comparison GUI. This note records what we found and decided; the harness is where it
gets quantified across the cohort.

---

## P1 — The detector keys on wide-band power, which the 1/f slope drags low

**Learned.** Detection thresholds a Teager-energy envelope of 80–250 Hz. Lower frequencies carry
far more power (the 1/f slope), so the envelope is dominated by the low edge — the detector
behaves like a low-frequency detector. It catches ripples riding on a large slow deflection and
misses clean fast ripples that lack one (the screenshot). Measured on two mice: the events'
whitened event-triggered spectrum peaks at ~170 Hz, but the detector's statistic peaks near the
band floor. Narrowing to 120–220 Hz on the best channel lifts multi-unit-spiking co-occurrence
(a label-free quality check) from 58→76% on raMCU1 and 77→84% on raMCU3, while finding *more*
events — fewer misses and fewer false alarms at once.

**Decided.** Detect on the ripple band proper (best single channel, ~120–220 Hz). Channel
averaging across shanks is a minor effect (sites are highly correlated); the band is the lever.

**Remains.** Confirm across the full lh* cohort via the screen, and pick the exact band from the
per-cohort whitened peak.

## P2 — Reported ripple frequency is biased low

**Learned.** `.freq`/`.freqEvent` use the Hilbert instantaneous frequency, which rides the 1/f
slope downward: median ~110 Hz reported vs a true spectral peak ~170 Hz. The bias depends on the
1/f slope, which can differ between genotypes — so the manuscript's ripple-frequency values are
both too low and potentially confounded.

**Decided.** Added `ripp.freqPeak` in `ripp_params` — the whitened event-PSD peak, the honest
frequency. It is additive: `.freq`/`.freqEvent` are unchanged, so stored files and the S5 table
do not shift silently. The screen reports both, so the size of the correction is visible.

**Remains.** Decide whether to switch the manuscript's frequency metric to `.freqPeak`.

## P3 — Normalization confound (both directions)

**Learned.** The threshold is z-scored to the NREM signal *including* the ripples. Two opposite
problems: (a) a mouse with bigger ripples gets a bigger SD and a higher bar, which can *hide* a
genuine group difference in rate (the original worry); (b) mice differing in 1/f slope can show a
*spurious* rate difference even with identical ripples (van Schalkwijk & Helfrich 2026, the "1/f
noise floor").

**Decided.** Address both in the screen. (a) A signal-independent baseline — new `zMet='nremBg'`
in `ripp_sigPrep` estimates the threshold from ripple-free NREM background, so ripple power no
longer inflates the bar; the `narrowAbs` method uses it, testing whether per-animal z was hiding
an effect. (b) A matched 1/f-noise surrogate per mouse (exponent from a robust low-frequency
fit → synthesise noise → run the identical detector → its rate is the noise floor); the group
rate difference is reported with the surrogate rate partialled out. A difference that survives
both is real.

**Remains.** Run the surrogate control across the cohort and read off which normalization is
defensible from the numbers.

## P4 — Detection channel

**Learned.** The 4×4 sites are highly correlated in the ripple band; averaging the tagged
channels vs the single best channel changed little. The band (P1) dominates.

**Decided.** Use the best single channel (most 120–220 Hz power in NREM). Low priority.

**Remains.** Confirm in the screen.

---

## Considered and set aside
Kept on record as future options, not adopted now — each is a larger build for marginal gain once
P1 + P3 are handled:
- Rhythmicity gate (cycle-count on the raw trace / spectral-peak prominence) to reject
  non-oscillatory transients.
- Multi-channel spatial consistency (require coherent power on ≥2 shanks).
- rippl-AI CNN as a second detector (realistic F1 ~0.65, ~= a plain band-pass; adds a Python/TF
  dependency).
- FOOOF as a detector front-end (whitening before detection). The van Schalkwijk paper does not
  endorse this; FOOOF's place is frequency/PSD measurement, not detection.

## FOOOF, briefly
Not for detection. It is the right tool for ripple frequency/PSD, but `lfp/FOOOF/fooof_calc` is a
1–100 Hz delta/theta/gamma tool (never fits above 100 Hz, no knee, hard-bins peaks into
delta/theta/gamma, and mis-reports R² as Pearson r) — so it structurally never saw ripples. A
separate, optional cleanup would make it ripple-capable. The harness uses a self-contained
whitened-peak measure instead, so it carries no Python dependency.

## Loading bug found and fixed (root cause of the "undefined states")
`basepaths2vars` matched files by a loose substring glob (`*sleep_states*.mat`), which also matches
`AccuSleep_states.mat`. On a session holding both, the alphabetical first (`AccuSleep_states.mat`)
won and loaded an `ss` without `.bouts`, silently breaking NREM z-scoring and state labels — why
stored raMCU `ripp.state` came out all-undefined. Fixed at the root: the glob now anchors a leading
dot and prefers the exact `.sleep_states.mat`, so it can no longer grab `AccuSleep_states.mat`
(verified on raMCU1: `ss.bouts` present).
