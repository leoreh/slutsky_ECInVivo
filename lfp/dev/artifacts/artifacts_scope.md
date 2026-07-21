# Where lfp_artifacts belongs, and where it must not go

Measured 260721 on the five raMCU sessions (the only ones with EDs). Figures in
this folder show the detected waveforms; the numbers below decide the scope.

## 1. Ripples: the baseline fix is the whole value

Fresh detection with the current `ripp_methods('default')`, restricted to NREM,
against the padded artifact mask on the ripple channel:

| session | NREM events | in artifact | of those, pass emg<1 | net leak |
|---|---|---|---|---|
| raMCU1 | 29835 | 0.3% | 96.5% | 0.27% |
| raMCU2 | 18758 | 0.0% | 100% | 0.02% |
| raMCU3 | 20625 | 0.1% | 94.7% | 0.09% |
| raMCU4 | 21387 | 0.4% | 13.4% | 0.05% |
| raMCU5 | 29498 | 2.3% | 12.4% | 0.29% |

Stored accepted sets are equally clean (0.2-0.8% co-occurrence). The mask is
baseline-only by design, so artifacts stay detectable - but almost none of them
survive into an accepted NREM event, and the residue is under 0.3% everywhere.

**No per-event artifact metric is needed in the ripple gate.**

## 2. ED: the detector flags the discharges, not the artifacts

`lfp_artifacts` thresholds absolute deviation from the median. A discharge is a
large sharp deflection; a movement step is a large sharp deflection. On the ED
channel they occupy the same amplitude range, so the rule cannot separate them.

| session | accepted EDs | thr (uV) | median abs(amp) | amp/thr | % amp>thr | % flagged |
|---|---|---|---|---|---|---|
| raMCU1 | 25 | 3060 | 3118 | 1.02 | 52.0 | 48.0 |
| raMCU2 | 147 | 1459 | 5771 | 3.96 | 84.4 | 85.7 |
| raMCU3 | 292 | 1173 | 4211 | 3.59 | 80.1 | 80.1 |
| raMCU4 | 75 | 1435 | 1607 | 1.12 | 66.7 | 64.0 |
| raMCU5 | 490 | 1340 | 3928 | 2.93 | 99.0 | 99.6 |

`% amp>thr` tracks `% flagged` session by session, and **98.4-100% of the
flagged discharges exceed the threshold on their own amplitude** - the flag is
the discharge itself, not an artifact it happens to sit inside.

No choice of `thrFactor` fixes this. raMCU2 would need >4x the current
threshold to clear its discharges; raMCU1 would need >1.02x, i.e. essentially
the same threshold that is meant to catch its artifacts.

Two consequences:

- **Do not gate EDs on `lfp_artifacts`.** It would delete 48-99.6% of the
  curated discharges.
- **The ED baseline does not need it either.** `ed_detect/sigPrep` scales by
  median + MAD over the recording, whose 50% breakdown point makes ~1%
  contamination irrelevant. The ripple pipeline was vulnerable because it used
  mean/SD of a quantity quartic in voltage, not because artifacts are large.

What separates a step from a discharge is the slow recovery (the +/-500 ms RC
tail in `raMCU5_..._artifacts.png`), not the size. `lfp_artifacts` deliberately
dropped its derivative criterion in favour of size alone - correct for a
baseline estimate, useless for this discrimination. `ed_clust` already keys on
the +/-50 ms waveform, which is where the distinction lives.

## 3. State labels: the same confound, one level up

Whole-recording epoch coverage (ED channel, thrFactor 8):

| session | all | WAKE | NREM |
|---|---|---|---|
| raMCU1 | 0.34% | 0.05% | 0.72% |
| raMCU2 | 1.08% | 1.64% | 0.77% |
| raMCU3 | 5.26% | 12.84% | 0.26% |
| raMCU4 | 0.92% | 1.73% | 0.19% |
| raMCU5 | 4.23% | 8.26% | 1.53% |

Three obstacles to injecting these as BIN by default:

1. **It would delete epileptiform activity from the spectrogram**, in exactly
   the mice that have it (section 2). A bias, not a cleanup.
2. **Bout architecture is a result, not a by-product.** raMCU5 NREM bouts go
   212 -> 405. `bouts.boutLen` is reported (`mcu_states_McuVsWt.m`,
   `mcu_lfp_wrapper.m`) and bout `Duration` is a regressor in the ripple
   density LME (`evt_states` -> `mcu_tblVivo`).
3. **The mask is channel-dependent.** raMCU3 NREM coverage is 0.058% on the
   ripple channel and 0.26% on the ED channel; WAKE reaches 12.84%. A
   state-level mask must be computed on `sSig.eeg` - the trace AccuSleep
   actually scored - or the mask and the labels describe different signals.

## 4. Channel picking: no change, measured

Both producer picks are amplitude statistics over a short probe, so both look
vulnerable. Reproducing each probe and scoring every channel with and without
its own artifact samples:

| session | ED probe flagged | ED pick raw / masked | RIPP probe flagged | RIPP pick raw / masked |
|---|---|---|---|---|
| raMCU1 | 0.0% | 5 / 5 | 0.0% | 5 / 5 |
| raMCU2 | 0.0% | 5 / 5 | 0.0% | 5 / 5 |
| raMCU3 | 1.2% | 14 / 14 | 0.1% | 15 / 15 |
| raMCU4 | 0.0% | 7 / 7 | 0.0% | 7 / 7 |
| raMCU5 | 0.3% | 12 / 12 | 0.0% | 8 / 8 |

**The pick never moves.** Artifacts are common-mode, so they inflate every
channel's score by a similar factor and the ranking survives. The ripple probe
is one NREM bout, the state where artifacts are rarest (0.0-0.1% of samples).

`ed_pickCh` step 1 - the Tukey crossing-rate drop that removes a broken channel
\- is also unchanged: raMCU1 drops [5 7] and raMCU5 drops 12, masked or not.
The prediction that masking would hide the broken channel was wrong: its excess
crossings sit at 8x the 60-150 Hz band scale, hundreds of uV, far below the
1.2-3 mV raw threshold, so the mask barely touches them.

(The step-2-only scores above pick 5 and 12 for raMCU1/raMCU5 where the shipped
function picks 6 and 3. That is step 1 doing its job - both raw picks are
exactly the channels it drops - and confirms the reproduction is faithful.)

## The rule that decides all of these

The mask is safe wherever the quantity of interest sits far below the
threshold, and unsafe wherever it does not:

| use | quantity | threshold | safe? |
|---|---|---|---|
| ripple baseline | ripples, ~100 uV | 0.9-3 mV | yes |
| ripple / ED channel pick | band RMS, percentile | same | yes, but changes nothing |
| ED gate, ED channel pick step 2 | discharges, 1.6-5.8 mV | 1.2-3.1 mV | no |
| state labels in ED mice | discharges | same | no |
