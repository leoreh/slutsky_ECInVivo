# ED pipeline — rebuild record and findings

`lfp/ed` rebuilt to the staged shape of the ripple pipeline
(`../../ripples/dev/ripp_pipeline_refactor.md`), then rebuilt a second time when
the first detector turned out to be detecting the wrong thing. This file records
what was measured, what was decided, and what is still open.

## The correction that matters most

The first version of this detector was calibrated on the EA cohort using the
legacy `.ied.mat` "accepted" flags, and validated by the fact that it returned
more events in EA than in `lh*` controls. **Both were wrong.**

- The legacy `accepted` flag is only real up to `last_mark`; beyond it it is the
  constructor default `true`. Using the full set as a target measured agreement
  with the old detector, not recall of discharges.
- "More than controls" is not validation. Only the **CAG-MCU-KO cohort
  (`raMCU*`)** carries epileptiform discharges, and only a few dozen per 24 h.
  Control and MCU-KO mice have essentially none. A detector reporting tens of
  EDs in an `lh*` mouse is broken, not interesting.
- Deriving the detector from first principles skipped two sources available the
  whole time: the legacy `lfp/+IED` implementation in this repo (over-built, but
  its defaults encode real choices — positive-going threshold, spike width, peak
  frequency) and the published IED-detection literature.

The ground truth used from here on is 47 discharges curated by hand in raMCU3,
raMCU4 and raMCU5, persisted alongside this file in `ed_curatedTimes.mat`.

## What a discharge actually looks like

A sharp **biphasic** complex ~20–30 ms long: a fast negative notch immediately
followed by a large **positive** peak, sitting on quiet, flat background, and
decaying back to baseline within ~50–100 ms. Amplitude 400–2000 units against a
background of ~100–500.

## The four things that are not discharges

Each defeats the test that catches the others, which is why the gate needs three
criteria and why three are enough.

| confusion | what it looks like | what catches it |
|---|---|---|
| sharp waves / slow LFP deflections | slow monophasic **negative** excursion, trough ~40 ms out, decaying over 300+ ms | `fastZ` (sharpness) |
| step artifacts | trace drops and never returns | `posZ` (positive excursion over a **pre**-event baseline) |
| noisy / oscillatory epochs | the "event" is the biggest of many comparable peaks | `isoZ` (isolation from a local event-excluded ring) |
| saturation artifacts (lh100) | amplitudes 10–25k, 10× any real discharge | an upper `posZ` bound was tried and removed nothing; left out |

## The detection band is the whole game

Peak amplitude in a band, over that band's own background, measured at the
curated discharges versus the events the curator rejected:

| band (Hz) | AUC | median at discharges | median at rejects |
|---|---|---|---|
| 2–25 (slow) | **0.54** | 7.7 / 7.2 / 9.0 | 7.0 / 7.9 / 8.6 |
| 25–80 | 0.958 | 11.8 / 18.3 / 20.1 | 3.7 / 6.3 / 6.2 |
| 60–250 | 0.985 | 37 / 52 / 61 | 10 / 7 / 6 |
| **60–150** | **0.989** | 32 / 34 / 47 | 8 / 6 / 6 |

Below ~30 Hz a discharge and an ordinary sharp wave are **not separable at
all**. The first detector ran at 10–100 Hz, dominated by that useless slow band,
so it proposed thousands of normal deflections — and in `ed_curate` the
kept-vs-removed mean was then computed over a set that was mostly ordinary LFP,
burying the discharge shape the view exists to show.

60–150 was taken over 60–250: statistically indistinguishable, but it sits above
50 Hz mains and clear of most of the ripple band.

Other measurements from the same labelled set: half-width 7–9 ms at discharges
vs 15–22 ms at rejects; polarity positive in 66% of discharges vs 17% of
rejects — informative but not decisive, because the complex is biphasic and the
largest excursion is sometimes the negative notch. Polarity is therefore
reported as the sign of `.amp` and never gated.

## Shipped defaults

Detection: `chMode 'eeg'`, `passband [60 150]`, `thr 8`, `limDur [4 200 40] ms`,
polarity-blind.

Gate: `fastZ ≥ 15`, `posZ ≥ 5`, `isoZ ≥ 20` — each near the 5th percentile of the
curated discharges (18, 8, 17), keeping 43 of 47.

**The gate produces the set worth REVIEWING, not the answer.** It leaves 15–80
events on a CAG mouse — the ~100 that gets curated down to ~20 — and single
digits on a control. The per-event pass in `guiPath` is the final word.

`emg` is computed and is a knob, but is not bounded: once the three criteria are
applied it removed nothing measurable.

## Architecture

`ed_detect` → `ed_curate` → `guiPath`, plus `ed_tbl`. Detect writes nothing and
seeds `.accepted` all-true; the gate is a separate stage over the saved struct;
curation MARKS and never removes; a save backs the file up first.

| file | role |
|---|---|
| `ed_methods` | the met config: every knob + the default `.qa` spec |
| `ed_sigLoad` | detection signal + EMG, windowed (`chMode` 'eeg' \| 'ripp') |
| `ed_detect` | stage 1: signal prep, candidates, every per-event feature |
| `ed_params` | the features |
| `ed_curate` | stage 2: the gate, headless or a bulk-threshold GUI |
| `ed_wrapper` | the thin chain |
| `ed_tbl` | counts and rates per session × state |

`ripp_gate` moved to `lfp/events/evt_gate.m` and is shared with the ripple
pipeline. `evt_qa` was deleted (`ed_wrapper` was its only caller). `evt_states`
gained a `basename` parameter, without which every session in a flat folder
overwrote one `<folder>.edStates.mat`.

## Defects fixed in the original ED code

An audit (5 readers, each finding adversarially verified; 48 confirmed) found: a
moving mean/SD baseline that each event inflated itself; an amplitude gate
comparing a signal LEVEL against a peak-to-peak DIFFERENCE, so its stringency
tracked the channel's DC offset; an inert twin-peak merge (`min(sig) > 0.2` on a
trace with sigma ~200, no distance cap, unreachable for negative-going events); a
merge chain picking survivors by sample order rather than magnitude; NaN widths
passing QA as valid; an fs mix-up scoring EMG at the wrong rate; and a stale
`'EDs'` preset token. All are gone with the rewrite.

## Simplification pass

Measured, then removed as earning nothing: a block-wise robust baseline (45
lines; changed no result versus one recording-wide scale), `.pol` (now the sign
of `.amp`), and the state filter in curation (ED asks how discharges distribute
over states, so restricting them is a question for `ed_tbl`). Also dropped as
unrequested: `flgNS`, `flgDetectOnly`, `mapDur`, `flgAll`.

## The detection signal is not comparable across mice

`met.chMode 'eeg'` reads `sleep_sig.eeg`, which is **an average of whichever
channels were picked for sleep scoring** — a different set in every mouse, and
in one case an average across two shanks:

| mouse | eegCh | nCh | spike groups | low-pass |
|---|---|---|---|---|
| raMCU1 | **[5 6 7 8]** — grp2 + grp3 | 15 | 4/3/4/4 | 450 |
| raMCU2 | [1 2 3 4] (ch3 is the weakest of all 15) | 15 | 4/3/4/4 | 450 |
| raMCU3 | [1 2 3 4] | 21 | 4×4 | **none** |
| raMCU4 | [7 8] | 12 | 1×12 | 450 |
| raMCU5 | [1 2 3 4] | 12 | 3×4 | 450 |

raMCU1 and raMCU2 share a recording date and probe (080621) and are exactly the
two mice whose discharges look atypical. Batch is confounded with phenotype.

Mean discharge amplitude also differs ~10× between mice (raMCU3 ~2000–2800 µV;
raMCU4/raMCU5 ~300–400 µV), so no absolute amplitude criterion transfers. The
one amplitude measure that does is **the ratio to that recording's own
ripples**: 3.3 / 4.4 / 4.7 in raMCU3/4/5, matching the ">5-fold over SPW-R"
of Maslarova et al. 2025. raMCU1's gated events reach only 2.4.

## Refuted: spatial extent as a criterion

Maslarova et al. report that an IED appears simultaneously on every hippocampal
channel while a SPW-R is confined to CA1, which suggests a shape-free ED test.
Measured here (`ed_laminar.m`) it **does not separate**: curated discharges and
the strongest gate-rejected candidates both hit every channel, with the same
amplitude-invariant spread (0.57–0.90 vs 0.67–0.92). Their criterion needs
their probe — 1024 channels at 30 µm pitch spanning all subfields. These
recordings give a nearly flat profile, so the detection channel cannot be
placed in a layer from the data, and the s.pyr-vs-dendritic morphology
distinction in their figures is not directly applicable here.

## Confirmed: firing suppression is the non-circular criterion

Every waveform measure describes the same trace the event was detected in, so
it can only agree with the detector. Spiking is a consequence, measured on
other electrodes (`ed_units.m`). Population rate 50–300 ms after the event,
over each event's own baseline 200–500 ms before:

| mouse | curated ED | gated | strong rejects | ripples | random |
|---|---|---|---|---|---|
| raMCU1 | — | **0.47** | 0.77 | 1.04 | 1.11 |
| raMCU2 | — | 0.82 | 0.98 | 1.00 | 0.97 |
| raMCU3 | **0.32** | 0.42 | 0.80 | 0.93 | 1.07 |
| raMCU4 | **0.62** | 0.75 | 1.00 | 0.80 | 0.94 |
| raMCU5 | **0.18** | 0.60 | 1.01 | 0.99 | 0.95 |

Curated discharges suppress firing to 18–62% of baseline; rejects, ripples and
random times do not move. raMCU1's gated events suppress as strongly as
raMCU3's curated ones — evidence that raMCU1 carries real discharges that do
not match the raMCU3/raMCU5 waveform. Caveats: n is 10–19 events, the traces
are noisy, and raMCU1 shows some pre-event dip that should not be there.

No population burst is visible at the event itself (raMCU5: zero spikes in
±20 ms), most likely because the large transient blocks spike detection.

## Open

- **`posZ` encodes the raMCU3/4/5 shape** and is the criterion most likely to
  be wrong. Polarity is only 66% positive among curated discharges, is
  layer-dependent, and the curated set came from the old detector's proposals.
- **raMCU4 runs high** (~80 through the gate, vs 15 and 45 in raMCU3/raMCU5).
- **lh100 has an artifact-heavy channel** — saturation steps of 10–25k that no
  physiological criterion should have to handle. Worth checking the recording.
- **raMCU1 and raMCU2 are uncurated**, so their counts are unvalidated.
- The gate rests on 47 discharges from 3 mice. It is a starting point per mouse,
  not a constant.
