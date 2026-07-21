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

## The channel (260721 rebuild)

Detection reads ONE auto-picked raw `.lfp` channel (`ed_pickCh`); the `sleep_sig`
eeg branch and `met.chMode` are gone. Every session in the MCU cohort has a
binary `.lfp` and a `session.mat`, so the fallback earned nothing. The EA cohort
has neither and is no longer runnable.

**Ground truth for the picker** (`ed_chSweep.m`): the real detector run on EVERY
neural channel of raMCU3/4/5, each scored by the AUC with which `fastZ`
separates the curated discharges. AUC varies hugely across channels of one
probe — 0.78 to 0.96 in raMCU3 — so the choice matters.

Two candidate criteria were refuted before the shipped one:

| criterion | raMCU3 | raMCU4 | raMCU5 | verdict |
|---|---|---|---|---|
| `p99.9 / background` | rank 3 | rank 6 | rank 6 | normalising rewards the quietest, usually marginal, channel |
| kurtosis | rank 3 | rank 6 | rank 6 | same failure |
| fewest candidates | AUC 0.810 | — | 0.978 | corr(−xRate, AUC) = **−0.69** in raMCU3: few crossings means the discharge itself barely crosses |

Shipped: **drop channels whose crossing rate is a Tukey outlier, then take the
largest transients** (`prctile(|filt|, 99.9)`). Step 1 exists because raMCU5's
ch12 produced 13457 candidates where every other channel gave 4760–6495 — a
broken electrode, not a better one.

| mouse | pick | its AUC | best possible | worst |
|---|---|---|---|---|
| raMCU3 | ch 14 | 0.957 | 0.960 | 0.780 |
| raMCU4 | ch 7 | **0.986** | 0.986 | 0.880 |
| raMCU5 | ch 3 | 0.978 | 0.983 | 0.951 |

**Averaging channels was tested and refuted.** A discharge appears on every
channel at once, so averaging looked like a matched filter for it. Measured
(`ed_sigSweep.m`), the mean of all neural channels is *worse* — AUC 0.862 vs
0.957 in raMCU3, 0.892 vs 0.986 in raMCU4. One channel wins.

## Shipped defaults

Detection: `passband [60 150]`, `thr 8`, `limDur [4 200 40] ms`, polarity-blind,
on one auto-picked channel.

Noise filter (`met.qa`): `fastZ ≥ 5`, `isoZ ≥ 5` — loosened by Leore from the
15 / 20 the sweep chose, which is what forced the cluster count to scale (below).
**`posZ` was dropped from the gate.** It required the event to rise above its
pre-event baseline, which encodes the shape of the discharges in three mice;
polarity is layer-dependent (Maslarova et al. report a sharp negative spike in
the dendritic layers and a positive slow wave in the pyramidal layer for the
*same* event), and only 66% of the curated discharges were positive-going. It is
still measured and reported.

Clustering (`met.clust`): `win [-0.05 0.05]`, `nPC 6`, `nClust []` (auto:
0.65·√n), plus the scalar measures.

`emg` is gone entirely — nothing read it once the shape criteria were applied.

## Cluster curation: what it does and does not do

`ed_clustSweep.m` sweeps window × normalisation × features × count on the three
curated mice, scored by the events you must review to find 80% of the
discharges, taking clusters in order of discharge density.

- **Scalar features are decisive.** Eleven of the top twelve configurations use
  waveform components *plus* `fastZ isoZ posZ amp dur`; shape alone is worse
  everywhere.
- **A longer window wins**, against the paper's 10–50 ms: ±50 and ±100 ms beat
  ±15 and ±25 by roughly threefold in review load. A discharge and a sharp wave
  differ most in the DECAY, which a ±15 ms window cannot see.
- **Unit-peak normalisation wins**, against the paper's absolute amplitudes. The
  pool spans two orders of magnitude and the largest events take the principal
  components with them. Scale returns through the ranked scalars.
- **A fixed count beats BIC**, which settled on ~7 where 12 measured better on
  every mouse. BIC maximises likelihood, which is not the objective.

**Honest limits.** No single cluster holds all the discharges — they spread over
two to four, so several boxes get ticked. Best-cluster purity tops out around
70–88%, and the 80% figure is 80%, not 100%. This makes curation a handful of
decisions instead of hundreds; it does not make it exact. A mouse whose result
matters still deserves the per-event pass in `guiPath`.

## The curation GUI is one pivotable view, not two fixed ones

First attempt drew its own per-cluster tiles (median + IQR) with a peri-event
MUA row, and dropped the threshold knobs. That was wrong three ways, and Leore
caught all three:

- **The knobs are needed.** They set the pool, and the pool is what gets
  clustered. They are back, and take effect on `Re-cluster` — the label updates
  live so you can see the pool size before paying for the fit.
- **The per-state view was lost.** It is back, and now it is the *same* view:
  `guiTbl_xy` over a table of every detected event carrying `lfp`, `cluster`,
  `state`, `status`, so "Plot By (Tiles)" pivots between per-cluster and
  per-state and "Group By (Colors)" overlays the other.
- **Hand-drawn tiles were a reinvention.** `guiTbl_xy` already does tiles,
  grouping, per-category show/hide, and a Dispersion + Median trace — which is
  the robust central waveform the hand-drawn version existed to provide.

The MUA row was removed rather than kept as a Y variable: a dozen events per
cluster is too few to read, and `guiTbl_xy` takes one x-axis, which a ±500 ms
suppression window cannot share with a ±100 ms waveform.

**Accept and view are separate controls.** Mixing them meant you could not
inspect what you rejected without changing what you kept.

- ACCEPT is checkboxes: which clusters, and which states. An event is accepted
  when its cluster AND its state are ticked. Both lists start **ticked** —
  curation is rejection, because the eye is much better at spotting the two or
  three tiles that are obviously not discharges than at confirming the ten that
  are. States are there to drop a stretch of recording wholesale — movement
  artifact in WAKE — without touching the shape decision.
- SHOW is a dropdown: `both | accepted | removed`. Rows drawn, nothing else.

**Re-cluster fits only the ACCEPTED events** (intersected with the thresholds),
so rejecting then re-clustering is a refinement loop. The accepted set does not
change when you press it. `Reset to filter` undoes the narrowing. `ed_clust`
also seeds its RNG (restored on exit), so the same pool and count give the same
partition — without that, `fitgmdist`'s random starts reshuffled the very
groups the user had just judged.

**Reopening resumes** the saved labels, ticks, state selection and thresholds
rather than re-fitting. The saved labels are the partition that was judged; a
re-fit would return a different one and leave the ticks pointing at groups
nobody looked at. A restore that does not line up with the event list is
refused, so a re-detection falls through to a fresh clustering.

## The cluster count has to scale with the pool

The filter was later loosened (`fastZ ≥ 5`, `isoZ ≥ 5`), which took a pool from
75–490 to 2700–8500 — and a fixed 12 groups over 8500 events leaves every group
a mixture. Re-measured at the loose filter (`ed_winSweep.m`), best-cluster
purity against the curated discharges:

| mouse | pool | k=12 | k=20 | k=30 | k=45 | k=60 |
|---|---|---|---|---|---|---|
| raMCU3 | 8460 | 3% | 6% | 7% | 14% | **80%** |
| raMCU4 | 2698 | 11% | 40% | **73%** | 39% | 35% |
| raMCU5 | 4158 | 7% | 13% | 50% | **71%** | 50% |

So `nClust` empty now means **0.65·√n**, which picks 60 / 34 / 42 for those
pools and 11–14 for the old strict ones, where a fixed 12 had measured best.

This is what makes the refinement loop work, and it is fragile in one specific
way. The GUI used to write the resolved count back into its `clusters (0 =
auto)` box, which reads as helpful — it shows what auto chose — but it turns
auto into a fixed number after the first fit. Measured on raMCU4: open at 2698
events → 34 groups; reject two clusters and WAKE → 512 accepted; `Re-cluster`
then split those 512 into **34** groups again instead of 15. The box now holds
what the user *asked* for and the label reports what came *back*; regression
test `test_curateAutoCountFollowsPool`.

Caveat worth keeping in view: purity is estimated from 10–19 curated events per
mouse, so individual cells above are noisy. What is robust is the direction —
k=12 is too small at these pool sizes on every mouse.

**The window is not the problem.** Swept at the loose filter, best review load
per window: ±10 ms 4.8%, ±20 ms 6.8%, ±30 ms 5.7%, **±50 ms 3.8%**, ±100 ms
4.6%. No penalty for the wider window; ±50 still wins.

Dropping `.amp` from the scalar set ('shape') beat keeping it (3.8% vs 4.5%),
and dropping the scalars entirely was far worse (13.7%) — confirming again that
they carry real information.

Cost, measured on the biggest session (lh100, 9531 events): 6.5 s to build the
view once, **1.1 s per interaction** thereafter, because the widget is fed rows
through `setDataFcn` instead of being rebuilt. That is also why a re-cluster
keeps whatever pivot the user set.

## End-to-end, whole cohort (`ed_verify.m`, ~8 s per session)

| session | grp | ch | candidates | pool | recall | review |
|---|---|---|---|---|---|---|
| lh100 | WT | 9 | 9531 | **16** | — | — |
| lh107 | WT | 14 | 355 | 9 | — | — |
| lh119 | WT | 2 | 1673 | 8 | — | — |
| lh122 | WT | 2 | 1432 | **239** | — | — |
| lh126 | WT | 15 | 1283 | 4 | — | — |
| lh132 | MCU-KO | 2 | 3388 | 12 | — | — |
| lh133 | MCU-KO | 5 | 1014 | 4 | — | — |
| lh134 | MCU-KO | 8 | 46 | 6 | — | — |
| lh136 | MCU-KO | 15 | 1820 | 5 | — | — |
| lh137 | MCU-KO | 15 | 5706 | 12 | — | — |
| lh140 | MCU-KO | 2 | 10931 | 32 | — | — |
| raMCU1 | CAG | 6 | 1219 | 25 | — | — |
| raMCU2 | CAG | 5 | 4917 | 147 | — | — |
| raMCU3 | CAG | 14 | 8533 | 292 | 80% | **18** |
| raMCU4 | CAG | 7 | 2984 | 75 | 89% | **20** |
| raMCU5 | CAG | 3 | 4760 | 490 | 94% | **18** |

Review load is ~20 events per mouse regardless of pool size, which is the point.
Controls land in the single digits to low tens and are rejected wholesale.

The channel picker also fixed lh100 without being asked: it used to yield 236
through the filter on the artifact-heavy channel, and yields 16 on the one it
picks now.

## Architecture

`ed_detect` → `ed_curate` → `ed_tbl`, with `guiPath` as an optional per-event
pass. Detect writes nothing and seeds `.accepted` all-true; the filter and the
clustering are a separate stage over the saved struct; curation MARKS and never
removes; a save backs the file up first.

| file | role |
|---|---|
| `ed_methods` | the met config: detection, `.qa` noise filter, `.clust` |
| `ed_pickCh` | the detection channel (outlier reject, then largest transients) |
| `ed_sigLoad` | one raw `.lfp` channel, windowed |
| `ed_detect` | stage 1: signal prep, candidates, every per-event feature |
| `ed_params` | the features |
| `ed_clust` | PCA + GMM over waveform shape; pure, no I/O |
| `ed_curate` | stage 2: filter, cluster, accept whole types |
| `ed_wrapper` | the thin chain |
| `ed_tbl` | counts and rates per session × state |

`tests/test_edCurate.m` covers `ed_clust` and the GUI on a synthetic session.

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
lines; changed no result versus one recording-wide scale) and `.pol` (now the
sign of `.amp`). Also dropped as unrequested: `flgNS`, `flgDetectOnly`,
`mapDur`, `flgAll`. The state filter was removed from the *detection* gate on
the same grounds — ED asks how discharges distribute over states, so
restricting them there is a question for `ed_tbl` — but came back in curation
as an acceptance criterion, which is a different job: dropping a stretch of
recording, not defining an event.

A second pass over `ed_curate` (260721) cut it from 620 to 530 lines with no
behaviour change: `buildChecks` merged into the adopt step it was the tail of;
`selectedClusters`/`selectedStates`/`scalarFeat` folded into their one or two
call sites; the `selRestore` field, whose lifetime was a single call, became an
argument; guards on `chkState` that could never fire removed; a six-entry
same-day HISTORY collapsed into the one design it converged on. The auto-count
defect above was found by that pass, not by the tests.

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

- **lh122 yields a pool of 239** where every other control gives 4–32. Either
  that channel is bad or the mouse has something; look before trusting its
  count.
- **raMCU1 and raMCU2 are uncurated.** raMCU1's pool is 25 and raMCU2's is 147,
  and their discharges are the ones that do NOT match the raMCU3/4/5 shape.
  They are the reason the shape criterion was removed, and they are still the
  test of whether that was enough.
- **The curated set is 47 discharges from 3 mice**, all proposed by the OLD
  detector, so recall against it is recall against a biased sample. It is the
  only ground truth there is.
- **Clustering does not isolate a single discharge cluster.** Best-cluster
  recall is 39–70%; the discharges spread over two to four clusters. The review
  load is ~20 events per mouse, not zero.
- The suppression measure (`ed_units.m`) is evidence, deliberately kept out of
  the pipeline. It is not in `ed_curate` at all — the peri-event MUA panel was
  dropped (too few events per cluster to read, and it cannot share an x-axis
  with the waveform), so nothing in the pipeline reads or gates on it.
