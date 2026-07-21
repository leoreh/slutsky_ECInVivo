# Ripple curation rebuilt on waveform clusters (260722)

## The problem

`ripp_curate` was three threshold knobs (`emg`, `spkGain`, `peakProm`) over a
kept-versus-removed mean waveform. Thresholds cannot see what an event looks
like. A step artifact with quiet muscle and a bystander burst passes every one
of them, and averaging a set that holds ripples, steps and spike-bleed
transients gives a curve that is none of them.

The ED pipeline had already solved this (260721): group the pool by waveform
shape, then accept or reject whole shapes. A ripple session needs the same tool
and needed it more — 30k–134k detected events per session against an ED pool of
3k–8k.

## What was done

`ed_curate`'s GUI was lifted into `lfp/events/evt_curate.m` and `ed_clust` into
`lfp/events/evt_clust.m`, following the `ripp_gate → evt_gate` precedent.
`ed_curate` and `ripp_curate` are now loaders over one engine: each finds its
own event struct and waveforms, hands them over with a config, and the shared
code does the pool, the clustering, the accept/reject model and the save.

Three things became generic in the move.

- **The metric knobs are built from `met.qa.ranges`.** One edit box per finite
  bound, labelled with its metric. The GUI never names a metric, so adding one
  to `ripp_methods` puts it in the GUI with no other change. ED gets `fastZ >=`
  and `isoZ >=`; ripples get `emg <=` and `spkGain >=`.
- **The waveforms are an input.** The engine does not know how to find a maps
  file. That is what let the two loaders keep their different rules — ED errors
  on a missing `edMaps`, ripples fall back to rebuilding the detection signal.
- **The state scope lives in the saved `info.qa`.** `info.clustStates` is still
  read for files written before this and is never written again, so the saved
  spec is one complete, replayable `evt_gate` input.

## What is new for ripples

`met.clust` in `ripp_methods`:

| field | value | why |
|---|---|---|
| `win` | `[-0.03 0.03]` | a ripple separates from a transient in the oscillation; a discharge separates from a sharp wave in the decay, which needs ED's ±50 ms |
| `nClust` | `20` | the sqrt rule is built for a RARE shape; a ripple pool's contaminant is a whole population, and `0.65*sqrt(30000)` = 110 tiles to read |
| `nFit` | `8000` | the mixture is estimated on 8000 events, the rest projected and assigned |
| `nView` | `300` | per cluster, for the view only |
| `scalar` | `peakProm freqPeak amp dur emg spkGain` | shape descriptors already computed |

`.qa` stopped being the answer and became the noise filter that sizes the pool.
The two criteria left in it ask only what no ripple can fail — was the muscle
quiet, did units fire. Everything about what a ripple *looks* like is now the
clustering's business.

## Cost, measured

Opening the GUI on real sessions (nothing saved):

| session | detected | pool | open | re-cluster | `UserData` |
|---|---|---|---|---|---|
| lh100 | 107,154 | 14,208 | 13 s | 9.1 s | 145 MB |
| lh137 | 32,871 | 28,543 | 11 s | 8.9 s | 43 MB |
| raMCU5 | 133,678 | 8,301 | 13 s | 6.7 s | 180 MB |

For scale, the ED GUI was measured at 6.5 s to open and 1.1 s per interaction on
its biggest session (9,531 events).

Open is dominated by the `rippMaps` read — the file stores six mapped signals
and `load` cannot fetch two of them selectively. Waveforms are kept `single`,
which halves what the figure holds; `evt_clust` casts its own window.

Neither cap touches the mask. Every event is labelled and every event is saved;
the counts on the cluster checkboxes are the true ones, and the view says so
when it has been thinned.

## Does it separate anything?

raMCU5, 20 clusters over its 8,301-event pool. Clusters 9, 17, 19 and 20 (n =
407, 237, 158, 155) carry median EMG scores of 0.75, 0.77, 0.30 and 0.64 against
~0.10 everywhere else, and the lowest spectral prominences in the session (7.6,
4.6, 5.5, 5.7 against 10–28). That is the movement-artifact population, and it
is four unticks.

lh100 and lh137 do not separate that way — their clusters are homogeneous in
prominence, gain and EMG, which is the same answer blind clustering gave when
the question was first asked of the manuscript cohort. The tool finds a
contaminant where there is one.

## Defects found by review, and fixed

A four-lens adversarial review over the diff produced 47 claims; six survived
refutation and were fixed, each with a test:

| what | why it mattered |
|---|---|
| headless nulled `clustId` / `clustSel` | one `ripp_wrapper` run would have destroyed a session's manual curation, silently. Headless makes no shape judgement, so it now leaves the partition alone |
| `restoreSaved` re-ticked everything when `clustSel` was empty | rejecting EVERY cluster is an answer; it was being read as "nothing was saved" and inverted into accept-all on reopen |
| `0 = auto` in the cluster box never reached the rule | `k < 2` fell back to `met.clust.nClust` (20 / 12), not to `[]`. The label and four comments described behaviour the code did not have |
| `nFit` below the fit floor became a 20-event fit | `max(MINEV, ...)` turned a misconfigured `0` into a partition estimated from 20 points that looks fine and means nothing. Now read as no cap |
| `info.qa` could claim a filter the mask disobeys | a knob typed into after a fit only previews, but Save read the live box. Save now records the ranges the fit used |
| edge events kept an unclearable "unsorted" hint | an all-NaN waveform (peak within the map window of the recording edge) re-entered the pool on every Re-cluster. Now recorded as rejected |

Also fixed in passing: the GUI's verbose line printed the SEED mask as if it
were a result, and `buildStates` did not pass `basename` (the old ripple version
never did, so every session in a flat folder overwrote one
`<folder>.rippStates.mat` — `ed_curate` had it right).

## Verification

`tests/test_edCurate.m` (32 tests) and `tests/test_rippCurate.m` (28 tests) both
pass. The ED suite is the regression harness for the extraction — it was updated
only where a field moved (`st.ed` → `st.evt`; the `fastZ` box is now looked up by
metric rather than by the name `edFast`). `test_clustNFitOffIsUnchanged` pins the
claim that an empty `nFit` leaves the ED partition bit-identical.

## Not done

- `rippMaps` is read whole. Reading two of its six fields would need `h5read`
  against MATLAB's `-v7.3` struct layout; 18 s to open did not justify it.
- `ed_methods('default').clust.nClust` is `12` while its own comment documents
  the auto rule. Left alone: it is a live tuning value in a pipeline that was
  just curated with it, and changing it is an ED decision, not a ripple one.
