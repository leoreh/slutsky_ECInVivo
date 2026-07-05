# Event pipelines: ripples + epileptiform discharges (ED)

Two detection pipelines — ripples (`lfp/ripples/`) and epileptiform discharges
(`lfp/ed/`) — run in parallel on a shared, event-agnostic layer (`lfp/events/`)
and produce a common curation-ready schema for `graphics/gui_curate.m`.

> **Path setup:** add `lfp/events/` to the MATLAB path (Set Path → Add Folder,
> or your startup) alongside `lfp/ripples/` and `lfp/ed/`.

## Canonical event schema

Both `<basename>.ripp.mat` (`ripp`) and `<basename>.ed.mat` (`ed`) carry the
same curation contract plus per-modality metric columns:

| Field       | Shape    | Meaning                                             |
|-------------|----------|-----------------------------------------------------|
| `.times`    | [N x 2]  | event start / stop (s, absolute)                    |
| `.peakTime` | [N x 1]  | event peak (s, absolute)                            |
| `.state`    | [N x 1]  | categorical vigilance state (`<undefined>` if none) |
| `.accepted` | [N x 1]  | logical curation mask (seeded all-true post-QA)     |
| `.ctrlTimes`| [N x 2]  | matched control intervals (analysis)                |
| `.info`     | struct   | detection params, fs, basename, win, runtime        |

Automatic QA (EMG / amplitude / duration / spike-gain / state) runs as a
**filter** at detection — only surviving candidates are saved, `.accepted`
starts all-true, and curation in `gui_curate` flips it. No `.idxQA` breakdown is
stored. Ripple metric columns: `.amp .freq .freqEvent .energy .dur .skew .emg
.spkGain`. ED metric columns: `.amp .ampZ .dur .width10 .emg`.

## File layout (parallel by construction)

Split **by dependency**: maps need only the LFP (always producible), spikes need
sorted units (skipped when absent). The old separate PETH file is gone — the 3D
peri-event spike maps now live inside the spikes file.

| Concern            | Ripples              | ED                 | Producer          |
|--------------------|----------------------|--------------------|-------------------|
| Events (curation)  | `.ripp.mat`          | `.ed.mat`          | wrapper           |
| LFP maps           | `.rippMaps.mat`      | `.edMaps.mat`      | `evt_maps`        |
| Spikes: stats+PETH | `.rippSpks.mat`      | `.edSpks.mat`      | `evt_spks`        |
| Spike raster (3D)  | `.rippSpkMaps.mat`   | (in `.edSpks`)     | `evt_spks`        |
| Per-bout states    | `.rippStates.mat`    | `.edStates.mat`    | `evt_states`      |
| Phase coupling     | `.rippSpkLfp.mat`    | —                  | `spklfp_phase`    |

`rippSpks`/`edSpks` carry per-unit scalar stats + the per-unit PETH (`.peth`, 2D
per-unit average) + `.tstamps`. The 3D raster (`.su`/`.mu`, each `.evt/.ctrl`) is
split into its own light-vs-heavy file — `rippSpkMaps` for ripples (so the
manuscript loader stays light); ED still carries it inline in `.edSpks.maps`
(pending the ED pass). Per-event population PETHs (RS/FS/MU) are **computed on
demand** from the 3D via `evt_pethPop` — never precomputed or saved.

## Shared layer (`lfp/events/`)

Event-agnostic; both wrappers call these with modality-specific inputs.

**Spikes.** Each wrapper's spike section is one prep call + one analysis call:
- `evt_spkPrep(v, win, sigDur, fsSpk)` → window-relative single-unit times, one
  pooled MUA vector, and unit types (empty when absent; guards the fields of the
  loader struct `v` internally).
- `evt_spks(spkTimes, muTimes, evtTimes, ctrlTimes, peakTime, ...)` → one struct:
  per-unit stats + per-event population metrics (`.events`), 3D raster maps
  (`.maps.su/.mu`), and the per-unit PETH (`.peth`/`.tstamps`). The spike entry
  point; orchestrates `evt_spksParams` (per-unit scalar stats, `winFxd` param) +
  `evt_spkPeth` + `evt_pethNorm`.
- `evt_pethNorm` (smooth + z-score against control; kernel from the PETH time
  base) and `evt_pethPop` (per-event population PETH from the 3D) are the reused
  reduction helpers. `evt_viewSpks` launches the shared interactive viewers.

**Events / maps / rate.** `evt_states`, `evt_ctrlTimes`, `evt_maps`, `evt_rate`,
`evt_rankOrder`, `evt_plotSpks`. `evt_states` / `evt_plotSpks` take a `name` param
that sets the saved filename/variable/figure (`'ripp'` or `'ed'`), preserving
downstream readers (e.g. `mcu_tblVivo` expects the variable `rippStates`), and a
`lbl` param that labels titles/legends (`'Ripple'` / `'ED'`).

## Modality-specific (deliberately not shared)

- **Detection** — different signal models: ripples use a narrowband envelope
  (`ripp_sigPrep` + `ripp_times`); ED uses a moving-z sharp-transient detector
  (`ed_detect`). Features differ too (`ripp_params` frequency/energy/skew from
  the analytic signal vs `ed_params` half-/10%-amplitude widths).
- **Signal loading** — `ripp_sigLoad` (embedded in `ripp_wrapper`) vs
  `ed_sigLoad`.
- **Orchestrators** — `ripp_wrapper`, `ed_wrapper` (thin).

## Robustness

Both wrappers degrade gracefully when sleep states, spikes, units, or EMG are
absent (real EA / RA sessions): missing spikes skip the spike analyses; missing
EMG / NREM relax the corresponding QA criteria; `'nrem'` z-scoring falls back to
`'adaptive'`.

## Deferred / notes

- `ripp_qa` stays embedded in `ripp_wrapper` (ED uses its own `ed_reject_emg`);
  extract a shared `evt_qa` only if a third consumer appears.
- ED detection quality is unvalidated; the ED maps/PETH double as a detection
  read-out to inform a future tuning pass.
- Not migrated by request: `lfp/+IED`, `reduct_displayer` (owned elsewhere).
