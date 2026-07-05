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

| Concern           | Ripples              | ED                 | Producer      |
|-------------------|----------------------|--------------------|---------------|
| Events (curation) | `.ripp.mat`          | `.ed.mat`          | wrapper       |
| LFP maps          | `.rippMaps.mat`      | `.edMaps.mat`      | `evt_maps`    |
| Spike modulation  | `.rippSpks.mat`      | `.edSpks.mat`      | `evt_spks`    |
| Spike PETH        | `.rippPeth.mat`      | `.edPeth.mat`      | `evt_spkPeth` |
| Per-bout states   | `.rippStates.mat`    | `.edStates.mat`    | `evt_states`  |
| Phase coupling    | `.rippSpkLfp.mat`    | —                  | `spklfp_phase`|

## Shared layer (`lfp/events/`)

Event-agnostic; both wrappers call these with modality-specific inputs:
`evt_states`, `evt_ctrlTimes`, `evt_maps`, `evt_spks` (`winFxd` param),
`evt_spkPeth`, `evt_rankOrder`, `evt_rate`, `evt_plotSpks`. `evt_states` /
`evt_plotSpks` take a `name` param that sets the saved filename/variable/figure
(`'ripp'` or `'ed'`), preserving downstream readers (e.g. `mcu_tblVivo` expects
the variable `rippStates`).

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

- `ripp_qa` and `peth_norm` remain embedded in `ripp_wrapper` (ED uses its own
  `ed_reject_emg` QA and `evt_plotSpks` self-normalizes, so a shared `evt_qa` /
  `evt_pethNorm` was not required). Extract them if a second consumer appears.
- Some internal variable names inside the moved `evt_*` functions still read
  `ripp*` (cosmetic; the function names and interfaces are generic).
- ED detection quality is unvalidated; the new ED maps/PETH double as a
  detection read-out to inform a future tuning pass.
- Not migrated by request: `lfp/+IED`, `reduct_displayer` (owned elsewhere).
