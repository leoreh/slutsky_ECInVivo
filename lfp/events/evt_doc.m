% EVT_DOC  Event pipelines (ripples + ED) — shared architecture and file map.
%
% Two detection pipelines — sharp-wave ripples (lfp/ripples) and epileptiform
% discharges (lfp/ed) — run in parallel on a shared, event-agnostic layer
% (lfp/events, the evt_* functions) and hand a common curation-ready struct to
% the GUI (guiPath). One session runs end to end through a wrapper; every
% result is a <basename>.<var>.mat file. This file is the shared-design
% reference for both pipelines.
%
% The two pipelines are parallel by construction: they share every
% detection-agnostic step (spike prep, bout times, QA, matched controls, state
% labelling, maps, spike analysis, saving, plotting) and differ only where the
% modality genuinely differs — the detection model, the per-event features, and
% the signal loaded. Only detection is modality-specific: ripples z-score a
% narrowband envelope; ED z-score a moving baseline to flag sharp transients.
% Every side effect is gated (flgSave writes disk, flgPlot draws figures,
% flgCurate opens the GUI), so a pipeline runs headless for verification.
%
% # Entry points
%
% ripp_wrapper  Ripple orchestrator: chains the three ripple stages for the
%               batch case - ripp_detect (features) -> ripp_curate (the .accepted
%               gate) -> ripp_analyze (spikes/phase/maps on accepted). Each stage
%               is also a standalone function, so a session can stop after detect
%               for manual curation in ripp_curate's GUI, then run ripp_analyze.
% ed_wrapper    ED orchestrator: the same shape with two stages -
%               ed_detect -> ed_curate. There is no analyze stage, because the
%               ED question (how many discharges, in which states) is answered
%               by the curate stage's edStates. The per-event pass that follows
%               the bulk gate is guiPath, preset 'ed'.
%
% Re-run model (both): if <basename>.<name>.mat exists and flgForce is false,
% the stored result is loaded (the cheap re-curate path); otherwise the pipeline
% runs. Signal loading lives inside the detect branch, so re-curation touches no
% signals. A forced re-detection backs the previous file up first, so a curated
% mask is never lost silently.
%
% # Canonical event schema
%
% Both ripp (<basename>.ripp.mat) and ed (<basename>.ed.mat) carry the same
% curation contract plus per-modality metric columns:
%   .times     [N x 2]  event start / stop (s, absolute)
%   .peakTime  [N x 1]  event peak (s, absolute)
%   .state     [N x 1]  categorical vigilance state (<undefined> if none)
%   .accepted  [N x 1]  logical curation mask (seeded all-true post-QA)
%   .ctrlTimes [N x 2]  matched control intervals (analysis)
%   .info      struct   detection params, fs, basename, win, runtime
% Ripple metric columns: .amp .freq .freqEvent .freqPeak .peakProm .energy .dur
%                        .skew .emg .spkGain.
% ED metric columns:     .fastZ .posZ .amp .dur .emg.
% QA is now the same on both sides: MARK, never remove. The detect stage seeds
% .accepted all-true and evt_gate turns a filter spec (states + per-metric
% [lo hi] ranges) into the mask, which the curate stage applies headless (the
% automatic gate) or in a GUI (per-mouse). Nothing is dropped, so a rejected
% event keeps its row and can be put back. The shipped specs differ in
% strictness on purpose: the ripple default is permissive because its gate is
% chosen per mouse over thousands of events, while the ED default is active
% because an ED session ends in a per-event pass that is only tractable once
% the artifact bulk is gone.
%
% # Shared layer (lfp/events)
%
% Event-agnostic; both wrappers call these with modality-specific inputs.
% evt_files     Standard output file paths for a modality ('ripp' / 'ed').
% evt_spkPrep   Window-relative single-unit + pooled-MUA times, unit types.
% evt_boutTimes Window-relative bout times + valid-state + NREM baseline sets.
% evt_emgScore  Per-event EMG z vs a baseline window (a QA metric).
% evt_spkGain   Per-event MUA spike-gain z (a QA metric; ripples).
% evt_gate      The QA engine shared by both pipelines: a filter spec (states +
%               per-metric [lo hi] ranges) -> the .accepted mask. Permissive by
%               construction - an empty state list keeps every state and a NaN
%               metric passes, so an unavailable criterion relaxes rather than
%               rejects.
% evt_subset    Subset every per-event field by a logical mask (used by readers
%               that want the accepted rows only, e.g. mcu_tblVivo).
% evt_ctrlTimes Duration-matched control intervals from valid states.
% evt_states    Per-event vigilance state + the per-bout rate/density table.
% evt_maps      Per-event signal maps around each event peak.
% evt_spks      Spike entry point: per-unit stats + per-event population metrics
%               + 3D raster + per-unit PETH (orchestrates evt_spksParams,
%               evt_spkPeth, evt_pethNorm, evt_rankOrder).
% evt_saveSpks  Split the spike result into a light stats file + heavy raster.
% evt_plotSpks  Spike-modulation summary figure.
% evt2ns        Export events to a NeuroScope .evt file (start / peak / stop
%               marks, labelled by acceptance). Ripples wire it via flgNS.
% guiPath Manual curation; reads and writes the .accepted mask.
%
% # Modality-specific (deliberately not shared)
%
% Ripples: ripp_pickCh (detection channel), ripp_sigLoad (channel + EMG),
%   ripp_sigPrep (filter + envelope + z-score), ripp_times (threshold
%   candidates), ripp_params (frequency / energy / skew), ripp_detect (stage 1:
%   detect + features), evt_gate (the filter-spec -> .accepted engine),
%   ripp_curate (stage 2: apply the gate headless or in a GUI), ripp_analyze
%   (stage 3: spikes/phase/maps on accepted), spklfp_phase (spike-LFP).
% QA marks .accepted (ripples, never removed) or subsets (ED); mcu_tblVivo
%   filters .accepted for the per-event ripple table.
% ED: ed_methods (the met config), ed_sigLoad (sleep_sig eeg or the ripple
%   channel + EMG), ed_sigPrep (discharge band + a robust block-wise scale +
%   the supra-physiological band), ed_detect (stage 1: candidates + features),
%   ed_curate (stage 2: the gate, headless or in a bulk GUI), ed_tbl (counts and
%   rates per session x state, the pipeline's cross-session product).
% The channel loaders share ripp_pickCh (the ripple channel) and evt_loadCh
% (binary .lfp read + averaging + bit2uv autodetect).
%
% # Outputs  (<basename>.<var>.mat, written when flgSave = true)
%
%   Concern             Ripples             ED                Producer
%   Events (curation)   .ripp.mat           .ed.mat           wrapper
%   Per-bout states     .rippStates.mat     .edStates.mat     evt_states
%   Signal maps         .rippMaps.mat       .edMaps.mat       evt_maps
%   Spikes: stats+PETH  .rippSpks.mat       .edSpks.mat       evt_saveSpks
%   Spike raster (3D)   .rippSpkMaps.mat    .edSpkMaps.mat    evt_saveSpks
%   Phase coupling      .rippSpkLfp.mat     —                 spklfp_phase
%
% The light spikes file (.rippSpks / .edSpks) carries per-unit scalar stats +
% the per-unit PETH + tstamps; the manuscript loads it hot. The 3D raster is
% split into its own heavy file. A per-type population PETH is not precomputed;
% derive it from the 3D raster if a table ever needs it.
%
% # Robustness
%
% Both wrappers degrade gracefully when sleep states, spikes, units, or EMG are
% absent (real EA / RA sessions): missing spikes skip the spike analyses;
% missing EMG / NREM relax the corresponding QA criteria; 'nrem' z-scoring falls
% back to 'adaptive'.
%
% Terminology: a "map" or "raster" is the 3D per-event array [unit x event x
% bin]; a "PETH" is its per-unit average over events (2D [unit x bin]).
